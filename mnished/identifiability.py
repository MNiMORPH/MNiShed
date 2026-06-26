"""
mnished.identifiability
~~~~~~~~~~~~~~~~~~~~~~~~~
Post-fit parameter-identifiability diagnostics for MNiShed calibration.

This module answers a single question about a *calibrated* model: **how
well does the data actually constrain each parameter — and which
*combinations* of parameters does it fail to constrain?**  It is a
diagnostic on the local shape of the calibration objective around an
optimum, not an optimiser or an uncertainty-quantification framework
(for the rigorous, global, posterior-based answer, drive Dakota's
``bayes_calibration`` from the same evaluator — see the module notes).

The core operates on two abstractions, both deliberately generic so the
machinery is independent of any particular parameter set:

* :class:`ParameterSet` — the free parameters at the optimum, each with
  its calibration-coordinate value and bounds.  *Calibration coordinate*
  means the coordinate the optimiser varies: log10(value) for a
  ``log__`` parameter, the natural value otherwise.  This matches the
  Dakota variables block exactly.
* an **objective callable** ``f(theta: dict) -> float`` mapping a
  ``{name: calibration_value}`` dict to a score (higher is better, e.g.
  KGE).  The caller supplies this closure; for an MNiShed calibration it
  wraps :func:`~mnished.run_and_score`, applying the same ``10**`` and
  argument-assembly its Dakota ``driver.py`` does.

What this first layer provides
------------------------------
* :func:`profile` / :func:`profile_all` — sweep one parameter across its
  bounds with the rest held at the optimum (the generalised "f_route
  sweep"); returns a :class:`Profile` exposing the Δ-flat half-width, a
  local-curvature width, and a bound-pegging flag.
* :class:`IdentifiabilityReport` — collects the per-parameter profiles,
  carries the calibration AIC, exposes a text :meth:`summary`, and holds
  an (initially empty) ``predicted`` slot where an *a-priori* structural
  prediction (diffusion-length / drainage-topology) will dock once those
  forward-measurement tools exist.

Single-parameter profiles see per-parameter flatness but **cannot see a
degenerate combination** (e.g. an ``f_route`` ↔ land-τ ridge): holding
one fixed while sweeping the other breaks their trade-off.  Naming those
combinations needs the curvature-eigenspectrum and 2-D ridge layers,
which build on this one.
"""

import numpy as np

__all__ = [
    "Parameter",
    "ParameterSet",
    "Profile",
    "Spectrum",
    "Ridge2D",
    "IdentifiabilityReport",
    "profile",
    "profile_all",
    "eigenspectrum",
    "ridge",
]


class _CachedObjective:
    """Memoised wrapper so Hessian corner-points are not re-evaluated.

    Keys on the rounded calibration-coordinate vector; the many repeated
    centre-point and shared-axis evaluations of a finite-difference
    Hessian then cost one model run apiece.
    """

    def __init__(self, objective, names, ndigits=10):
        self._f = objective
        self._names = names
        self._nd = ndigits
        self._cache = {}
        self.n_calls = 0

    def __call__(self, theta):
        key = tuple(round(float(theta[n]), self._nd) for n in self._names)
        if key not in self._cache:
            try:
                val = self._f(theta)
                val = float(val) if np.isfinite(val) else np.nan
            except Exception:
                val = np.nan
            self._cache[key] = val
            self.n_calls += 1
        return self._cache[key]


class Parameter:
    """A single free parameter at the calibration optimum.

    All of ``value``, ``lower`` and ``upper`` are in the *calibration
    coordinate* (log10 for a ``log__`` parameter, natural otherwise), so
    that bounds, sweeps and widths are all expressed in the coordinate
    the optimiser actually sees.  ``log`` is carried for labelling only;
    the objective closure is responsible for any ``10**`` back-transform.

    Parameters
    ----------
    name : str
        Parameter name, matching the key the objective closure expects.
    value : float
        Calibrated value (the optimum) in calibration coordinate.
    lower, upper : float
        Calibration-coordinate bounds.
    log : bool, optional
        True if ``name`` is a ``log__`` parameter (display only).
    description : str, optional
        Human-readable description for the summary.
    """

    def __init__(self, name, value, lower, upper, log=None, description=""):
        if upper <= lower:
            raise ValueError(
                f"Parameter '{name}': upper ({upper}) must exceed lower "
                f"({lower})."
            )
        self.name        = name
        self.value       = float(value)
        self.lower       = float(lower)
        self.upper       = float(upper)
        # Infer the log flag from the conventional name prefix if unset.
        self.log         = name.startswith("log__") if log is None else bool(log)
        self.description = description

    @property
    def range(self):
        """Width of the calibration-coordinate bound interval."""
        return self.upper - self.lower

    def normalize(self, x):
        """Map a calibration-coordinate value to a [0, 1] fraction of range."""
        return (np.asarray(x, dtype=float) - self.lower) / self.range

    def __repr__(self):
        return (f"Parameter({self.name!r}, value={self.value:.4g}, "
                f"bounds=[{self.lower:.4g}, {self.upper:.4g}], log={self.log})")


class ParameterSet:
    """An ordered collection of :class:`Parameter` at one optimum."""

    def __init__(self, parameters):
        self.parameters = list(parameters)
        names = [p.name for p in self.parameters]
        if len(names) != len(set(names)):
            raise ValueError("Duplicate parameter names in ParameterSet.")
        self._by_name = {p.name: p for p in self.parameters}

    def __len__(self):
        return len(self.parameters)

    def __iter__(self):
        return iter(self.parameters)

    def __getitem__(self, key):
        if isinstance(key, str):
            return self._by_name[key]
        return self.parameters[key]

    @property
    def names(self):
        return [p.name for p in self.parameters]

    def optimum(self):
        """The optimum as a ``{name: value}`` dict in calibration coordinate."""
        return {p.name: p.value for p in self.parameters}

    @classmethod
    def from_params_yml(cls, param_cfg, optimum=None):
        """Build from a ``params.yml`` ``parameters:`` block.

        Only ``active: true`` parameters become free :class:`Parameter`
        objects (fixed parameters are not identifiability questions).
        Each entry's ``lower``/``upper`` are taken as the calibration
        bounds; the value is taken from ``optimum[name]`` if given, else
        the entry's ``initial``.

        Parameters
        ----------
        param_cfg : dict
            The ``parameters`` mapping from a parsed ``params.yml``.
        optimum : dict, optional
            ``{name: calibration_value}`` of the calibrated optimum.
            Falls back to each parameter's ``initial`` when absent.
        """
        optimum = optimum or {}
        params = []
        for name, spec in param_cfg.items():
            if not spec.get("active", True):
                continue
            params.append(Parameter(
                name=name,
                value=optimum.get(name, spec["initial"]),
                lower=spec["lower"],
                upper=spec["upper"],
                description=spec.get("description", ""),
            ))
        return cls(params)


class Profile:
    """A 1-D objective profile for one parameter, others held at optimum.

    Attributes
    ----------
    name : str
        Parameter that was swept.
    x : ndarray
        Calibration-coordinate sweep values.
    score : ndarray
        Objective (higher better) at each sweep value.
    x_opt, score_opt : float
        The held optimum value and the score there.
    parameter : Parameter
        The swept parameter (for bounds / labelling).
    """

    def __init__(self, parameter, x, score, x_opt, score_opt):
        self.parameter = parameter
        self.name      = parameter.name
        self.x         = np.asarray(x, dtype=float)
        self.score     = np.asarray(score, dtype=float)
        self.x_opt     = float(x_opt)
        self.score_opt = float(score_opt)

    # -- identifiability summaries ---------------------------------------

    def half_width(self, delta=0.01):
        """Δ-flat half-width: how far the parameter can move, in fraction
        of its bound range, before the score degrades by ``delta``.

        Measured as the half-extent of the contiguous interval around the
        optimum over which ``score >= score_opt - delta``, expressed as a
        fraction of the bound range so it is comparable across parameters.
        A value near (or above) the full range means the data barely
        constrains the parameter; a small value means it is well pinned.

        Returns ``np.nan`` if the optimum sample cannot be located.
        """
        s = self.score
        thresh = self.score_opt - delta
        i_opt = int(np.argmin(np.abs(self.x - self.x_opt)))
        if not np.isfinite(s[i_opt]):
            return np.nan
        # walk outward while the profile stays within delta of the peak
        lo = i_opt
        while lo - 1 >= 0 and np.isfinite(s[lo - 1]) and s[lo - 1] >= thresh:
            lo -= 1
        hi = i_opt
        n = len(s)
        while hi + 1 < n and np.isfinite(s[hi + 1]) and s[hi + 1] >= thresh:
            hi += 1
        span = self.x[hi] - self.x[lo]
        # half-width as a fraction of the bound range
        return 0.5 * span / self.parameter.range

    def curvature(self):
        """Local quadratic curvature of the score at the optimum.

        Fits ``score ≈ a + b·u + c·u²`` in the normalized coordinate
        ``u = (x - x_opt)/range`` over the samples nearest the optimum,
        and returns ``c`` (negative at a well-identified maximum).  The
        magnitude is, up to the unknown error scale, the observed
        information in that parameter — the larger ``|c|``, the sharper
        the constraint.  Returns ``np.nan`` if too few finite samples.
        """
        u = (self.x - self.x_opt) / self.parameter.range
        s = self.score
        m = np.isfinite(s)
        if m.sum() < 3:
            return np.nan
        # central window: the 5 samples closest to the optimum (or all)
        order = np.argsort(np.abs(u[m]))
        uu, ss = u[m][order], s[m][order]
        k = min(len(uu), 5)
        c2, c1, c0 = np.polyfit(uu[:k], ss[:k], 2)
        return float(c2)

    def at_bound(self, tol=1e-3, delta=0.01):
        """Whether the optimum sits at (or wants to run past) a bound.

        True if either (a) the optimum is within ``tol`` of a bound — the
        signature of a parameter the optimiser has pushed against its
        prior, the way a stretched matrix-τ pegs at its upper bound — or
        (b) the profile's best sample is at an endpoint *and* beats the
        held optimum by more than ``delta`` (the true optimum lies at or
        beyond the bound, and the supplied value is merely interior to
        it).  A flat-but-slightly-tilted profile, whose endpoint barely
        exceeds the optimum, is reported by the half-width instead and is
        *not* flagged here.
        """
        p = self.parameter
        f = p.normalize(self.x_opt)
        if f <= tol or f >= 1.0 - tol:
            return True
        if not np.any(np.isfinite(self.score)):
            return False
        best = int(np.nanargmax(self.score))
        return bool((best == 0 or best == len(self.score) - 1)
                    and self.score[best] > self.score_opt + delta)

    def __repr__(self):
        return (f"Profile({self.name!r}, n={len(self.x)}, "
                f"half_width={self.half_width():.3f}, "
                f"at_bound={self.at_bound()})")


def profile(objective, pset, name, n=21, span="bounds", width=0.5):
    """Sweep one parameter, holding the rest at the optimum.

    Parameters
    ----------
    objective : callable
        ``f(theta: dict) -> float`` mapping a ``{name: calibration_value}``
        dict to a score (higher better).  Non-finite returns are allowed
        (model failures) and are carried through as ``nan``.
    pset : ParameterSet
        The optimum; every parameter except ``name`` is held at its value.
    name : str
        Parameter to sweep.
    n : int, optional
        Number of sweep samples (default 21).
    span : {'bounds', 'local'}, optional
        ``'bounds'`` sweeps the full ``[lower, upper]`` (reveals
        bound-pegging; matches the f_route sweep).  ``'local'`` sweeps a
        symmetric window of half-width ``width`` × range about the
        optimum (sharper curvature resolution).
    width : float, optional
        Half-width of the ``'local'`` window as a fraction of range.

    Returns
    -------
    Profile
    """
    p = pset[name]
    if span == "bounds":
        x = np.linspace(p.lower, p.upper, n)
    elif span == "local":
        half = width * p.range
        x = np.linspace(max(p.lower, p.value - half),
                        min(p.upper, p.value + half), n)
    else:
        raise ValueError("span must be 'bounds' or 'local'")

    theta = pset.optimum()
    scores = np.empty(n, dtype=float)
    for i, xi in enumerate(x):
        theta[name] = xi
        try:
            s = objective(theta)
            scores[i] = s if np.isfinite(s) else np.nan
        except Exception:
            scores[i] = np.nan
    theta[name] = p.value  # restore

    score_opt = objective(pset.optimum())
    return Profile(p, x, scores, p.value, score_opt)


def profile_all(objective, pset, n=21, span="bounds", width=0.5, aic=None,
                verbose=False):
    """Profile every parameter in ``pset`` and assemble a report.

    See :func:`profile` for the sweep parameters.  ``aic`` (e.g.
    ``CalibResult.aic`` at the optimum) is carried into the report.
    """
    profiles = {}
    for p in pset:
        if verbose:
            print(f"  profiling {p.name} ...", flush=True)
        profiles[p.name] = profile(objective, pset, p.name, n=n,
                                   span=span, width=width)
    return IdentifiabilityReport(pset, profiles, aic=aic)


class IdentifiabilityReport:
    """Per-parameter identifiability profiles plus calibration metadata.

    The ``predicted`` attribute is the docking point for the *a-priori*
    structural prediction (diffusion-length for matrix-τ, drainage
    topology for ``f_route``).  It is ``None`` until those
    forward-measurement tools exist; when they do, the reconciliation of
    predicted-silent against empirically-flat parameters becomes the
    payoff — and tells you which field measurements are worth making.
    """

    def __init__(self, pset, profiles, aic=None, predicted=None, spectrum=None):
        self.pset      = pset
        self.profiles  = profiles
        self.aic       = aic
        self.predicted = predicted   # a-priori prediction seam (deferred)
        self.spectrum  = spectrum    # Spectrum of named combinations (optional)

    def summary(self, delta=0.01):
        """Return a human-readable identifiability table as a string."""
        lines = []
        lines.append("MNiShed identifiability report")
        lines.append("=" * 64)
        if self.aic is not None:
            lines.append(f"calibration AIC: {self.aic:.2f}")
        lines.append(
            f"Δ-flat half-width = fraction of bound range within "
            f"{delta:g} of peak score")
        lines.append("(large half-width or 'pegged' ⇒ poorly constrained by data)")
        lines.append("-" * 64)
        lines.append(f"{'parameter':<28}{'half-width':>12}{'curvature':>12}"
                     f"{'flag':>10}")
        lines.append("-" * 64)
        for name in self.pset.names:
            pr = self.profiles[name]
            hw = pr.half_width(delta)
            cv = pr.curvature()
            flag = "PEGGED" if pr.at_bound() else (
                "flat" if (np.isfinite(hw) and hw > 0.4) else "")
            hw_s = f"{hw:>12.3f}" if np.isfinite(hw) else f"{'nan':>12}"
            cv_s = f"{cv:>12.2g}" if np.isfinite(cv) else f"{'nan':>12}"
            lines.append(f"{name:<28}{hw_s}{cv_s}{flag:>10}")
        lines.append("-" * 64)
        if self.spectrum is not None:
            lines.append("")
            lines.append(self.spectrum.summary())
        if self.predicted is None:
            lines.append("")
            lines.append("a-priori prediction: (not available — awaits "
                         "Ksat / drainage-topology tools)")
        return "\n".join(lines)

    def __repr__(self):
        return (f"IdentifiabilityReport(n_params={len(self.pset)}, "
                f"aic={self.aic})")


# ======================================================================
# Curvature eigenspectrum — naming the degenerate *combinations*
# ======================================================================

class Spectrum:
    """Eigendecomposition of the objective's local curvature at the optimum.

    Built from the finite-difference Hessian of the score in *normalized*
    coordinates ``u_i = (x_i - x_opt_i) / range_i`` (so every parameter is
    measured as a fraction of its bound range and the eigenvectors are
    directly comparable).  At a maximum the score-Hessian is negative
    (semi-)definite, so this works with the **information matrix**
    ``M = -H``, whose eigenvalues are non-negative.

    A large eigenvalue is a *stiff* direction — a parameter combination
    the data pins sharply.  A small eigenvalue is a *sloppy* direction —
    a combination the data barely constrains; its eigenvector **names**
    the degenerate combination (the ``f_route`` ↔ land-τ ridge would
    appear as a sloppy eigenvector with weight on both, which no 1-D
    profile can reveal).

    Attributes
    ----------
    names : list of str
        Parameter order of the rows/columns.
    M : ndarray (n, n)
        The information matrix in normalized coordinates.
    eigvals : ndarray (n,)
        Eigenvalues (stiffness), sorted descending.
    eigvecs : ndarray (n, n)
        Columns are the corresponding eigenvectors, in ``names`` order.
    step : float
        Normalized finite-difference step used.
    """

    def __init__(self, names, M, step):
        full_names = list(names)
        self.M = np.asarray(M, dtype=float)     # full matrix, retained
        self.step = step
        # A parameter pegged at a bound has no room for a finite-difference
        # step, so its whole row/column is nan.  Drop those parameters
        # (they carry no local curvature information) and decompose the
        # remaining submatrix, rather than letting one pegged parameter —
        # e.g. a matrix-τ against its upper bound — void the entire
        # spectrum.  Excluded names are reported separately.
        diag = np.diag(self.M)
        keep = np.where(np.isfinite(diag))[0]
        self.names = [full_names[i] for i in keep]
        self.excluded = [full_names[i] for i in range(len(full_names))
                         if i not in keep]
        sub = self.M[np.ix_(keep, keep)] if keep.size else np.empty((0, 0))
        # residual nans (e.g. a model failure at one corner) → 0: treat that
        # off-diagonal coupling as unresolved rather than discarding the row.
        if sub.size and not np.all(np.isfinite(sub)):
            sub = np.where(np.isfinite(sub), sub, 0.0)
        sub = 0.5 * (sub + sub.T)               # symmetrise FD noise
        if sub.size:
            w, V = np.linalg.eigh(sub)
            order = np.argsort(w)[::-1]          # descending: stiff first
            self.eigvals = w[order]
            self.eigvecs = V[:, order]
        else:
            self.eigvals = np.array([])
            self.eigvecs = np.empty((0, 0))

    @property
    def condition_number(self):
        """Stiffest / sloppiest eigenvalue — overall identifiability span.

        Large condition numbers are the quantitative face of
        equifinality: a wide spread means a few combinations do all the
        constraining while others float.
        """
        w = self.eigvals
        w = w[np.isfinite(w)]
        w = w[w > 0]
        if w.size < 2:
            return np.nan
        return w.max() / w.min()

    def flat_directions(self, rel_tol=1e-3):
        """Sloppy eigen-directions: eigenvalue < ``rel_tol`` × the stiffest.

        Returns a list of ``(eigenvalue, {name: loading})`` for each
        direction below the threshold, loadings sorted by magnitude — the
        named degenerate combinations the data cannot separate.
        """
        w = self.eigvals
        finite = w[np.isfinite(w)]
        if finite.size == 0:
            return []
        cutoff = rel_tol * np.nanmax(w)
        out = []
        for k in range(len(w)):
            if np.isfinite(w[k]) and w[k] < cutoff:
                v = self.eigvecs[:, k]
                loading = {n: float(v[i]) for i, n in enumerate(self.names)}
                loading = dict(sorted(loading.items(),
                                      key=lambda kv: -abs(kv[1])))
                out.append((float(w[k]), loading))
        return out

    def describe_direction(self, k, top=3):
        """One-line label for eigen-direction ``k`` (its top loadings)."""
        v = self.eigvecs[:, k]
        idx = np.argsort(-np.abs(v))[:top]
        terms = [f"{v[i]:+.2f}·{self.names[i]}" for i in idx]
        return "  ".join(terms)

    def summary(self, rel_tol=1e-3):
        lines = ["curvature eigenspectrum (normalized coords)",
                 "-" * 64]
        cn = self.condition_number
        lines.append(f"condition number (stiff/sloppy): "
                     f"{cn:.3g}" if np.isfinite(cn) else
                     "condition number: nan")
        if self.excluded:
            lines.append("excluded (pegged at bound, no curvature): "
                         + ", ".join(self.excluded))
        lines.append(f"{'#':>2} {'stiffness':>12}  direction (top loadings)")
        for k in range(len(self.eigvals)):
            w = self.eigvals[k]
            tag = "  <-- SLOPPY" if (np.isfinite(w) and np.isfinite(cn)
                                     and w < rel_tol * np.nanmax(self.eigvals)) else ""
            ws = f"{w:>12.3g}" if np.isfinite(w) else f"{'nan':>12}"
            lines.append(f"{k:>2} {ws}  {self.describe_direction(k)}{tag}")
        return "\n".join(lines)

    def __repr__(self):
        return (f"Spectrum(n={len(self.names)}, "
                f"condition_number={self.condition_number:.3g})")


def _hessian_normalized(cobj, pset, step):
    """Finite-difference information matrix M = -H(score) in u-coords.

    Per-direction steps are shrunk to stay inside the bounds, so a
    near-bound parameter degrades gracefully (its row/column goes ``nan``)
    rather than stepping outside its prior.
    """
    names = pset.names
    n = len(names)
    opt = pset.optimum()
    ranges = {p.name: p.range for p in pset}

    # per-direction normalized step, clipped to 0.9 * room-to-nearest-bound
    h = np.empty(n)
    for i, p in enumerate(pset):
        f = p.normalize(p.value)             # optimum as fraction of range
        room = 0.9 * min(f, 1.0 - f)
        h[i] = step if room >= step else (room if room > 1e-6 else 0.0)

    def at(deltas):
        """Evaluate score at optimum + sum_i deltas[i]*range_i (u-shift)."""
        theta = dict(opt)
        for i, name in enumerate(names):
            theta[name] = opt[name] + deltas[i] * ranges[name]
        return cobj(theta)

    def shift(*axes):
        """Normalized-shift vector from (index, signed_step) pairs."""
        d = np.zeros(n)
        for idx, s in axes:
            d[idx] = s
        return d

    f0 = at(np.zeros(n))
    M = np.full((n, n), np.nan)

    # diagonal: second derivative of score, negated
    for i in range(n):
        if h[i] == 0.0:
            continue
        fp = at(shift((i, h[i])))
        fm = at(shift((i, -h[i])))
        M[i, i] = -(fp - 2.0 * f0 + fm) / (h[i] ** 2)

    # off-diagonal: mixed partial, negated
    for i in range(n):
        for j in range(i + 1, n):
            if h[i] == 0.0 or h[j] == 0.0:
                continue
            mij = -(at(shift((i, h[i]), (j, h[j])))
                    - at(shift((i, h[i]), (j, -h[j])))
                    - at(shift((i, -h[i]), (j, h[j])))
                    + at(shift((i, -h[i]), (j, -h[j])))) \
                / (4.0 * h[i] * h[j])
            M[i, j] = mij
            M[j, i] = mij
    return names, M


def eigenspectrum(objective, pset, step=0.05):
    """Curvature eigenspectrum of the objective at the optimum.

    Parameters
    ----------
    objective : callable
        ``f(theta: dict) -> float`` score (higher better); see
        :func:`profile`.
    pset : ParameterSet
        The optimum.
    step : float, optional
        Finite-difference step as a fraction of each bound range
        (default 0.05).  Because MNiShed's thresholds make the objective
        non-smooth, the spectrum is step-sensitive — compare two steps
        (e.g. ``0.05`` and ``0.025``) and trust the 2-D ridge grids for
        any direction the two disagree on.

    Returns
    -------
    Spectrum
    """
    cobj = _CachedObjective(objective, pset.names)
    names, M = _hessian_normalized(cobj, pset, step)
    return Spectrum(names, M, step)


# ======================================================================
# 2-D ridge grids — robust, non-smooth-safe view of a flagged combination
# ======================================================================

class Ridge2D:
    """A 2-D score grid over two parameters, the rest held at the optimum.

    Visualises a degenerate combination directly — without trusting a
    finite-difference Hessian — and estimates the ridge orientation as
    the principal axis of the near-optimal region (the set of grid points
    within ``delta`` of the best score).  That orientation, in normalized
    coordinates, can be read straight against the sloppy eigenvector from
    :class:`Spectrum` as a cross-check.
    """

    def __init__(self, pi, pj, xi, xj, Z):
        self.pi, self.pj = pi, pj
        self.name_i, self.name_j = pi.name, pj.name
        self.xi = np.asarray(xi, dtype=float)     # axis i values
        self.xj = np.asarray(xj, dtype=float)     # axis j values
        self.Z = np.asarray(Z, dtype=float)       # Z[a, b] at (xi[a], xj[b])

    def ridge_axis(self, delta=0.01):
        """Principal axis (in normalized i,j coords) of the near-optimal band.

        Returns ``(vi, vj)`` a unit vector, or ``None`` if the band is too
        small.  A band elongated along, say, ``(+0.7, -0.7)`` says the two
        parameters trade off one-for-one — the ridge.
        """
        zmax = np.nanmax(self.Z)
        A, B = np.meshgrid(self.pi.normalize(self.xi),
                           self.pj.normalize(self.xj), indexing="ij")
        mask = np.isfinite(self.Z) & (self.Z >= zmax - delta)
        pts = np.column_stack([A[mask], B[mask]])
        if len(pts) < 3:
            return None
        pts = pts - pts.mean(axis=0)
        _, _, Vt = np.linalg.svd(pts, full_matrices=False)
        v = Vt[0]
        return float(v[0]), float(v[1])

    def __repr__(self):
        ax = self.ridge_axis()
        ax_s = (f"({ax[0]:+.2f}, {ax[1]:+.2f})" if ax else "n/a")
        return (f"Ridge2D({self.name_i!r}, {self.name_j!r}, "
                f"grid={self.Z.shape}, ridge_axis={ax_s})")


def ridge(objective, pset, name_i, name_j, n=21, span="bounds", width=0.5):
    """Compute a 2-D score grid over two parameters at the optimum.

    See :func:`profile` for ``span``/``width``.  The off-axis parameters
    are held at their optimum values throughout.

    Returns
    -------
    Ridge2D
    """
    pi, pj = pset[name_i], pset[name_j]

    def axis(p):
        if span == "bounds":
            return np.linspace(p.lower, p.upper, n)
        elif span == "local":
            half = width * p.range
            return np.linspace(max(p.lower, p.value - half),
                               min(p.upper, p.value + half), n)
        raise ValueError("span must be 'bounds' or 'local'")

    xi, xj = axis(pi), axis(pj)
    theta = pset.optimum()
    Z = np.empty((n, n), dtype=float)
    for a, va in enumerate(xi):
        for b, vb in enumerate(xj):
            theta[name_i] = va
            theta[name_j] = vb
            try:
                s = objective(theta)
                Z[a, b] = s if np.isfinite(s) else np.nan
            except Exception:
                Z[a, b] = np.nan
    theta[name_i], theta[name_j] = pi.value, pj.value
    return Ridge2D(pi, pj, xi, xj, Z)
