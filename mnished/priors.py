"""
mnished.priors
~~~~~~~~~~~~~~~~~~
Data-driven prior estimation for MNiShed calibration.

:func:`suggest_priors` wraps :class:`~mnished.BrutsaertNieber` and
:class:`~mnished.HydrographSeparation` to produce a coherent set of
parameter starting points from an observed discharge record, before any
model run or calibration is attempted.

The returned :class:`Priors` object holds:

* Reservoir e-folding timescales (from spectral + recession decomposition).
* Power-law recession exponents (from the B–N recession cloud).
* Initial reservoir storage depths (from hydrograph separation).

It also exposes the underlying analysis objects for deeper inspection and
prints a human-readable summary.

Example
-------
>>> import pandas as pd
>>> from mnished import suggest_priors
>>> df = pd.read_csv('input.csv', parse_dates=['Date'])
>>> Q = df['Specific Discharge [mm/day]'].values
>>> P = df['Precipitation [mm/day]'].values
>>> pr = suggest_priors(Q, P=P, n_reservoirs=3)
>>> pr.summary()
"""

import warnings

import numpy as np

from .recession import BrutsaertNieber
from .hydrograph_separation import HydrographSeparation

# Theoretical Brutsaert & Nieber (1977) long-time baseflow recession exponent
# (b_HR in MNiShed notation, corresponding to B-N slope n ≈ 1.55).
# Appropriate as a fixed prior for the slow (karst/groundwater) reservoir.
_B_KARST_THEORETICAL = 2.203


class Priors:
    """
    Data-driven parameter priors for MNiShed.

    Produced by :func:`suggest_priors`; not intended to be instantiated
    directly.

    Attributes
    ----------
    t_recession : list of float
        Suggested e-folding residence times [days], fastest reservoir first.
        ``None`` entries indicate that the timescale could not be estimated
        from the data (fall back to calibration defaults).
    recession_exponents : list of float
        Suggested power-law recession exponents, fastest reservoir first.
        The fastest reservoir uses the catchment-integrated B–N estimate;
        the slow (karst) reservoir uses the theoretical value
        (Brutsaert & Nieber 1977, b ≈ 2.203); deeper reservoirs default
        to 1.0 (linear).
    initial_depths : list of float
        Estimated initial storage depths [mm], fastest reservoir first.
    log_t_recession_bounds : dict
        Calibration bounds in log10(days) for each ``log__t_recession_*``
        parameter, as returned by
        :meth:`~mnished.HydrographSeparation.get_parameter_priors`.
    bn : BrutsaertNieber
        Fitted Brutsaert & Nieber object; inspect ``bn.n_``, ``bn.a_``,
        or call ``bn.plot()`` for the recession cloud.
    hs : HydrographSeparation
        Fitted HydrographSeparation object; call ``hs.summary()`` for
        full spectral decomposition detail.
    n_reservoirs : int
        Number of reservoirs these priors are intended for.
    """

    def __init__(self, t_recession, recession_exponents, initial_depths,
                 log_t_recession_bounds, bn, hs, n_reservoirs):
        """
        Parameters
        ----------
        t_recession : list of float or None
            Recession timescales [days], fastest reservoir first. ``None``
            entries indicate the timescale could not be estimated.
        recession_exponents : list of float
            Power-law recession exponents, fastest reservoir first.
        initial_depths : list of float
            Initial storage depths [mm], fastest reservoir first.
        log_t_recession_bounds : dict
            Calibration bounds in log10(days) keyed by parameter name.
        bn : BrutsaertNieber
            Fitted recession analysis object.
        hs : HydrographSeparation
            Fitted hydrograph separation object.
        n_reservoirs : int
            Number of reservoirs these priors cover.
        """
        self.t_recession              = t_recession
        self.recession_exponents  = recession_exponents
        self.initial_depths       = initial_depths
        self.log_t_recession_bounds   = log_t_recession_bounds
        self.bn                   = bn
        self.hs                   = hs
        self.n_reservoirs         = n_reservoirs

    def summary(self):
        """
        Print a human-readable summary of all suggested priors.

        Includes timescales, recession exponents, initial depths, and
        calibration bounds, alongside guidance on which values to fix
        versus calibrate.
        """
        n = self.n_reservoirs
        labels = (['soil'] + ['karst'] + ['deep'] * (n - 2))[:n]

        print("=" * 60)
        print("MNiShed data-driven priors")
        print("=" * 60)

        print("\nE-folding residence times (fastest → slowest):")
        for label, tau in zip(labels, self.t_recession):
            if tau is None:
                print(f"  {label:8s}: could not be estimated — use calibration default")
            else:
                print(f"  {label:8s}: {tau:.1f} days")

        print("\nPower-law recession exponents (fastest → slowest):")
        b_note = {
            labels[0]:  "B–N data-driven estimate — calibrate",
            labels[-2] if n > 1 else labels[0]:
                        f"theoretical B–N 2.203 — consider fixing",
        }
        if n > 2:
            b_note[labels[-1]] = "linear (b=1) — deep reservoir default"
        for label, b in zip(labels, self.recession_exponents):
            note = b_note.get(label, "")
            print(f"  {label:8s}: {b:.3f}  ({note})")

        print("\nInitial storage depths [mm] (fastest → slowest):")
        for label, h0 in zip(labels, self.initial_depths):
            print(f"  {label:8s}: {h0:.1f} mm")

        print("\nCalibration bounds (log10 days) for params.yml:")
        if self.log_t_recession_bounds:
            for key, bounds in self.log_t_recession_bounds.items():
                if bounds is None:
                    print(f"  {key}: not estimated — keep params.yml defaults")
                else:
                    print(f"  {key}:  initial={bounds['initial']:.3f},"
                          f"  lower={bounds['lower']:.3f},"
                          f"  upper={bounds['upper']:.3f}")
        else:
            print("  (not available)")

        print("\nBrutsaert–Nieber recession cloud:")
        print(f"  slope n   = {self.bn.n_:.4f}")
        print(f"  coeff a   = {self.bn.a_:.4g}")
        print(f"  → b_HR    = {self.bn.to_reservoir_exponent():.3f}  "
              f"(used as b_{labels[0]} prior)")
        print("=" * 60)

    def to_yaml_snippet(self):
        """
        Return a YAML string snippet for the ``reservoirs:`` and
        ``initial_conditions:`` sections, populated with these priors.

        The snippet is a starting point only; review and adjust before
        running calibration.

        Returns
        -------
        str
        """
        n = self.n_reservoirs
        labels = (['soil'] + ['karst'] + ['deep'] * (n - 2))[:n]

        def _fmt(v):
            return f"{v:.1f}" if v is not None else "null  # could not be estimated"

        tau_lines  = "\n".join(f"        - {_fmt(t)}  # {l}" for t, l in
                                zip(self.t_recession, labels))
        exfilt_lines = "\n".join(
            f"        - {0.8 if i == 0 else (0.5 if i < n - 1 else 1.0)}"
            f"  # {l} — placeholder; calibrate"
            for i, l in enumerate(labels)
        )
        hmax_lines = "\n".join(f"        - .inf  # {l}" for l in labels)
        b_lines    = "\n".join(f"        - {b:.3f}  # {l}" for b, l in
                                zip(self.recession_exponents, labels))
        h0_lines   = "\n".join(f"        - {h:.1f}  # {l}"
                                for h, l in zip(self.initial_depths, labels))

        snippet = (
            f"reservoirs:\n"
            f"    recession_timescales__days:\n{tau_lines}\n"
            f"    exfiltration_fractions:\n{exfilt_lines}\n"
            f"    maximum_effective_depths__mm:\n{hmax_lines}\n"
            f"    recession_exponents:\n{b_lines}\n"
            f"\ninitial_conditions:\n"
            f"    water_reservoir_effective_depths__mm:\n{h0_lines}\n"
            f"    snowpack__mm_SWE: 0\n"
        )
        return snippet


def suggest_priors(Q, P=None, n_reservoirs=3, dt=1.0,
                   min_recession_days=3):
    """
    Estimate MNiShed parameter priors from an observed discharge record.

    Combines :class:`~mnished.BrutsaertNieber` recession analysis with
    :class:`~mnished.HydrographSeparation` to produce timescale
    estimates, recession exponents, and initial storage depths without
    running any model.

    Parameters
    ----------
    Q : array-like
        Observed specific discharge time series [mm/day]. Must be
        non-negative and at daily resolution.
    P : array-like, optional
        Observed precipitation time series [mm/day], same length as *Q*.
        Improves the HydrographSeparation spectral decomposition when
        provided. Default ``None``.
    n_reservoirs : int, optional
        Number of reservoirs in the intended model structure (2 or 3).
        Default ``3``.
    dt : float, optional
        Timestep [days]. Default ``1.0``.
    min_recession_days : int, optional
        Minimum recession length passed to
        :class:`~mnished.BrutsaertNieber`. Default ``3``.

    Returns
    -------
    Priors
        Object containing timescales, recession exponents, initial depths,
        calibration bounds, and the underlying analysis objects.

    Notes
    -----
    **Recession exponents** are assigned as follows:

    * *Fastest reservoir* (soil): the catchment-integrated B–N estimate
      ``b = 1 / (2 − n)``, where *n* is the slope of the log(−dQ/dt) vs
      log(Q) cloud.  Treat as a calibration starting point.
    * *Slow reservoir* (karst): the theoretical Brutsaert & Nieber (1977)
      long-time value, b ≈ 2.203.  Consider fixing rather than calibrating.
    * *Deep reservoir* (if present): b = 1.0 (linear).

    **Timescales** come from the spectral + recession decomposition in
    :class:`~mnished.HydrographSeparation`.  ``None`` entries in
    ``Priors.t_recession`` mean the timescale could not be resolved from the
    data; fall back to calibration defaults in that case.

    The B–N slope and the overall b estimate reflect the *catchment-
    integrated* recession, not individual reservoir responses.  They
    anchor the soil exponent but cannot distinguish the karst component
    directly — that is why the karst exponent defaults to the theoretical
    value.

    Examples
    --------
    >>> import pandas as pd
    >>> from mnished import suggest_priors
    >>> df = pd.read_csv('input.csv', parse_dates=['Date'])
    >>> Q  = df['Specific Discharge [mm/day]'].values
    >>> P  = df['Precipitation [mm/day]'].values
    >>> pr = suggest_priors(Q, P=P, n_reservoirs=3)
    >>> pr.summary()
    >>> print(pr.to_yaml_snippet())
    """
    Q = np.asarray(Q, dtype=float)
    if P is not None:
        P = np.asarray(P, dtype=float)

    # --- Brutsaert & Nieber recession analysis ---------------------------
    bn = BrutsaertNieber(Q, dt=dt, min_recession_days=min_recession_days)
    try:
        bn.fit()
        b_fast = bn.to_reservoir_exponent()
        if not np.isfinite(b_fast):
            warnings.warn(
                "B–N slope n ≥ 2; cannot convert to a finite recession "
                "exponent.  Defaulting b_soil prior to 2.0.",
                UserWarning, stacklevel=2,
            )
            b_fast = 2.0
    except ValueError as e:
        warnings.warn(
            f"BrutsaertNieber fit failed ({e}); defaulting b_soil prior to 2.0.",
            UserWarning, stacklevel=2,
        )
        b_fast = 2.0

    # Build recession_exponents list: soil / karst / deep(linear)
    if n_reservoirs == 1:
        recession_exponents = [b_fast]
    elif n_reservoirs == 2:
        recession_exponents = [b_fast, _B_KARST_THEORETICAL]
    else:
        recession_exponents = ([b_fast, _B_KARST_THEORETICAL]
                               + [1.0] * (n_reservoirs - 2))

    # --- Hydrograph separation ------------------------------------------
    hs = HydrographSeparation(Q, n_reservoirs=n_reservoirs, precip=P)
    try:
        hs.fit()
        ic     = hs.get_initial_conditions()
        bounds = hs.get_parameter_priors()
        h0_list = ic['H0']                        # fastest first
        # Convert log-scale bounds back to linear for t_recession display
        t_recession = []
        for i in range(n_reservoirs):
            key = list(bounds.keys())[i] if i < len(bounds) else None
            if key and bounds[key] is not None:
                t_recession.append(round(10 ** bounds[key]['initial'], 1))
            else:
                t_recession.append(None)
    except Exception as e:
        warnings.warn(
            f"HydrographSeparation fit failed ({e}); "
            f"timescales and initial depths will be None.",
            UserWarning, stacklevel=2,
        )
        h0_list = [None] * n_reservoirs
        t_recession = [None] * n_reservoirs
        bounds  = {}

    return Priors(
        t_recession             = t_recession,
        recession_exponents = recession_exponents,
        initial_depths      = [round(float(h), 1) if h is not None else None
                                for h in h0_list],
        log_t_recession_bounds  = bounds,
        bn                  = bn,
        hs                  = hs,
        n_reservoirs        = n_reservoirs,
    )
