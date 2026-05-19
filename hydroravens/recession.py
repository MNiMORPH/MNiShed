"""
hydroravens.recession
~~~~~~~~~~~~~~~~~~~~~
Brutsaert & Nieber (1977) recession analysis.

Fits the power-law relationship between recession rate and discharge::

    -dQ/dt = a * Q^n

on a log–log plot of all recession segments extracted from an observed
discharge time series.

The slope *n* characterises the nonlinearity of the storage–discharge
relationship.  In the HydroRaVENS power-law reservoir (Q ∝ H^b), the
two exponents are related by::

    n = (2b - 1) / b   →   b = 1 / (2 - n)

Special cases:

* n = 1, b = 1  — linear reservoir (exponential recession).
* n = 3/2, b = 2  — long-time Boussinesq groundwater solution
  (Brutsaert & Nieber 1977).
* n → 2, b → ∞  — extreme nonlinearity; storage empties abruptly.

The conversion is valid only for n < 2.

References
----------
Brutsaert, W. and Nieber, J. L. (1977). Regionalized drought flow
hydrographs from a mature glaciated plateau. Water Resources Research,
13(3), 637–643. https://doi.org/10.1029/WR013i003p00637

Kirchner, J. W. (2009). Catchments as simple dynamical systems:
Catchment characterization, rainfall-runoff modeling, and doing
hydrology backward. Water Resources Research, 45, W02429.
https://doi.org/10.1029/2008WR006912
"""

import warnings

import numpy as np


class BrutsaertNieber:
    """
    Brutsaert & Nieber (1977) recession analysis.

    Extracts all recession segments from an observed specific-discharge
    record, computes –dQ/dt vs Q for each consecutive pair within a
    segment, and fits a power law on the resulting log–log cloud::

        -dQ/dt = a * Q^n

    The fitted slope *n* can be converted to a HydroRaVENS
    ``recession_exponent`` *b* via :meth:`to_reservoir_exponent`.

    Parameters
    ----------
    Q : array-like
        Observed specific discharge time series [mm/day, or any consistent
        units]. Must be non-negative; zero values are skipped.
    dt : float, optional
        Timestep length [days]. Default ``1.0``.
    min_recession_days : int, optional
        Minimum number of consecutive declining days required for a
        segment to be included in the fit. Default ``3``.

    Attributes
    ----------
    n_ : float
        Fitted B–N slope (exponent in –dQ/dt = a·Q^n). Set by :meth:`fit`.
    a_ : float
        Fitted coefficient. Set by :meth:`fit`.
    r2_ : float
        Coefficient of determination of the power-law (linear log–log) fit.
        Set by :meth:`fit`.  Values below ~0.4 indicate a poor fit; a
        warning is issued automatically.
    r2_quadratic_ : float
        R² of a degree-2 polynomial fit in log–log space. Set by
        :meth:`fit`.  When substantially larger than :attr:`r2_`, the
        recession cloud is systematically curved, suggesting the power-law
        storage–discharge assumption does not hold; a warning is issued
        automatically.

    Examples
    --------
    >>> import numpy as np
    >>> from hydroravens import BrutsaertNieber
    >>> Q = np.array([10.0, 8.0, 6.5, 5.2, 4.0, 3.1, 2.4, 1.8, 1.3])
    >>> bn = BrutsaertNieber(Q).fit()
    >>> print(f"n = {bn.n_:.3f},  b_HR = {bn.to_reservoir_exponent():.3f}")
    """

    def __init__(self, Q, dt=1.0, min_recession_days=3):
        self.Q = np.asarray(Q, dtype=float)
        self.dt = dt
        self.min_recession_days = int(min_recession_days)
        self.n_ = None
        self.a_ = None
        self.r2_ = None
        self.r2_quadratic_ = None
        self._Q_mid = None
        self._neg_dQdt = None

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _extract_pairs(self):
        """
        Return (Q_mid, –dQ/dt) arrays for all qualifying recession steps.

        A recession segment is a maximal run of consecutive time steps
        where Q[t+1] < Q[t] and both values are positive.  Segments
        shorter than ``min_recession_days`` are discarded.
        """
        Q = self.Q
        n = len(Q)
        Q_mid_acc = []
        neg_dQ_acc = []

        i = 0
        while i < n - 1:
            if Q[i] > 0 and Q[i + 1] < Q[i] and Q[i + 1] > 0:
                # Walk to the end of this recession segment.
                j = i + 1
                while j < n - 1 and Q[j + 1] < Q[j] and Q[j + 1] > 0:
                    j += 1
                seg_len = j - i  # number of declining steps

                if seg_len >= self.min_recession_days - 1:
                    for k in range(i, j):
                        q_mid = 0.5 * (Q[k] + Q[k + 1])
                        neg_dq = (Q[k] - Q[k + 1]) / self.dt
                        if q_mid > 0 and neg_dq > 0:
                            Q_mid_acc.append(q_mid)
                            neg_dQ_acc.append(neg_dq)
                i = j  # skip to end of segment
            else:
                i += 1

        return np.array(Q_mid_acc), np.array(neg_dQ_acc)

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def fit(self):
        """
        Extract recession pairs, fit the log–log power law, and assess fit
        quality.

        Performs ordinary least-squares regression of log(–dQ/dt) on
        log(Q).  The resulting slope is stored in :attr:`n_` and the
        coefficient in :attr:`a_`.

        Two goodness-of-fit diagnostics are computed and stored:

        * :attr:`r2_` — R² of the linear (power-law) log–log fit.
        * :attr:`r2_quadratic_` — R² of a degree-2 polynomial log–log fit.

        A ``UserWarning`` is issued if the recession cloud shows systematic
        curvature (quadratic fit substantially better than linear, ΔR² > 0.05),
        or if the power-law fit is poor overall (R² < 0.4).  Curvature
        indicates that the power-law storage–discharge assumption may not
        hold; following Kirchner (2009), a non-parametric characterisation of
        the catchment sensitivity function should be considered before
        committing to a reservoir functional form.

        Returns
        -------
        self
            Allows method chaining (e.g. ``BrutsaertNieber(Q).fit()``).

        Raises
        ------
        ValueError
            If fewer than three recession pairs are found (cannot fit).
        """
        Q_mid, neg_dQdt = self._extract_pairs()
        if len(Q_mid) < 3:
            raise ValueError(
                f"Only {len(Q_mid)} recession pair(s) found "
                f"(min_recession_days={self.min_recession_days}); "
                f"need at least 3 to fit."
            )
        self._Q_mid = Q_mid
        self._neg_dQdt = neg_dQdt

        log_Q = np.log(Q_mid)
        log_R = np.log(neg_dQdt)

        # Linear (power-law) fit
        slope, intercept = np.polyfit(log_Q, log_R, 1)
        self.n_ = float(slope)
        self.a_ = float(np.exp(intercept))

        # R² of the power-law fit
        log_R_pred = slope * log_Q + intercept
        ss_res = np.sum((log_R - log_R_pred) ** 2)
        ss_tot = np.sum((log_R - np.mean(log_R)) ** 2)
        self.r2_ = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else 1.0

        # Quadratic fit in log–log space (curvature test)
        coeffs2 = np.polyfit(log_Q, log_R, 2)
        log_R_pred2 = np.polyval(coeffs2, log_Q)
        ss_res2 = np.sum((log_R - log_R_pred2) ** 2)
        self.r2_quadratic_ = float(1.0 - ss_res2 / ss_tot) if ss_tot > 0 else 1.0
        delta_r2 = self.r2_quadratic_ - self.r2_

        # Goodness-of-fit warnings
        if delta_r2 > 0.05:
            warnings.warn(
                f"B–N recession cloud shows systematic curvature in log–log "
                f"space (power-law R² = {self.r2_:.2f}, quadratic R² = "
                f"{self.r2_quadratic_:.2f}, ΔR² = {delta_r2:.2f}). "
                f"The power-law storage–discharge assumption may not hold for "
                f"this catchment. Following Kirchner (2009), consider "
                f"characterising the catchment sensitivity function "
                f"non-parametrically before choosing a reservoir functional "
                f"form. Inspect the recession plot (bn.plot()) for guidance.",
                UserWarning,
                stacklevel=2,
            )
        elif self.r2_ < 0.4:
            warnings.warn(
                f"B–N recession cloud is poorly described by a power law "
                f"(R² = {self.r2_:.2f}). "
                f"Inspect the recession plot (bn.plot()) to assess whether a "
                f"different functional form may be more appropriate for this "
                f"catchment.",
                UserWarning,
                stacklevel=2,
            )

        return self

    def to_reservoir_exponent(self):
        """
        Convert the fitted B–N slope *n* to a HydroRaVENS
        ``recession_exponent`` *b*.

        The relationship derives from substituting Q ∝ H^b and
        dH/dt = –Q into the B–N form –dQ/dt = a·Q^n::

            n = (2b - 1) / b   →   b = 1 / (2 - n)

        This conversion is valid only for n < 2.  For n ≥ 2, the
        equivalent power-law reservoir exponent is infinite (or
        ill-defined), and a ``UserWarning`` is issued.

        Returns
        -------
        float
            HydroRaVENS ``recession_exponent`` *b*.  Returns
            ``numpy.inf`` when n ≥ 2.

        Raises
        ------
        RuntimeError
            If :meth:`fit` has not been called.
        """
        if self.n_ is None:
            raise RuntimeError("Call fit() before to_reservoir_exponent().")
        if self.n_ >= 2.0:
            warnings.warn(
                f"Fitted B–N slope n = {self.n_:.3f} ≥ 2; the conversion "
                f"b = 1/(2 − n) is undefined.  Consider fitting only the "
                f"lower envelope of the recession cloud (long-time behaviour).",
                UserWarning,
                stacklevel=2,
            )
            return np.inf
        return 1.0 / (2.0 - self.n_)

    def summary(self):
        """
        Print a short summary of the fitted parameters.

        Raises
        ------
        RuntimeError
            If :meth:`fit` has not been called.
        """
        if self.n_ is None:
            raise RuntimeError("Call fit() before summary().")
        b_hr = self.to_reservoir_exponent()
        b_str = f"{b_hr:.3f}" if np.isfinite(b_hr) else "∞"
        delta_r2 = self.r2_quadratic_ - self.r2_
        curv_note = (f"  ** curvature detected (ΔR² = {delta_r2:.2f})"
                     if delta_r2 > 0.05 else "")
        print(
            f"Brutsaert–Nieber recession analysis\n"
            f"  Recession pairs used : {len(self._Q_mid)}\n"
            f"  Fitted slope  n      : {self.n_:.4f}\n"
            f"  Fitted coeff  a      : {self.a_:.4g}\n"
            f"  HydroRaVENS   b      : {b_str}\n"
            f"  Power-law R²         : {self.r2_:.3f}\n"
            f"  Quadratic R²         : {self.r2_quadratic_:.3f}"
            + (f"\n{curv_note}" if curv_note else "") +
            f"\n  Reference (long-time Boussinesq): n = 1.5, b = 2.0"
        )

    def plot(self, ax=None, show_fit=True, label=None, **scatter_kwargs):
        """
        Plot –dQ/dt vs Q on a log–log scale.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw on.  Created if not supplied.
        show_fit : bool, optional
            Overlay the fitted power-law line. Default ``True``.
        label : str, optional
            Label for the scatter points (appears in legend when
            ``show_fit=True``).
        **scatter_kwargs
            Passed to ``ax.scatter``.

        Returns
        -------
        matplotlib.axes.Axes

        Raises
        ------
        RuntimeError
            If :meth:`fit` has not been called.
        """
        if self._Q_mid is None:
            raise RuntimeError("Call fit() before plot().")

        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots()

        scatter_kwargs.setdefault('s', 4)
        scatter_kwargs.setdefault('alpha', 0.35)
        scatter_kwargs.setdefault('color', 'steelblue')
        ax.scatter(self._Q_mid, self._neg_dQdt, label=label, **scatter_kwargs)

        if show_fit and self.n_ is not None:
            Q_line = np.logspace(
                np.log10(self._Q_mid.min()),
                np.log10(self._Q_mid.max()),
                200,
            )
            b_hr = self.to_reservoir_exponent()
            b_str = f"{b_hr:.2f}" if np.isfinite(b_hr) else "∞"
            ax.plot(
                Q_line,
                self.a_ * Q_line ** self.n_,
                color='firebrick',
                linewidth=1.5,
                label=f"fit: n = {self.n_:.3f},  b = {b_str}",
            )

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Q  [mm day⁻¹]')
        ax.set_ylabel('−dQ/dt  [mm day⁻²]')
        ax.set_title('Brutsaert & Nieber (1977) recession analysis')
        ax.legend()
        return ax
