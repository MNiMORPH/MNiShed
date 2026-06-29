"""
mnished.priors
~~~~~~~~~~~~~~~~~~
Data-driven prior estimation for MNiShed calibration.

:func:`suggest_priors` wraps :class:`~mnished.BrutsaertNieber` and
:class:`~mnished.HydrographSeparation` to produce a coherent set of
parameter starting points from an observed discharge record, before any
model run or calibration is attempted.

The returned :class:`Priors` object holds:

* Reservoir e-folding timescales (from spectral, time-domain, and recession decomposition).
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
import pandas as pd

from .hydrograph_separation import HydrographSeparation
from .recession import BrutsaertNieber

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
    recession_coeff : list of float
        Suggested recession coefficients [days], fastest reservoir first.
        These are e-folding timescales from hydrograph separation (b=1
        linear fit); for b>1 they are not residence times but still serve
        as initial values for calibration.  ``None`` entries indicate that
        the coefficient could not be estimated from the data (fall back to
        calibration defaults).
    recession_exponents : list of float
        Suggested power-law recession exponents, fastest reservoir first.
        The fastest reservoir uses the catchment-integrated B–N estimate;
        the slow (karst) reservoir uses the theoretical value
        (Brutsaert & Nieber 1977, b ≈ 2.203); deeper reservoirs default
        to 1.0 (linear).
    initial_depths : list of float
        Estimated initial storage depths [mm], fastest reservoir first.
    log_recession_coeff_bounds : dict
        Calibration bounds in log10(days) for each ``log__recession_coeff_*``
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

    def __init__(self, recession_coeff, recession_exponents, initial_depths,
                 log_recession_coeff_bounds, bn, hs, n_reservoirs):
        """
        Parameters
        ----------
        recession_coeff : list of float or None
            Recession coefficients [days], fastest reservoir first. ``None``
            entries indicate the coefficient could not be estimated.
        recession_exponents : list of float
            Power-law recession exponents, fastest reservoir first.
        initial_depths : list of float
            Initial storage depths [mm], fastest reservoir first.
        log_recession_coeff_bounds : dict
            Calibration bounds in log10(days) keyed by parameter name.
        bn : BrutsaertNieber
            Fitted recession analysis object.
        hs : HydrographSeparation
            Fitted hydrograph separation object.
        n_reservoirs : int
            Number of reservoirs these priors cover.
        """
        self.recession_coeff              = recession_coeff
        self.recession_exponents      = recession_exponents
        self.initial_depths           = initial_depths
        self.log_recession_coeff_bounds   = log_recession_coeff_bounds
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

        print("\nRecession coefficients (fastest → slowest):")
        for label, tau in zip(labels, self.recession_coeff):
            if tau is None:
                print(f"  {label:8s}: could not be estimated — use calibration default")
            else:
                print(f"  {label:8s}: {tau:.1f} days")

        print("\nPower-law recession exponents (fastest → slowest):")
        b_note = {
            labels[0]:  "B–N data-driven estimate — calibrate",
            labels[-2] if n > 1 else labels[0]:
                        "theoretical B–N 2.203 — consider fixing",
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
        if self.log_recession_coeff_bounds:
            for key, bounds in self.log_recession_coeff_bounds.items():
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
                                zip(self.recession_coeff, labels))
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
            f"    recession_coefficients:\n{tau_lines}\n"
            f"    exfiltration_fractions:\n{exfilt_lines}\n"
            f"    maximum_effective_depths__mm:\n{hmax_lines}\n"
            f"    recession_exponents:\n{b_lines}\n"
            f"\ninitial_conditions:\n"
            f"    water_reservoir_effective_depths__mm:\n{h0_lines}\n"
            f"    snowpack__mm_SWE: 0\n"
        )
        return snippet


def leafout_GDD_from_date(forcing, month, day, base_temperature__C=5.0,
                          tmax_col='Maximum Temperature [C]',
                          tmin_col='Minimum Temperature [C]',
                          date_col='Date', return_years=False):
    """
    Growing-degree-day leaf-out threshold from a regional green-up date.

    Converts a known canopy leaf-out *date* (``month``, ``day``) into the
    ``leafout_GDD`` phenology prior for a basin, by accumulating the basin's
    own daily growing-degree-days (GDD, base ``base_temperature__C``) from
    January 1 to that date each year and averaging across complete years.

    This grounds the green-up prior in the basin's actual thermal climate
    rather than a fabricated latitude-to-GDD curve. The leaf-out date is what
    spring-index phenology climatologies provide — e.g. the USA National
    Phenology Network Extended Spring Indices (the SI-x first-leaf index;
    Schwartz, Ault & Betancourt, 2013) — and it carries the latitude
    dependence implicitly: a later, more-northern green-up date integrates to a
    larger GDD against the same temperature forcing, so the resulting
    ``leafout_GDD`` is automatically larger at higher latitude. (Crow Wing,
    ~46–47°N, calibrates near 200 GDD against a ~late-May leaf-out, versus the
    generic ~100-GDD default tuned to a basin further south.)

    Parameters
    ----------
    forcing : pandas.DataFrame or str
        The basin temperature forcing, or a path to a CSV of it. Must contain
        a date column and daily maximum and minimum temperature columns.
    month, day : int
        The regional leaf-out date (e.g. ``5, 20`` for May 20). GDD is
        accumulated from January 1 to this date inclusive, each year.
    base_temperature__C : float, optional
        GDD base temperature [°C]. Default ``5.0`` — matches the phenology
        ``base_temperature__C`` default, so the returned threshold is on the
        same scale the model accumulates. Pass the value you use in the
        ``phenology:`` block if it differs.
    tmax_col, tmin_col, date_col : str, optional
        Column names in *forcing*. Defaults match the model's input columns.
    return_years : bool, optional
        If ``True``, also return the per-year GDD totals. Default ``False``.

    Returns
    -------
    float
        The mean accumulated GDD at the leaf-out date — use as the
        ``leafout_GDD`` prior in the ``phenology:`` config block.
    dict, optional
        ``{year: GDD}`` for each complete year, returned only when
        ``return_years`` is ``True``.

    Notes
    -----
    A year contributes only if its record reaches the leaf-out date; partial
    leading/trailing years are skipped (with a warning if any are). The early
    pre-spring days carry essentially no GDD (temperatures below the base), so
    a record that starts a little after January 1 changes the total
    negligibly. Use the *same* ``base_temperature__C`` here as in the
    ``phenology:`` block, or the threshold will not line up with the model's
    accumulation.

    Examples
    --------
    >>> import pandas as pd
    >>> from mnished import leafout_GDD_from_date
    >>> df = pd.read_csv('crow_wing_forcing.csv', parse_dates=['Date'])
    >>> # north-central MN canopy leaf-out ~ May 20 (USA-NPN SI-x)
    >>> leafout_GDD_from_date(df, 5, 20)            # doctest: +SKIP
    207.4
    """
    if isinstance(forcing, str):
        forcing = pd.read_csv(forcing, parse_dates=[date_col])
    if date_col not in forcing.columns:
        raise KeyError(
            f"leafout_GDD_from_date: forcing is missing column {date_col!r}.")
    dates = pd.DatetimeIndex(forcing[date_col])
    # Use the same daily mean the model would: the 'Mean Temperature [C]' column
    # when present (it may be a true integrated mean), else the min/max midpoint
    # the model synthesizes at load — so the prior lands on the model's GDD scale.
    if 'Mean Temperature [C]' in forcing.columns:
        tmean = np.asarray(forcing['Mean Temperature [C]'], dtype=float)
    else:
        for col in (tmax_col, tmin_col):
            if col not in forcing.columns:
                raise KeyError(
                    f"leafout_GDD_from_date: forcing needs 'Mean Temperature [C]' "
                    f"or both {tmax_col!r} and {tmin_col!r}; missing {col!r}.")
        tmean = 0.5 * (np.asarray(forcing[tmax_col], dtype=float)
                       + np.asarray(forcing[tmin_col], dtype=float))
    gdd_day = np.maximum(tmean - base_temperature__C, 0.0)
    years = dates.year.to_numpy()

    per_year = {}
    skipped = []
    for y in np.unique(years):
        try:
            target = pd.Timestamp(year=int(y), month=month, day=day)
        except ValueError as e:      # e.g. Feb 29 in a non-leap year
            raise ValueError(
                f"leafout_GDD_from_date: invalid leaf-out date "
                f"month={month}, day={day} for year {y} ({e}).")
        in_year = years == y
        # The record must reach the leaf-out date for the accumulation to be
        # complete; otherwise the partial sum under-counts that year.
        if dates[in_year].max() < target:
            skipped.append(int(y))
            continue
        mask = in_year & np.asarray(dates <= target)
        per_year[int(y)] = float(np.nansum(gdd_day[mask]))

    if not per_year:
        raise ValueError(
            "leafout_GDD_from_date: no year in the forcing reaches the "
            f"leaf-out date (month={month}, day={day}).")
    if skipped:
        warnings.warn(
            f"leafout_GDD_from_date: skipped {len(skipped)} year(s) whose "
            f"record does not reach the leaf-out date: {skipped}.",
            UserWarning, stacklevel=2)

    mean_gdd = float(np.mean(list(per_year.values())))
    if return_years:
        return mean_gdd, per_year
    return mean_gdd


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

    **Timescales** come from the spectral, time-domain, and recession decomposition in
    :class:`~mnished.HydrographSeparation`.  ``None`` entries in
    ``Priors.recession_coeff`` mean the coefficient could not be resolved from the
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
        # Convert log-scale bounds back to linear for recession_coeff display
        recession_coeff = []
        for i in range(n_reservoirs):
            key = list(bounds.keys())[i] if i < len(bounds) else None
            if key and bounds[key] is not None:
                recession_coeff.append(round(10 ** bounds[key]['initial'], 1))
            else:
                recession_coeff.append(None)
    except Exception as e:
        warnings.warn(
            f"HydrographSeparation fit failed ({e}); "
            f"timescales and initial depths will be None.",
            UserWarning, stacklevel=2,
        )
        h0_list = [None] * n_reservoirs
        recession_coeff = [None] * n_reservoirs
        bounds  = {}

    return Priors(
        recession_coeff             = recession_coeff,
        recession_exponents     = recession_exponents,
        initial_depths          = [round(float(h), 1) if h is not None else None
                                    for h in h0_list],
        log_recession_coeff_bounds  = bounds,
        bn                  = bn,
        hs                  = hs,
        n_reservoirs        = n_reservoirs,
    )
