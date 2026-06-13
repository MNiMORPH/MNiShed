"""
hydroravens.calibration
~~~~~~~~~~~~~~~~~~~~~~~
Run hydroRaVENS with a given parameter set and return a CalibResult
named tuple containing the goodness-of-fit score, AIC, baseflow index,
flow duration curve, end-of-run reservoir states, and the Buckets object.

Intended use: call run_and_score() from a Dakota driver or any other
optimizer. Set spin_up_cycles: 0 in the YAML config; run_and_score()
manages spin-up itself using the calibrated parameters. If the YAML
requests spin-up cycles, initialize() will run them before the parameter
overrides are applied, so those cycles use uncalibrated values.

Supported metrics
-----------------
'NSE'    Nash-Sutcliffe Efficiency.  Biased toward high flows because its
         denominator is the variance of observed discharge (dominated by
         peaks).  Use as a baseline or when peaks are the primary concern.

'KGE'    Kling-Gupta Efficiency.  Decomposes fit into correlation (r),
         variability ratio (alpha = std_mod/std_obs), and bias ratio
         (beta = mean_mod/mean_obs), weighting all three equally.
         Better balanced across the full flow range than NSE.

'logKGE' KGE applied to log-transformed flows.  Shifts sensitivity toward
         low flows and base flow; useful when base flow matters as much as
         peaks.  A small epsilon (1 % of mean observed flow) is added
         before taking logs to avoid log(0).

'KGE_logKGE'
         Equal-weight average of KGE and logKGE: 0.5*KGE + 0.5*logKGE.
         Balances sensitivity to peaks (KGE) and low flows (logKGE).
         Recommended when neither regime should dominate calibration.
         Following Yilmaz et al. (2008, WRR), who show that no single
         metric captures both high- and low-flow behaviour; this composite
         is the single-objective analogue of their multi-segment FDC
         approach.

AIC
---
AIC = N * ln(SS_res_log / N) + 2k, where SS_res_log is the sum of squared
residuals on log-transformed flows and k is the number of free parameters
passed to run_and_score().  Log-transforming flows makes the Gaussian
residual assumption more defensible for discharge data.  AIC is intended
for comparing models with different numbers of reservoirs; lower is better.

Baseflow index (BFI)
--------------------
Computed with the Eckhardt (2005) recursive digital filter applied to both
observed and modelled specific discharge within the scoring window.
BFI = baseflow volume / total flow volume.  alpha=0.98 and bfi_max=0.80
are standard values for perennial daily streamflow.

Flow duration curve (FDC)
--------------------------
Discharge at 99 evenly-spaced exceedance probabilities (0.5–99.5 %).
Stored as pd.Series indexed by exceedance probability (%) in both
CalibResult.fdc_obs and CalibResult.fdc_mod.

Chaining decades
----------------
run_and_score() returns final_states, a dict of reservoir water depths and
snowpack SWE at the end of the scored run.  Pass this as initial_states to
the next decade's run_and_score() call with spin_up_cycles=0 so that water
storage is physically continuous across decade boundaries.
"""

from collections import namedtuple

import math
import warnings

import numpy as np
import pandas as pd

from .hydroravens import Buckets, Reservoir


# ---------------------------------------------------------------------------
# Named tuple for return value
# ---------------------------------------------------------------------------

CalibResult = namedtuple('CalibResult', [
    'score',        # float: KGE / NSE / logKGE (higher is better)
    'aic',          # float: AIC on log-transformed flows (lower is better)
    'bfi_obs',      # float: observed baseflow index
    'bfi_mod',      # float: modelled baseflow index
    'kge_logfdc',   # float: KGE on log-transformed FDC ordinates
    'fdc_obs',      # pd.Series: observed flow at exceedance probabilities
    'fdc_mod',      # pd.Series: modelled flow at exceedance probabilities
    'final_states', # dict: {'reservoirs': [...], 'snowpack': float, 'fgi': float}
    'buckets',      # Buckets object after the final run
])
CalibResult.__doc__ = """
Named tuple returned by :func:`run_and_score`.

Attributes
----------
score : float
    Goodness-of-fit score (higher is better). Metric is one of KGE,
    NSE, or logKGE as requested; see :func:`run_and_score`.
aic : float
    Akaike Information Criterion computed on log-transformed flows
    (lower is better). Useful for comparing models with different
    numbers of free parameters.
bfi_obs : float
    Baseflow index of the observed discharge, computed with the
    Eckhardt (2005) recursive digital filter.
bfi_mod : float
    Baseflow index of the modelled discharge.
fdc_obs : pd.Series
    Observed flow duration curve: discharge at 99 evenly-spaced
    exceedance probabilities (0.5–99.5 %), indexed by exceedance %.
fdc_mod : pd.Series
    Modelled flow duration curve (same format as fdc_obs).
final_states : dict
    End-of-run reservoir states suitable for passing as
    ``initial_states`` to the next call::

        {'reservoirs': [H_shallow, H_deep, ...],  # [mm]
         'snowpack':    H_snow_SWE,               # [mm SWE]
         'fgi':         frozen_ground_index}       # [°C·day]

buckets : Buckets
    The :class:`~hydroravens.Buckets` object after the final run,
    including the full ``hydrodata`` DataFrame with modelled discharge.

All scalar fields are ``np.nan`` if the scoring window contains no
valid overlapping data.
"""


# ---------------------------------------------------------------------------
# Metric helpers – operate on plain numpy arrays
# ---------------------------------------------------------------------------

def _nse(m, o):
    return float(1.0 - np.sum((m - o) ** 2) / np.sum((o - o.mean()) ** 2))


def _kge(m, o):
    r     = np.corrcoef(m, o)[0, 1]
    alpha = m.std() / o.std()
    beta  = m.mean() / o.mean()
    return float(1.0 - np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2))


def _log_kge(m, o):
    eps = 0.01 * o.mean()
    return _kge(np.log(m + eps), np.log(o + eps))


def _kge_logkge(m, o):
    return 0.5 * _kge(m, o) + 0.5 * _log_kge(m, o)


def _kge_logfdc(m, o):
    """KGE between log-transformed FDC ordinates (rank-sorted flows)."""
    eps  = 0.01 * o.mean()
    m_fdc = np.sort(m)[::-1]
    o_fdc = np.sort(o)[::-1]
    return _kge(np.log(m_fdc + eps), np.log(o_fdc + eps))


def _kge_logkge_logfdc(m, o):
    return (  _kge(m, o)
            + _log_kge(m, o)
            + _kge_logfdc(m, o)) / 3.0


def _kge_logkge_logfdc_bfi(m, o):
    bfi_score = 1.0 - abs(_eckhardt_bfi(m) / _eckhardt_bfi(o) - 1.0)
    return (  _kge(m, o)
            + _log_kge(m, o)
            + _kge_logfdc(m, o)
            + bfi_score) / 4.0


def _logkge_logfdc_bfi(m, o):
    bfi_score = 1.0 - abs(_eckhardt_bfi(m) / _eckhardt_bfi(o) - 1.0)
    return (  _log_kge(m, o)
            + _kge_logfdc(m, o)
            + bfi_score) / 3.0


def _aic(m, o, k):
    eps        = 0.01 * o.mean()
    ss_res_log = np.sum((np.log(m + eps) - np.log(o + eps)) ** 2)
    n          = len(o)
    return float(n * np.log(ss_res_log / n) + 2 * k)


def _eckhardt_bfi(q, alpha=0.98, bfi_max=0.80):
    """
    Eckhardt (2005) recursive digital filter for baseflow separation.
    Returns baseflow index = baseflow volume / total flow volume.
    alpha   : recession constant (~0.98 for daily perennial streams)
    bfi_max : maximum BFI (0.80 perennial; 0.50 ephemeral)
    """
    b     = np.empty_like(q, dtype=float)
    b[0]  = q[0] * bfi_max
    denom = 1.0 - alpha * bfi_max
    for t in range(1, len(q)):
        b[t] = ((1.0 - bfi_max) * alpha * b[t - 1]
                + (1.0 - alpha) * bfi_max * q[t]) / denom
        if b[t] > q[t]:
            b[t] = q[t]
    return float(b.sum() / q.sum())


_FDC_PROBS = np.arange(0.5, 100.0, 1.0)   # exceedance probabilities (%)


def _fdc(q):
    """pd.Series of discharge at standard exceedance probabilities."""
    flows = np.percentile(q, 100.0 - _FDC_PROBS)
    return pd.Series(flows, index=_FDC_PROBS, name='Specific discharge [mm/day]')


def _nash_cascade(q, N, K, dt=1.0):
    """
    Route a runoff time series through a Nash cascade of N identical
    linear reservoirs, each with storage time constant K [days].

    The cascade impulse response (instantaneous unit hydrograph) is a
    two-parameter gamma distribution:

        h(t) = t^(N-1) * exp(-t/K) / (K^N * Gamma(N))

    with mean travel time N*K [days] and variance N*K^2 [days^2].  For
    N = 1 the response reduces to a single exponential.

    Each reservoir is updated with the exact analytical solution for
    piecewise-constant inflow over a timestep dt:

        S_i(t+dt) = S_i(t) * exp(-dt/K)  +  K * Q_{i-1}(t) * (1 - exp(-dt/K))
        Q_i(t+dt) = S_i(t+dt) / K

    This form is unconditionally stable for any K > 0 and dt > 0, unlike
    the forward-Euler discretisation which requires dt/K < 2.

    Parameters
    ----------
    q : array-like, shape (T,)
        Runoff input time series [mm/day].
    N : int
        Number of linear reservoirs in series (shape parameter).
        N = 2 is typical for medium-sized catchments (area ~ 10^3 km^2).
    K : float
        Storage time constant of each reservoir [days] (scale parameter).
        Mean travel time through the cascade is N * K.
    dt : float, optional
        Timestep [days].  Default 1.0.

    Returns
    -------
    np.ndarray, shape (T,)
        Routed discharge time series [mm/day].

    References
    ----------
    Nash, J. E. (1957). The form of the instantaneous unit hydrograph.
        IAHS Publ. 45, 114–121.
        (Introduced the N-reservoir cascade and its gamma IUH.)

    Dooge, J. C. I. (1959). A general theory of the unit hydrograph.
        J. Geophys. Res., 64(2), 241–256.
        https://doi.org/10.1029/JZ064i002p00241

    Rodriguez-Iturbe, I. and Valdés, J. B. (1979). The geomorphologic
        structure of hydrologic response. Water Resour. Res., 15(6),
        1409–1420.
        https://doi.org/10.1029/WR015i006p01409
        (Shows that the gamma IUH arises naturally from Horton scaling
        laws, justifying its use without explicit flow-path geometry.)
    """
    q     = np.asarray(q, dtype=float)
    N     = int(round(N))
    alpha = np.exp(-dt / K)          # exact decay factor over one timestep
    beta  = K * (1.0 - alpha)        # input gain:  K*(1 - exp(-dt/K))

    # NaN in q (missing-forcing days) must not corrupt routing state.
    # Treat missing land-surface inflow as zero (channel drains normally);
    # restore NaN in the output so those days are excluded from scoring.
    # Caveat: substituting 0 drains the routing reservoir during long NaN
    # gaps; if the gap abuts the scoring window and K is large, early scored
    # discharge will be biased low.
    nan_mask = np.isnan(q)
    q_safe   = np.where(nan_mask, 0.0, q)

    S   = np.zeros(N)                # initial storage in each reservoir [mm]
    out = np.empty_like(q)

    for t in range(len(q_safe)):
        inflow = q_safe[t]
        for i in range(N):
            S[i]   = alpha * S[i] + beta * inflow
            inflow = S[i] / K        # outflow from reservoir i → inflow to i+1
        out[t] = inflow              # outflow from the final reservoir

    out[nan_mask] = np.nan
    return out


_METRICS = {'NSE': _nse, 'KGE': _kge, 'logKGE': _log_kge,
            'KGE_logKGE': _kge_logkge,
            'KGE_logKGE_logFDC': _kge_logkge_logfdc,
            'KGE_logKGE_logFDC_BFI': _kge_logkge_logfdc_bfi,
            'logKGE_logFDC_BFI': _logkge_logfdc_bfi}


def _steady_state_depths(reservoirs, mean_q):
    """
    Analytical steady-state water depth for each reservoir in a cascade.

    For a linear reservoir with constant mean inflow Q̄_in and dt = 1 day:

        H_eq = Q̄_in / (exp(1/τ) − 1)

    derived from the exact update H_{t+1} = (H_t + Q̄_in) · exp(−1/τ) and
    the steady-state condition H_{t+1} = H_t.

    Mean recharge to the top reservoir equals the long-run mean streamflow
    (mass conservation at steady state).  Recharge to each deeper reservoir
    is the fraction (1 − f_i) of the drainage from the reservoir above.

    Parameters
    ----------
    reservoirs : list of Reservoir
        Ordered shallowest to deepest; each has .t_efold, .f_to_discharge,
        and .Hmax attributes.
    mean_q : float
        Long-run mean specific discharge [mm/day], used as the steady-state
        mean recharge to the top reservoir.

    Returns
    -------
    list of float
        Steady-state water depth [mm] for each reservoir, capped at Hmax.

    Notes
    -----
    Using these depths as initial conditions serves two purposes: (1) it
    gives physically correct starting storage without relying on spin-up to
    fill reservoirs whose timescale exceeds the record length, and (2) it
    allows spin-up cycles to converge faster because short-timescale
    reservoirs start near equilibrium too -- only transient variability
    (seasonal cycles, wet/dry year sequences) needs to be resolved by
    spin-up rather than the slow drift from an arbitrary starting depth.
    """
    depths = []
    q_in = float(mean_q)
    for res in reservoirs:
        H_eq = q_in / (np.exp(1.0 / res.t_efold) - 1.0)
        depths.append(min(H_eq, res.Hmax))
        # Tile drainage reduces the recharge reaching the next reservoir.
        q_in *= (1.0 - res.f_to_discharge) * (1.0 - res.f_tile)
    return depths


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_and_score(cfg, t_efold=None, f_to_discharge=None, Hmax=None,
                  pdm_H0=None, f_tile=None, tau_tile=None,
                  melt_factor=None, fdd_threshold=None, snow_insulation_k=None,
                  et_scale=None, et_alpha=None,
                  wp_soil=None, wp_soil_sigma=None,
                  recession_exponents=None, recession_exponents_calibrated=0,
                  direct_runoff_fraction=None, baseflow_Q=None,
                  modules=None,
                  initial_states=None,
                  start=None, end=None, spin_up_cycles=None,
                  metric='KGE', routing_N=2, routing_K=None,
                  enforce_water_balance='water-year'):
    """
    Run hydroRaVENS and return a CalibResult named tuple.

    Parameters
    ----------
    cfg : str
        Path to a YAML configuration file. Should have spin_up_cycles: 0
        so that spin-up is performed here with the calibrated parameters.
    t_efold : list of float, optional
        E-folding residence times [days], one per reservoir (shallowest
        first). Overrides the values in cfg.
    f_to_discharge : list of float, optional
        Exfiltration fractions to stream, one per reservoir except the
        deepest (which is always 1.0). Overrides the values in cfg.
    Hmax : list of float, optional
        Maximum effective water depths [mm], one per reservoir. Overrides
        the values in cfg.
    melt_factor : float, optional
        Degree-day snowmelt factor [mm SWE per degC per day]. Overrides
        the value in cfg.
    direct_runoff_fraction : float or None, optional
        Fraction (0–1) of positive daily recharge that bypasses the
        reservoir cascade and exits directly as runoff.  Conceptually
        inspired by Hortonian (infiltration-excess) overland flow, but
        at a daily timestep rainfall intensity is unavailable, so the
        fraction is not a rigorous physical representation -- except in
        extreme events where intense rainfall dominates the daily total.
        In practice it acts as a calibrated fast-bypass fraction.
        None (default) leaves the value from the YAML config (itself
        defaulting to 0; off by default).
    fdd_threshold : float or None, optional
        Frozen ground index threshold [°C·day].  The frozen ground index
        accumulates freezing degree-days and decays during warming
        (Molnau & Bissell 1983, https://westernsnowconference.org/bibliography/1983Molnau.pdf).  When the index
        exceeds fdd_threshold, infiltration from the top reservoir to
        deeper layers is set to zero for that timestep (all drainage
        becomes direct runoff).  None (default) disables the effect.
    snow_insulation_k : float or None, optional
        Snow insulation decay constant [mm⁻¹ SWE].  Scales the effective
        air temperature reaching the soil as T_eff = T · exp(-k · SWE),
        reducing FGI accumulation under a deep snowpack.  Applied to
        both freezing and thawing temperature forcing; excess degree-days
        from meltwater (excess_dd) are not scaled because meltwater
        delivers heat directly to the soil surface.  None (default)
        leaves the value from cfg (itself defaulting to 0.0, i.e. no
        insulation).  Literature values: LISFLOOD uses ~0.057 mm⁻¹;
        GSSHA default is 0.5 (units ambiguous in original source).
    initial_states : dict, optional
        Starting reservoir water depths, snowpack SWE, and frozen ground
        index, as returned by a previous call's CalibResult.final_states.
        Format::

            {'reservoirs': [H_shallow, H_deep, ...],
             'snowpack':    H_snow_SWE,
             'fgi':         frozen_ground_index}

        When provided, these override the H0 values from cfg.  Use with
        spin_up_cycles=0 when chaining consecutive decades so that water
        storage and frozen-ground state are physically continuous.
        When None (default), reservoirs are initialised to their
        analytical steady-state depths before spin-up (see
        :func:`_steady_state_depths`).
    start : str or datetime-like, optional
        Start of the scoring window (inclusive). Score, AIC, BFI, and FDC
        are all computed within this window. Spin-up still uses the full
        record.
    end : str or datetime-like, optional
        End of the scoring window (inclusive). Same as start.
    spin_up_cycles : int or None, optional
        Number of times to loop the full record before the scored run.
        ``None`` (default) auto-computes as
        ``ceil(tau_max / record_length)``, where ``tau_max`` is the
        longest reservoir e-folding time after parameter overrides and
        ``record_length`` is the number of days in the input record.
        Because initial conditions are set to analytical steady-state
        depths, one e-folding time is sufficient to resolve the seasonal
        and inter-annual climate memory.  Set to 0 when providing
        initial_states for chained decade runs.
    metric : {'KGE', 'NSE', 'logKGE', 'KGE_logKGE', 'KGE_logKGE_logFDC', 'KGE_logKGE_logFDC_BFI'}, optional
        Goodness-of-fit metric.  Default is 'KGE'.
        ``'KGE_logKGE'`` returns 0.5*KGE + 0.5*logKGE, balancing
        peak and low-flow performance (Yilmaz et al. 2008).
        ``'KGE_logKGE_logFDC_BFI'`` adds a BFI bias-ratio score
        (1 - |BFI_mod/BFI_obs - 1|) as a fourth equal-weight component.
    routing_N : int, optional
        Number of identical linear reservoirs in the Nash cascade used
        for channel routing (shape parameter of the gamma IUH).
        Default is 2.  Set routing_K to enable routing; routing_N is
        counted as a free parameter in k only when it is explicitly
        calibrated (the caller must add 1 to the k count via
        run_and_score if routing_N is varied).
    routing_K : float or None, optional
        Storage time constant [days] of each Nash-cascade reservoir
        (scale parameter).  Mean travel time through the cascade is
        routing_N * routing_K.  When None (default), no routing is
        applied and the model output is compared directly to observed
        discharge.  When provided, routing_K is counted as one free
        parameter.
    modules : dict or None, optional
        Enable/disable optional process modules.  Keys: ``'snowpack'``,
        ``'frozen_ground'``, ``'rain_on_snow'``, ``'direct_runoff'``.
        Values: bool.  Overrides the ``modules`` block in the YAML config.
        Absent keys leave the YAML/default value unchanged.
    enforce_water_balance : str, optional
        'water-year' (default): scale ET by a per-water-year multiplier so
        that P - Q - ET = 0 over each water year.
        'global': scale ET by a single multiplier computed from
        sum(P - Q_obs) / sum(ET_raw) over the full record — no per-year
        overfitting, does not add hidden degrees of freedom to AIC.
        'none': use raw ET without correction (only for measured ET).
        Legacy boolean True/'False are silently mapped to 'water-year'/'none'.
        Overrides the enforce_water_balance key in the YAML config.

    Returns
    -------
    CalibResult
        Named tuple with fields:
        score        : goodness-of-fit (higher is better)
        aic          : AIC on log flows (lower is better; compare across
                       models with different reservoir counts)
        bfi_obs      : observed baseflow index (Eckhardt filter)
        bfi_mod      : modelled baseflow index (Eckhardt filter)
        fdc_obs      : pd.Series, observed FDC indexed by exceedance %
        fdc_mod      : pd.Series, modelled FDC indexed by exceedance %
        final_states : dict suitable for use as initial_states next decade
        buckets      : Buckets object after the final run

        All scalar fields are np.nan if the scoring window is empty.

    Notes
    -----
    Spin-up runs the full record regardless of start/end so that the
    ET multiplier (computed per water year during initialize()) remains
    valid and the initial storage state reflects long-run climatology.
    """
    if metric not in _METRICS:
        raise ValueError(f"metric must be one of {list(_METRICS)}; got {metric!r}")

    b = Buckets()
    b.initialize(cfg, enforce_water_balance=enforce_water_balance)

    # --- Module flag overrides (must precede parameter overrides that depend on them) ---
    if modules is not None:
        _MATTR = {'snowpack':          'use_snowpack',
                  'frozen_ground':     'use_frozen_ground',
                  'rain_on_snow':      'use_rain_on_snow',
                  'direct_runoff':     'use_direct_runoff',
                  'dtr_fgi_decay':     'use_dtr_fgi_decay',
                  'et_water_stress':   'use_et_water_stress',
                  'et_reservoir_draw': 'use_et_reservoir_draw'}
        for mod, val in modules.items():
            if mod in _MATTR:
                setattr(b, _MATTR[mod], val)
        if not b.use_snowpack:
            b.has_snowpack = False
        if not b.use_dtr_fgi_decay:
            b._has_trange = False
        if b.use_et_water_stress or b.use_et_reservoir_draw:
            # Rebuild ET column in the correct mode (global scale or et_scale,
            # no per-year multiplier). et_scale/et_alpha are still at their
            # defaults here; they will be overridden below if provided.
            b.compute_ET()

    # --- Parameter overrides and free-parameter count ---
    k = 0

    if t_efold is not None:
        for i, val in enumerate(t_efold):
            b.reservoirs[i].t_efold = val
        k += len(t_efold)

    if f_to_discharge is not None:
        for i, val in enumerate(f_to_discharge):
            b.reservoirs[i].f_to_discharge = val
        k += len(f_to_discharge)

    if Hmax is not None:
        for i, val in enumerate(Hmax):
            b.reservoirs[i].Hmax = val
        k += len(Hmax)

    if pdm_H0 is not None:
        for i, val in enumerate(pdm_H0):
            if val is not None:
                b.reservoirs[i].pdm_H0 = val
        k += sum(1 for v in pdm_H0 if v is not None)

    if f_tile is not None:
        any_tile = False
        for i, ft in enumerate(f_tile[:len(b.reservoirs)]):
            b.reservoirs[i].f_tile = ft
            if ft > 0.0 and tau_tile is not None:
                b.reservoirs[i].tile_res = Reservoir(tau_tile, f_to_discharge=1.0)
                any_tile = True
            else:
                b.reservoirs[i].tile_res = None
        k += len({ft for ft in f_tile if ft > 0.0})  # unique values = independent params
        if any_tile:
            k += 1  # tau_tile counted once across all tiled reservoirs

    if et_scale is not None and (b.use_et_water_stress or b.use_et_reservoir_draw):
        if et_scale != 1.0 and b.enforce_water_balance != 'none':
            warnings.warn(
                f"et_scale={et_scale:.4g} with enforce_water_balance="
                f"'{b.enforce_water_balance}' and dynamic ET active: "
                "the ET multiplier is computed assuming et_scale=1, so a "
                "non-unity et_scale amplifies that correction and breaks "
                "the water balance. Set enforce_water_balance='none' when "
                "calibrating et_scale.",
                UserWarning, stacklevel=2,
            )
        b.et_scale = et_scale
        b.compute_ET()   # re-build 'ET for model' column with new scale
        k += 1

    if et_alpha is not None and b.use_et_reservoir_draw:
        b.et_alpha = et_alpha
        k += 1

    if wp_soil is not None and b.use_et_reservoir_draw:
        b.wp_soil = wp_soil
        k += 1

    if wp_soil_sigma is not None and b.use_et_reservoir_draw:
        b.wp_soil_sigma = wp_soil_sigma
        k += 1

    if recession_exponents is not None:
        # H_ref per reservoir: τ is the e-folding time at H = H_ref [mm].
        # Keeps τ bounds interpretation-stable across different exponent values.
        _H_REFS = [50.0, 100.0, 1000.0]
        for i, b_exp in enumerate(recession_exponents):
            if i >= len(b.reservoirs):
                break
            b.reservoirs[i].recession_exponent = float(b_exp)
            b.reservoirs[i].recession_H_ref    = _H_REFS[i] if i < len(_H_REFS) else 100.0
        k += recession_exponents_calibrated

    if melt_factor is not None and b.has_snowpack:
        b.snowpack.melt_factor = melt_factor
        b.melt_factor = melt_factor  # keep Buckets-level attribute in sync
        k += 1

    if fdd_threshold is not None:
        b.fdd_threshold = fdd_threshold
        k += 1

    if snow_insulation_k is not None:
        b.snow_insulation_k = snow_insulation_k
        k += 1

    if direct_runoff_fraction is not None:
        b.direct_runoff_fraction = direct_runoff_fraction
        k += 1

    if baseflow_Q is not None:
        b.baseflow_Q = baseflow_Q
        k += 1

    if routing_K is not None:
        k += 1

    # --- Set initial storage states ---
    if initial_states is not None:
        # Chained decade: restore end-of-previous-decade states exactly.
        for i, h in enumerate(initial_states['reservoirs']):
            b.reservoirs[i].Hwater = h
        if b.has_snowpack:
            b.snowpack.Hwater = initial_states.get('snowpack', 0.0)
        b._fgi = initial_states.get('fgi', 0.0)
    else:
        # Analytical steady-state initialization: correct for reservoirs
        # whose timescale exceeds the record length and accelerates spin-up
        # for all others.
        q_obs  = b.hydrodata['Specific Discharge [mm/day]'].dropna()
        mean_q = float(q_obs.mean())
        if np.isfinite(mean_q) and mean_q > 0:
            mean_q_eff = (mean_q - b.baseflow_Q) * (1.0 - b.direct_runoff_fraction)
            for res, h in zip(b.reservoirs,
                              _steady_state_depths(b.reservoirs, mean_q_eff)):
                res.Hwater = h

    # --- Spin up on the full record with calibrated parameters ---
    if spin_up_cycles is None:
        tau_max = max(r.t_efold for r in b.reservoirs)
        spin_up_cycles = math.ceil(tau_max / len(b.hydrodata))
    for _ in range(spin_up_cycles):
        b.run()

    # --- Final scored run ---
    b.run()

    # --- Capture end-of-run states for chaining to next decade ---
    final_states = {
        'reservoirs': [res.Hwater for res in b.reservoirs],
        'snowpack':   b.snowpack.Hwater if b.has_snowpack else 0.0,
        'fgi':        b._fgi,
    }

    # --- Optional: route total runoff through Nash cascade ---
    # Routing is applied to the full simulation output before the scoring
    # window is applied, so that routing-reservoir state is correct at the
    # window boundaries.  The routed series is written back into the Buckets
    # hydrodata frame so that CalibResult.buckets reflects routed discharge
    # for downstream plotting.
    q_mod = b.hydrodata['Specific Discharge (modeled) [mm/day]']
    if routing_K is not None:
        routed = _nash_cascade(q_mod.to_numpy(), routing_N, routing_K)
        q_mod  = pd.Series(routed, index=q_mod.index, name=q_mod.name)
    # Add constant regional baseflow after routing (external to reservoir cascade).
    if b.baseflow_Q != 0.0:
        q_mod = q_mod + b.baseflow_Q
    b.hydrodata['Specific Discharge (modeled) [mm/day]'] = q_mod

    # --- Mask to scoring window ---
    q_obs = b.hydrodata['Specific Discharge [mm/day]']

    mask = q_mod.notna() & q_obs.notna()
    if start is not None:
        mask &= b.hydrodata['Date'] >= pd.Timestamp(start)
    if end is not None:
        mask &= b.hydrodata['Date'] <= pd.Timestamp(end)

    nan_result = CalibResult(
        score=np.nan, aic=np.nan,
        bfi_obs=np.nan, bfi_mod=np.nan, kge_logfdc=np.nan,
        fdc_obs=pd.Series(dtype=float), fdc_mod=pd.Series(dtype=float),
        final_states=final_states, buckets=b,
    )
    if mask.sum() == 0:
        return nan_result

    m = np.asarray(q_mod[mask], dtype=float)
    o = np.asarray(q_obs[mask], dtype=float)

    return CalibResult(
        score        = _METRICS[metric](m, o),
        aic          = _aic(m, o, k),
        bfi_obs      = _eckhardt_bfi(o),
        bfi_mod      = _eckhardt_bfi(m),
        kge_logfdc   = _kge_logfdc(m, o),
        fdc_obs      = _fdc(o),
        fdc_mod      = _fdc(m),
        final_states = final_states,
        buckets      = b,
    )
