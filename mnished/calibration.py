"""
mnished.calibration
~~~~~~~~~~~~~~~~~~~~~~~
Run MNiShed with a given parameter set and return a CalibResult
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

import copy
import math
import re
import warnings
from collections import namedtuple

import numpy as np
import pandas as pd
import yaml

from .mnished import Buckets, Reservoir

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
    End-of-run storage states suitable for passing as ``initial_states`` to
    the next call. For a single sub-catchment, flat and scalar::

        {'reservoirs': [H_shallow, H_deep, ...],  # [mm]
         'snowpack':    H_snow_SWE,               # [mm SWE]
         'fgi':         frozen_ground_index,      # [°C·day]
         'H_deficit_carry': carried_deficit}      # [mm]

    For several parallel sub-catchments, nested one level deeper — a
    ``'sub_catchments'`` list with one such dict per sub-catchment (config
    order), each carrying that zone's reservoir depths, snowpack, FGI, and
    carried deficit::

        {'sub_catchments': [
            {'reservoirs': [...], 'snowpack': ..., 'fgi': ...,
             'H_deficit_carry': ...},  # sub-catchment 0
            ...]}

buckets : Buckets
    The :class:`~mnished.Buckets` object after the final run,
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
        Ordered shallowest to deepest; each has .recession_coeff, .f_to_discharge,
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
        H_eq = q_in / (np.exp(1.0 / res.recession_coeff) - 1.0)
        depths.append(min(H_eq, res.Hmax))
        # Tile drainage reduces the recharge reaching the next reservoir.
        q_in *= (1.0 - res.f_to_discharge) * (1.0 - res.f_tile)
    return depths


def _apply_reservoir_overrides(reservoirs, over):
    """
    Apply per-reservoir parameter overrides from a dict to one cascade.

    Used for per-sub-catchment overrides in run_and_score; mirrors the flat
    per-reservoir argument semantics but reads from a dict and auto-counts free
    parameters (the flat path uses explicit ``*_calibrated`` counts instead).

    Parameters
    ----------
    reservoirs : list of Reservoir
        The cascade to mutate, shallowest first.
    over : dict
        Override values keyed by parameter name (``'recession_coeff'``,
        ``'f_to_discharge'``, ``'Hmax'``, ``'multipath_threshold'``,
        ``'multipath_timescale'``, ``'leakance_R'``, ``'H_threshold'``,
        ``'recession_exponents'``, ``'f_tile'``, ``'tau_tile'``, ``'pdm_H0'``).

    Returns
    -------
    int
        Number of free parameters set, for AIC k-counting.
    """
    k = 0

    rc = over.get('recession_coeff')
    if rc is not None:
        for i, val in enumerate(rc):
            reservoirs[i].recession_coeff = val
        k += len(rc)

    fd = over.get('f_to_discharge')
    if fd is not None:
        for i, val in enumerate(fd):
            if val is not None:
                reservoirs[i].f_to_discharge = val
        k += sum(1 for v in fd if v is not None)

    lk = over.get('leakance_R')
    if lk is not None:
        for i, val in enumerate(lk):
            if val is not None:
                reservoirs[i].leakance_R = val
                reservoirs[i].junction_type = 'leakance'
        k += sum(1 for v in lk if v is not None)

    ht = over.get('H_threshold')
    if ht is not None:
        for i, val in enumerate(ht):
            if val is not None:
                reservoirs[i].H_threshold = val
                reservoirs[i].junction_type = 'threshold'
        k += sum(1 for v in ht if v is not None)

    mthr = over.get('multipath_threshold')
    mtau = over.get('multipath_timescale')
    if mthr is not None or mtau is not None:
        if mthr is None or mtau is None:
            raise ValueError(
                "multipath_threshold and multipath_timescale must be "
                "provided together (or both omitted).")
        if len(mthr) != len(mtau):
            raise ValueError(
                "multipath_threshold and multipath_timescale must have "
                "the same length.")
        for i, (thr, tau) in enumerate(zip(mthr, mtau)):
            if i >= len(reservoirs):
                break
            if (thr is None) ^ (tau is None):
                raise ValueError(
                    f"Reservoir {i}: multipath_threshold and "
                    "multipath_timescale must be both None or both set.")
            reservoirs[i].multipath_threshold = thr
            reservoirs[i].multipath_timescale = tau
        k += sum(1 for v in mthr if v is not None)

    hm = over.get('Hmax')
    if hm is not None:
        for i, val in enumerate(hm):
            reservoirs[i].Hmax = val
        k += sum(1 for v in hm if np.isfinite(v))

    pd_h0 = over.get('pdm_H0')
    if pd_h0 is not None:
        for i, val in enumerate(pd_h0):
            if val is not None:
                reservoirs[i].pdm_H0 = val
        k += sum(1 for v in pd_h0 if v is not None)

    ft = over.get('f_tile')
    tt = over.get('tau_tile')
    if ft is not None:
        any_tile = False
        for i, ftv in enumerate(ft[:len(reservoirs)]):
            reservoirs[i].f_tile = ftv
            if ftv > 0.0 and tt is not None:
                reservoirs[i].tile_res = Reservoir(tt, f_to_discharge=1.0)
                any_tile = True
            else:
                reservoirs[i].tile_res = None
        k += len({v for v in ft if v > 0.0})
        if any_tile:
            k += 1

    re = over.get('recession_exponents')
    if re is not None:
        for i, b_exp in enumerate(re):
            if i >= len(reservoirs):
                break
            reservoirs[i].recession_exponent = float(b_exp)
        k += len(re)

    return k


def _capture_states(b):
    """
    Capture end-of-run storage state for chaining to a following run.

    For a single sub-catchment the state is flat and scalar
    (``{'reservoirs': [...], 'snowpack': float, 'fgi': float,
    'H_deficit_carry': float}``), identical to earlier versions. For several
    sub-catchments it is nested one level deeper — a ``'sub_catchments'`` list
    with one such dict per sub-catchment (in config order), each carrying that
    zone's reservoir depths plus its own snowpack, FGI, and carried deficit.
    """
    def _one(sc):
        return {
            'reservoirs':      [r.Hwater for r in sc.reservoirs],
            'snowpack':        (sc.snowpack.Hwater
                                if (b.has_snowpack and sc.snowpack is not None)
                                else 0.0),
            'fgi':             sc.fgi,
            'H_deficit_carry': sc.H_deficit_carry,
        }

    if b.n_sub_catchments == 1:
        return _one(b.sub_catchments[0])
    return {'sub_catchments': [_one(sc) for sc in b.sub_catchments]}


def _restore_initial_states(b, states):
    """
    Restore full initial storage state from a (possibly nested) state dict.

    Accepts either the flat single-sub-catchment form or the nested
    ``{'sub_catchments': [...]}`` form produced by :func:`_capture_states`;
    missing snowpack/FGI/deficit entries default to zero (a full restore).
    """
    if 'sub_catchments' in states:
        for sc, s in zip(b.sub_catchments, states['sub_catchments']):
            for r, h in zip(sc.reservoirs, s['reservoirs']):
                r.Hwater = h
            if b.has_snowpack and sc.snowpack is not None:
                sc.snowpack.Hwater = s.get('snowpack', 0.0)
            sc.fgi             = s.get('fgi', 0.0)
            sc.H_deficit_carry = s.get('H_deficit_carry', 0.0)
    else:
        for i, h in enumerate(states['reservoirs']):
            b.reservoirs[i].Hwater = h
        if b.has_snowpack:
            b.snowpack.Hwater = states.get('snowpack', 0.0)
        b._fgi             = states.get('fgi', 0.0)
        b.H_deficit_carry  = states.get('H_deficit_carry', 0.0)


def _inject_post_spinup_states(b, states):
    """
    Inject post-spin-up states, leaving anything not provided at its current
    (spun-up) value.

    Reservoir entries may be ``None`` to keep that reservoir's spun-up depth;
    absent ``snowpack``/``fgi`` keep their current values; ``H_deficit_carry``
    defaults to 0. Accepts the flat or nested form, matching
    :func:`_capture_states`.
    """
    if 'sub_catchments' in states:
        for sc, s in zip(b.sub_catchments, states['sub_catchments']):
            for r, h in zip(sc.reservoirs, s.get('reservoirs', [])):
                if h is not None:
                    r.Hwater = h
            if b.has_snowpack and sc.snowpack is not None:
                sc.snowpack.Hwater = s.get('snowpack', sc.snowpack.Hwater)
            sc.fgi             = s.get('fgi', sc.fgi)
            sc.H_deficit_carry = s.get('H_deficit_carry', 0.0)
    else:
        for i, h in enumerate(states.get('reservoirs', [])):
            if h is not None and i < len(b.reservoirs):
                b.reservoirs[i].Hwater = h
        if b.has_snowpack:
            b.snowpack.Hwater = states.get('snowpack', b.snowpack.Hwater)
        b._fgi             = states.get('fgi', b._fgi)
        b.H_deficit_carry  = states.get('H_deficit_carry', 0.0)


def _validate_finite_states(states, arg_name):
    """
    Raise ValueError if a chained state dict carries non-finite (NaN/inf)
    values.

    Chaining initial conditions from a partial-data or failed decade can yield
    NaN end-states; passed on unchecked they propagate silently through the run
    (every modelled flow becomes NaN and the score looks merely poor rather
    than broken). Surfacing it at the boundary turns silent corruption into an
    immediate, actionable error.

    Accepts the flat single-sub-catchment form or the nested
    ``{'sub_catchments': [...]}`` form. ``None`` reservoir entries (used by
    ``post_spinup_states`` to mean "keep the spin-up value") are skipped.
    """
    if states is None:
        return

    _hint = ("Chained initial states from a partial-data or failed decade can "
             "be NaN; do not chain from it — use initial_states=None for "
             "analytical steady-state initialisation, or skip that decade.")

    def _check(d, where):
        for i, h in enumerate(d.get('reservoirs', [])):
            if h is not None and not math.isfinite(h):
                raise ValueError(
                    f"{where}['reservoirs'][{i}] is non-finite ({h}). {_hint}")
        for key in ('snowpack', 'fgi', 'H_deficit_carry'):
            v = d.get(key)
            if v is not None and not math.isfinite(v):
                raise ValueError(
                    f"{where}['{key}'] is non-finite ({v}). {_hint}")

    if 'sub_catchments' in states:
        for k, sc in enumerate(states['sub_catchments']):
            _check(sc, f"{arg_name}['sub_catchments'][{k}]")
    else:
        _check(states, arg_name)


def _warn_if_flat_states(states, arg_name):
    """
    Emit a DeprecationWarning when a *flat* (single-sub-catchment) state dict is
    passed to run_and_score.

    The flat form (``{'reservoirs': [...], 'snowpack': ..., 'fgi': ...}``) is
    the backward-compatible single-sub-catchment shape; it will be removed in
    v4.0 in favour of the uniform per-sub-catchment form
    (``{'sub_catchments': [...]}``) — see MNiMORPH/MNiShed#18. Still accepted
    for now; the nested form is not deprecated and does not warn.
    """
    if states is not None and 'sub_catchments' not in states:
        warnings.warn(
            f"Passing a flat `{arg_name}` "
            f"({{'reservoirs': [...], 'snowpack': ..., 'fgi': ...}}) is "
            f"deprecated and will be removed in v4.0 (MNiMORPH/MNiShed#18); "
            f"use the per-sub-catchment form "
            f"{{'sub_catchments': [{{'reservoirs': [...], 'snowpack': ..., "
            f"'fgi': ..., 'H_deficit_carry': ...}}]}} instead. (run_and_score "
            f"still returns the flat form at K=1 in v3.x; it becomes nested in "
            f"v4.0.)",
            DeprecationWarning, stacklevel=3)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_and_score(cfg, recession_coeff=None, f_to_discharge=None, Hmax=None,
                  pdm_H0=None, f_tile=None, tau_tile=None,
                  multipath_threshold=None, multipath_timescale=None,
                  multipath_calibrated=0,
                  leakance_R=None, leakance_R_calibrated=0,
                  H_threshold=None, H_threshold_calibrated=0,
                  melt_factor=None, fdd_threshold=None, snow_insulation_k=None,
                  et_scale=None, et_alpha=None,
                  wp_soil=None, wp_soil_sigma=None,
                  recession_exponents=None, recession_exponents_calibrated=0,
                  direct_runoff_fraction=None, baseflow_Q=None,
                  modules=None,
                  sub_catchments=None,
                  initial_states=None,
                  post_spinup_states=None, post_spinup_k=0,
                  start=None, end=None, spin_up_cycles=None,
                  metric='KGE', routing_N=2, routing_K=None,
                  enforce_water_balance='water-year', store_fluxes=False,
                  _model=None):
    """
    Run MNiShed and return a CalibResult named tuple.

    Parameters
    ----------
    cfg : str
        Path to a YAML configuration file. Should have spin_up_cycles: 0
        so that spin-up is performed here with the calibrated parameters.
    recession_coeff : list of float, optional
        Recession coefficients [days], one per reservoir (shallowest
        first). Overrides the values in cfg.
    f_to_discharge : list of float, optional
        Exfiltration fractions to stream, one per reservoir except the
        deepest (which is always 1.0). Overrides the values in cfg.
        Used only for reservoirs whose junction_type is 'fraction' or 'threshold'.
    leakance_R : list of float or None, optional
        Leakance resistance [days] for each reservoir, one per reservoir.
        ``None`` entries leave the reservoir at its config-defined junction type.
        Non-None entries set junction_type to 'leakance' and assign the resistance.
        ``Q_leak = max(H_this - H_next, 0) / R``.  Default None (all fraction).
    leakance_R_calibrated : int, optional
        Number of entries in ``leakance_R`` that are free calibration
        parameters (for AIC k-counting).  Default 0.
    H_threshold : list of float or None, optional
        Dead-storage threshold depth [mm] for each reservoir, one per reservoir.
        ``None`` entries leave the reservoir unchanged.  Non-None entries set
        junction_type to 'threshold': only ``max(H - H_threshold, 0)`` drains.
        Models a stream-aquifer connection that activates above a threshold head.
        Default None (no threshold).
    H_threshold_calibrated : int, optional
        Number of entries in ``H_threshold`` that are free calibration
        parameters (for AIC k-counting).  Default 0.
    multipath_threshold : list of float or None, optional
        Storage depth [mm] above which a parallel fast drain activates,
        per reservoir.  ``None`` entries leave the reservoir's multipath
        configuration unchanged from cfg.  Non-None entries enable a
        threshold-activated parallel drain that adds
        ``max(0, H - thr)/multipath_timescale`` to discharge above the
        threshold.  Distinct from ``f_tile``/``tau_tile`` (which is a
        constant-fraction bypass through a downstream sub-reservoir);
        see :class:`mnished.Reservoir` for the contrast.
        Requires the matching ``multipath_timescale`` entry to be non-None.
        Default None.
    multipath_timescale : list of float or None, optional
        E-folding timescale [days] of the parallel multipath drain, per
        reservoir.  Required paired with ``multipath_threshold``.
        Default None.
    multipath_calibrated : int, optional
        Number of free parameters contributed by multipath_threshold +
        multipath_timescale combined (for AIC k-counting).  Default 0.
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
        (Molnau & Bissell 1983,
        https://westernsnowconference.org/bibliography/1983Molnau.pdf).  When
        the index exceeds fdd_threshold, infiltration from the top reservoir to
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
        Starting reservoir water depths, snowpack SWE, frozen ground index,
        and carried deficit, as returned by a previous call's
        CalibResult.final_states.  Flat for a single sub-catchment::

            {'reservoirs': [H_shallow, H_deep, ...],
             'snowpack':    H_snow_SWE,
             'fgi':         frozen_ground_index,
             'H_deficit_carry': carried_deficit}

        or nested under ``'sub_catchments'`` for several (see
        ``CalibResult.final_states``).  When provided, these override the H0
        values from cfg.  Use with spin_up_cycles=0 when chaining consecutive
        decades so that water storage and frozen-ground state are physically
        continuous.  When None (default), reservoirs are initialised to their
        analytical steady-state depths before spin-up (see
        :func:`_steady_state_depths`).
    post_spinup_states : dict, optional
        Reservoir water depths injected *after* spin-up completes and
        *before* the scored run begins.  Same format as ``initial_states``
        (flat, or nested under ``'sub_catchments'``); individual reservoir
        entries may be ``None`` to leave that reservoir at its spin-up end
        state, and absent snowpack/FGI keep their spun-up values.
        Intended for calibrating decade-specific initial storage (e.g.
        ``log__H0_deep``) when the spin-up equilibrium is poorly
        constrained by sparse pre-decade forcing.  Only applied in decade
        mode (``start`` is not None); ignored in full-record mode.
    post_spinup_k : int, optional
        Number of free parameters contributing to ``post_spinup_states``
        (for AIC counting).  Default 0.
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
    metric : {'KGE', 'NSE', 'logKGE', 'KGE_logKGE', 'KGE_logKGE_logFDC',
              'KGE_logKGE_logFDC_BFI'}, optional
        Goodness-of-fit metric.  Default is 'KGE'.
        ``'KGE_logKGE'`` returns 0.5*KGE + 0.5*logKGE, balancing
        peak and low-flow performance (Yilmaz et al. 2008).
        ``'KGE_logKGE_logFDC_BFI'`` adds a BFI bias-ratio score
        (1 - abs(BFI_mod/BFI_obs - 1)) as a fourth equal-weight component.
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
    sub_catchments : list of dict or None, optional
        Per-sub-catchment parameter overrides for a basin partitioned into
        parallel sub-catchments (see :class:`mnished.SubCatchment`). The YAML
        config defines the *structure* (how many sub-catchments and reservoirs);
        this argument overrides values by position, one dict per sub-catchment
        in config order. Each dict may contain ``'area_fraction'`` and any of
        the per-reservoir override lists (``'recession_coeff'``,
        ``'f_to_discharge'``, ``'Hmax'``, ``'multipath_threshold'``,
        ``'multipath_timescale'``, ``'leakance_R'``, ``'H_threshold'``,
        ``'recession_exponents'``, ``'f_tile'``, ``'tau_tile'``, ``'pdm_H0'``),
        scoped to that sub-catchment's cascade with the same semantics as the
        flat arguments. A lake sub-catchment's dict may also contain
        ``'f_route_lake'`` (the channelized-routing fraction in ``[0, 1]``;
        requires the lake's ``gw_partner``), letting it be calibrated directly
        without rewriting the config; the land zone's routed-away fraction is
        kept consistent automatically. Mutually exclusive with the flat per-reservoir
        arguments. When ``'area_fraction'`` is given it must still sum to 1
        across sub-catchments. For AIC, every overridden value counts as one
        free parameter, plus n_sub_catchments - 1 when area fractions are set.
        Snow and ET parameters remain basin-level (the flat ``melt_factor``,
        ``et_scale``, etc.). Default None (single-cascade calibration).
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
        Named tuple with fields ``score``, ``aic``, ``bfi_obs``,
        ``bfi_mod``, ``fdc_obs``, ``fdc_mod``, ``final_states``,
        and ``buckets``.  See :class:`CalibResult` for full field
        descriptions.  All scalar fields are ``np.nan`` if the scoring
        window contains no valid overlapping data.

    Notes
    -----
    When ``start`` is provided, the function operates in decade mode:

    * ``enforce_water_balance='global'`` closes the water balance over
      the full record (computed once at initialization), not re-fitted to
      the decade window.  Re-fitting would assume ΔS = 0 over the decade,
      which is incorrect for transitional wet/dry decades.
    * Spin-up runs only the pre-decade data (record start through the day
      before ``start``), so reservoir states at the beginning of the decade
      reflect pre-decade climatology rather than an arbitrary end-of-record
      state.
    * The scored run covers ``[start, end]`` only, continuing directly from
      the spin-up end state.

    When ``start`` is None (default), the original full-record behaviour is
    preserved: spin-up and scoring both run the complete hydrodata record.
    """
    if metric not in _METRICS:
        raise ValueError(f"metric must be one of {list(_METRICS)}; got {metric!r}")

    # Deprecation: the flat single-sub-catchment state shape is removed in v4.0
    # (#18). Warn on flat input; the nested form is the forward contract.
    _warn_if_flat_states(initial_states, 'initial_states')
    _warn_if_flat_states(post_spinup_states, 'post_spinup_states')
    # Fail loudly on NaN/inf chained states rather than silently propagating
    # them through the run (a partial-data/failed decade can produce them).
    _validate_finite_states(initial_states, 'initial_states')
    _validate_finite_states(post_spinup_states, 'post_spinup_states')

    if _model is not None:
        # Reuse a pre-built model (forcing read, construction, and base ET
        # already done once): a fresh copy of the template, no CSV re-read.
        b = _model._fresh_buckets()
    else:
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

    if sub_catchments is not None:
        # Per-sub-catchment overrides. The flat per-reservoir arguments are
        # mutually exclusive with this; basin-level (snow/ET/routing) overrides
        # below still apply. The cfg owns the structure (number of
        # sub-catchments and reservoirs); this overrides values by position.
        _flat = {
            'recession_coeff': recession_coeff, 'f_to_discharge': f_to_discharge,
            'Hmax': Hmax, 'pdm_H0': pdm_H0, 'f_tile': f_tile, 'tau_tile': tau_tile,
            'multipath_threshold': multipath_threshold,
            'multipath_timescale': multipath_timescale,
            'leakance_R': leakance_R, 'H_threshold': H_threshold,
            'recession_exponents': recession_exponents,
        }
        _bad = [name for name, val in _flat.items() if val is not None]
        if _bad:
            raise ValueError(
                "When 'sub_catchments' is given, per-reservoir parameters must "
                "be supplied inside each sub-catchment dict, not as flat "
                f"arguments; remove: {_bad}.")
        if len(sub_catchments) != b.n_sub_catchments:
            raise ValueError(
                f"'sub_catchments' has {len(sub_catchments)} entries but the "
                f"config defines {b.n_sub_catchments} sub-catchments.")
        _set_area = False
        _routed_overridden = False
        for _sc_obj, _sc_over in zip(b.sub_catchments, sub_catchments):
            if 'area_fraction' in _sc_over:
                _sc_obj.area_fraction = _sc_over['area_fraction']
                _set_area = True
            if 'f_route_lake' in _sc_over:
                if _sc_obj.kind != 'lake':
                    raise ValueError(
                        f"Sub-catchment '{_sc_obj.name}': 'f_route_lake' "
                        "override only applies to a lake sub-catchment.")
                _fr = float(_sc_over['f_route_lake'])
                if not 0.0 <= _fr <= 1.0:
                    raise ValueError(
                        f"Lake '{_sc_obj.name}': f_route_lake override must be "
                        f"in [0, 1]; got {_fr}.")
                if _fr > 0.0 and _sc_obj.gw_partner is None:
                    raise ValueError(
                        f"Lake '{_sc_obj.name}': f_route_lake > 0 needs a "
                        "routing source (the lake's 'gw_partner').")
                _sc_obj.f_route_lake = _fr
                _routed_overridden = True
                k += 1
            k += _apply_reservoir_overrides(_sc_obj.reservoirs, _sc_over)
        if _routed_overridden:
            # Keep each land zone's routed-away fraction in step with the
            # overridden f_route_lake (the parser does this at config time).
            b._resolve_routed_away_fractions()
        if _set_area:
            _asum = sum(s.area_fraction for s in b.sub_catchments)
            if abs(_asum - 1.0) > 1e-6:
                raise ValueError(
                    f"Sub-catchment area_fraction values must sum to 1.0; "
                    f"got {_asum:.6f}.")
            k += b.n_sub_catchments - 1   # last fraction is determined

    if recession_coeff is not None:
        for i, val in enumerate(recession_coeff):
            b.reservoirs[i].recession_coeff = val
        k += len(recession_coeff)

    if f_to_discharge is not None:
        for i, val in enumerate(f_to_discharge):
            if val is not None:
                b.reservoirs[i].f_to_discharge = val
        k += sum(1 for v in f_to_discharge if v is not None)

    if leakance_R is not None:
        for i, val in enumerate(leakance_R):
            if val is not None:
                b.reservoirs[i].leakance_R = val
                b.reservoirs[i].junction_type = 'leakance'
        k += leakance_R_calibrated

    if H_threshold is not None:
        for i, val in enumerate(H_threshold):
            if val is not None:
                b.reservoirs[i].H_threshold = val
                b.reservoirs[i].junction_type = 'threshold'
        k += H_threshold_calibrated

    if multipath_threshold is not None or multipath_timescale is not None:
        # Both lists must be supplied if either is, and must be the same length.
        if multipath_threshold is None or multipath_timescale is None:
            raise ValueError(
                "multipath_threshold and multipath_timescale must be "
                "provided together (or both omitted).")
        if len(multipath_threshold) != len(multipath_timescale):
            raise ValueError(
                "multipath_threshold and multipath_timescale must have "
                "the same length.")
        for i, (thr, tau) in enumerate(zip(multipath_threshold,
                                            multipath_timescale)):
            if i >= len(b.reservoirs):
                break
            if (thr is None) ^ (tau is None):
                raise ValueError(
                    f"Reservoir {i}: multipath_threshold and "
                    "multipath_timescale must be both None or both set.")
            b.reservoirs[i].multipath_threshold = thr
            b.reservoirs[i].multipath_timescale = tau
        k += multipath_calibrated

    if Hmax is not None:
        for i, val in enumerate(Hmax):
            b.reservoirs[i].Hmax = val
        # Count only finite caps as free parameters; an .inf entry means
        # "no saturation-excess cap" (a structural choice, not a calibrated
        # value), consistent with the non-None counts for f_to_discharge
        # and pdm_H0.
        k += sum(1 for v in Hmax if np.isfinite(v))

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

    if et_scale is not None:
        if et_scale != 1.0 and b.enforce_water_balance != 'none':
            warnings.warn(
                f"et_scale={et_scale:.4g} with enforce_water_balance="
                f"'{b.enforce_water_balance}': et_scale is applied on top of "
                "the water-balance multiplier, so exact WB closure is not "
                "guaranteed. Set enforce_water_balance='none' to use et_scale "
                "as the primary water-balance parameter.",
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
        # recession_H_ref stays at the Reservoir default of 1.0, so the
        # recession coefficient is the raw drainage constant.  A non-unit
        # H_ref is a redundant gauge: only the product recession_coeff *
        # H_ref^(b-1) is identifiable, so it just rescales the coefficient
        # (see Reservoir.discharge / mean_residence_time).
        for i, b_exp in enumerate(recession_exponents):
            if i >= len(b.reservoirs):
                break
            b.reservoirs[i].recession_exponent = float(b_exp)
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

    k += post_spinup_k

    # --- Set initial storage states ---
    if initial_states is not None:
        # Analytical or chained initial conditions supplied by the caller.
        # H_deficit_carry can accumulate a large phantom deficit during
        # b.initialize()'s internal spin-up (especially with
        # enforce_water_balance='none'); _restore_initial_states resets it
        # (default 0) so the decade run starts cleanly from the given depths.
        _restore_initial_states(b, initial_states)
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

    # Decade mode: the global ET multiplier is computed once from the full
    # record during b.initialize() above.  It is NOT recomputed for the
    # decade window here; doing so would assume ΔS = 0 over the decade,
    # which is wrong for transitional decades (e.g. wet→Dust Bowl) where
    # reservoirs are net draining and Q_obs is elevated by that drainage.

    # --- Spin up, then final scored run ---
    if spin_up_cycles is None:
        tau_max = max(r.recession_coeff for r in b.reservoirs)
        spin_up_cycles = math.ceil(tau_max / len(b.hydrodata))

    if start is not None:
        # Decade mode: spin up on pre-decade data only, then run the decade.
        pre_decade_end = (pd.Timestamp(start)
                          - pd.Timedelta(days=1)).strftime('%Y-%m-%d')
        # H_deficit_carry may have accumulated a large phantom value during
        # initialize()'s internal spin-up (which runs without et_reservoir_draw).
        # Reset it before the pre-decade spin-up so the spin-up is physically
        # clean and its end states are usable as decade initial conditions.
        b.H_deficit_carry = 0.0
        for _ in range(spin_up_cycles):
            b.run(end=pre_decade_end)
        # Inject post-spin-up reservoir depths if provided (e.g. log__H0_deep).
        # Allows calibrating decade-specific initial storage independently of
        # the spin-up equilibrium, which may be biased by sparse pre-decade
        # forcing or by candidate parameters far from the decade optimum.
        if post_spinup_states is not None:
            _inject_post_spinup_states(b, post_spinup_states)
        b.run(start=start, end=end, store_fluxes=store_fluxes)
    else:
        # Full-record mode: spin up and score on the complete hydrodata.
        # initialize() runs an internal spin-up that can leave H_deficit_carry at
        # a large phantom value; reset so the calibration spin-up starts clean.
        b.H_deficit_carry = 0.0
        for _ in range(spin_up_cycles):
            b.run()
        b.run(store_fluxes=store_fluxes)

    # --- Capture end-of-run states for chaining to next decade ---
    # Flat/scalar for a single sub-catchment (back-compatible); nested per
    # sub-catchment when there are several.
    final_states = _capture_states(b)

    # --- Optional: route total runoff through Nash cascade ---
    # Routing is applied to the full simulation output before the scoring
    # window is applied, so that routing-reservoir state is correct at the
    # window boundaries.  The routed series is written back into the Buckets
    # hydrodata frame so that CalibResult.buckets reflects routed discharge
    # for downstream plotting.
    # pd.to_numeric converts pd.NA (unrun rows in decade mode) to np.nan so
    # that the numpy-based Nash cascade and downstream code handle them cleanly.
    q_mod = pd.to_numeric(
        b.hydrodata['Specific Discharge (modeled) [mm/day]'], errors='coerce')
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


def log_flow_residual_terms(result, start=None, end=None, eps=None):
    """Per-day scored log-flow terms for Bayesian calibration.

    Reproduces :func:`run_and_score`'s scoring mask (observed *and*
    modelled discharge both present, within ``[start, end]``) on a
    :class:`CalibResult`'s ``buckets``, and returns the log-transformed
    observed and modelled flows whose difference is the Gaussian log-flow
    residual that the AIC — and a Dakota ``bayes_calibration`` — use.

    The mask is **independent of the calibrated parameters** (it depends
    only on observation availability), so across a calibration the day set
    and the ``log_obs`` column are fixed and only ``log_mod`` changes.
    That lets a Bayesian driver return ``log_mod`` as ``calibration_terms``
    against a once-computed ``log_obs`` ``calibration_data_file``.

    Parameters
    ----------
    result : CalibResult
        A result whose ``.buckets`` holds the run's ``hydrodata`` frame
        (i.e. the modelled discharge has been written back, as
        :func:`run_and_score` does before returning).
    start, end : str or pandas.Timestamp, optional
        Scoring-window bounds; pass the same values given to
        :func:`run_and_score` so the day set matches the scored window.
    eps : float, optional
        Stabilising offset added before the log.  Defaults to 1 % of the
        mean scored observed flow, matching :func:`run_and_score`.

    Returns
    -------
    pandas.DataFrame
        One row per scored day, with columns ``date``, ``obs``, ``mod``,
        ``log_obs``, ``log_mod`` and ``residual`` (= ``log_mod -
        log_obs``).
    """
    h = result.buckets.hydrodata
    q_obs = pd.to_numeric(h['Specific Discharge [mm/day]'], errors='coerce')
    q_mod = pd.to_numeric(
        h['Specific Discharge (modeled) [mm/day]'], errors='coerce')
    mask = q_mod.notna() & q_obs.notna()
    if start is not None:
        mask &= h['Date'] >= pd.Timestamp(start)
    if end is not None:
        mask &= h['Date'] <= pd.Timestamp(end)
    o = q_obs[mask].to_numpy(dtype=float)
    m = q_mod[mask].to_numpy(dtype=float)
    if eps is None:
        eps = 0.01 * o.mean()
    log_o = np.log(o + eps)
    log_m = np.log(m + eps)
    return pd.DataFrame({
        'date':     h['Date'][mask].to_numpy(),
        'obs':      o,
        'mod':      m,
        'log_obs':  log_o,
        'log_mod':  log_m,
        'residual': log_m - log_o,
    })


class ScoringModel:
    """Build-once, score-many wrapper around :func:`run_and_score`.

    Building a model from a config reads the forcing CSV, constructs the
    sub-catchment/reservoir cascade, and computes the base ET — work that is
    independent of the calibrated parameters but that :func:`run_and_score`
    otherwise repeats on every call. For an in-process optimiser or sampler
    (thousands of evaluations) that rebuild dominates the wall-clock, and a
    multi-window objective re-reads the same CSV once per window per eval.

    ``ScoringModel`` does the build once; :meth:`score` then reuses a fresh
    copy of the template for each evaluation (a deep copy of the built model is
    ~free relative to a re-initialise) and is **bit-identical** to the
    equivalent :func:`run_and_score` call — that equivalence is the contract.

    Parameters
    ----------
    cfg : str
        Path to the YAML config, as for :func:`run_and_score`.
    enforce_water_balance : str, optional
        Water-balance closure mode. Fixed for the model's lifetime (it is a
        structural choice, not a calibrated parameter). Default 'water-year'.

    Examples
    --------
    >>> sm = ScoringModel('config.yml', enforce_water_balance='none')
    >>> r = sm.score(recession_coeff=[14, 500], et_scale=0.8,
    ...              start='2001-01-01', end='2010-12-31', metric='KGE')
    >>> r.score                       # identical to the run_and_score call
    """

    def __init__(self, cfg, enforce_water_balance='water-year'):
        self.cfg = cfg
        self.enforce_water_balance = enforce_water_balance
        self._template = Buckets()
        self._template.initialize(
            cfg, enforce_water_balance=enforce_water_balance)

    def _fresh_buckets(self):
        """A clean copy of the built template (no CSV re-read / rebuild)."""
        return copy.deepcopy(self._template)

    def score(self, **kwargs):
        """Score one parameter set, reusing the pre-built model.

        Takes the same keyword arguments as :func:`run_and_score` other than
        ``cfg`` and ``enforce_water_balance`` (fixed at construction) and
        returns the same :class:`CalibResult`.
        """
        return run_and_score(self.cfg,
                             enforce_water_balance=self.enforce_water_balance,
                             _model=self, **kwargs)


def target_kwargs(parameters, theta, n_sub=1):
    """Map declarative parameter ``target``\\ s to :func:`run_and_score` keywords.

    ``parameters`` is a config ``parameters`` block: a dict of
    ``{name: {target, fixed, ...}}``. ``theta`` gives the current value of each
    free parameter by name (a parameter absent from ``theta`` uses its
    ``fixed`` value). Each parameter's ``target`` declares where its value goes:

    * ``name``                    -> a scalar keyword
    * ``name[i]``                 -> a list keyword at position i
    * ``sub_catchments[I].key``   -> a scalar override on sub-catchment(s) I
    * ``sub_catchments[I].key[j]``-> a list override, position j

    ``I`` is an index or comma-list (``0,1`` for a parameter shared across
    zones). A ``log__`` name prefix applies ``10**`` to the value. Untargeted
    list positions are left ``None`` (run_and_score keeps the config value), so
    a parameter need only name the elements it calibrates. ``n_sub`` is the
    number of sub-catchments (sizes the positional ``sub_catchments`` override).

    Parameters without a ``target`` are ignored, so non-model entries (e.g. a
    sampled error-model nuisance parameter) can coexist in the config.
    """
    flat_lists, kw = {}, {}
    sub_lists, sub_scalar = {}, {}
    for name, spec in parameters.items():
        target = spec.get('target')
        if not target:
            continue
        val = theta.get(name, spec.get('fixed'))
        if name.startswith('log__'):
            val = 10.0 ** val
        nested = re.fullmatch(r'sub_catchments\[([\d,]+)\]\.(.+)', target)
        if nested:
            rest = re.fullmatch(r'(\w+)\[(\d+)\]', nested.group(2))
            for i in (int(x) for x in nested.group(1).split(',')):
                if rest:
                    sub_lists.setdefault(
                        (i, rest.group(1)), {})[int(rest.group(2))] = val
                else:
                    sub_scalar[(i, nested.group(2))] = val
            continue
        flat = re.fullmatch(r'(\w+)\[(\d+)\]', target)
        if flat:
            flat_lists.setdefault(flat.group(1), {})[int(flat.group(2))] = val
        else:
            kw[target] = val
    for key, idx in flat_lists.items():
        kw[key] = [idx.get(i) for i in range(max(idx) + 1)]
    if sub_lists or sub_scalar:
        sub = [{} for _ in range(n_sub)]
        for (i, key), idx in sub_lists.items():
            sub[i][key] = [idx.get(j) for j in range(max(idx) + 1)]
        for (i, key), v in sub_scalar.items():
            sub[i][key] = v
        kw['sub_catchments'] = sub
    return kw


class Calibrator:
    """A declarative, build-once calibration problem read from a config.

    This is MNiShed's standard model-setup for calibration: the *run method is
    defined in config, not code*. A config (the ``parameters``, ``driver``, and
    ``modules`` blocks of a ``params.yml``) names the free parameters, their
    bounds, and — via each parameter's ``target`` (see :func:`target_kwargs`) —
    where each maps in the model. ``Calibrator`` ties that together with a
    build-once :class:`ScoringModel`, so any sampler or optimiser needs only to
    call :meth:`score`; there is no per-basin Python.

    It is sampler-agnostic (no optimiser dependency): point SciPy, SPOTPY
    (SCE-UA / DREAM), or any driver at :meth:`score` and
    :attr:`parameter_set` (the free parameters and bounds).

    Parameters
    ----------
    parameters : dict
        The ``parameters`` block (each entry with bounds and a ``target``).
    driver : dict
        Run settings: ``config_template`` (the model YAML), ``metric``,
        ``spin_up_cycles``, ``routing_N``, ``enforce_water_balance``, and the
        calibration window(s). Windows are given as ``decades:`` — a list of
        ``{start, end}`` dicts (multi-window objective; ``None`` = full record).
        ``decade_start`` / ``decade_end`` are a shorthand for a single such
        window. (The name ``decades:`` will become ``windows:`` in v4.0; see
        MNiMORPH/MNiShed#24.)
    modules : dict, optional
        Process-module toggles, applied per evaluation.

    Examples
    --------
    >>> cal = Calibrator.from_yaml('params.yml')
    >>> cal.score({'log__t_recession_shallow': 1.4, ...}).score
    >>> cal.score(vector).score        # vector ordered as cal.names
    """

    def __init__(self, parameters, driver, modules=None):
        from .identifiability import ParameterSet
        self.parameters = parameters
        self.driver = driver
        self.modules = modules or {}
        self.parameter_set = ParameterSet.from_params_yml(parameters)
        self.names = self.parameter_set.names
        self.model = ScoringModel(
            driver['config_template'],
            enforce_water_balance=driver.get('enforce_water_balance',
                                             'water-year'))
        with open(driver['config_template']) as f:
            self.n_sub = len(yaml.safe_load(f).get('sub_catchments', [])) or 1
        # Calibration windows: the one window mechanism. ``decades:`` is a list
        # of ``{start, end}`` dicts; ``decade_start`` / ``decade_end`` is the
        # shorthand for a single window (``None`` = full record). Both land in
        # ``self.windows``, which score() / score_windows() read uniformly.
        self.windows = driver.get('decades') or [
            {'start': driver.get('decade_start'),
             'end':   driver.get('decade_end')}]

    @classmethod
    def from_yaml(cls, path):
        """Build a :class:`Calibrator` from a ``params.yml`` path."""
        with open(path) as f:
            cfg = yaml.safe_load(f)
        return cls(cfg['parameters'], cfg['driver'], cfg.get('modules'))

    def run_kwargs(self, theta):
        """The :func:`run_and_score` keywords for a parameter dict."""
        return target_kwargs(self.parameters, theta, self.n_sub)

    def score(self, theta, start=None, end=None, metric=None):
        """Score one parameter set, reusing the build-once model.

        ``theta`` is a ``{name: value}`` dict or a vector ordered as
        :attr:`names`. ``start`` / ``end`` / ``metric`` default to the driver
        config. Returns a :class:`CalibResult` — identical to the equivalent
        :func:`run_and_score` call.
        """
        if not isinstance(theta, dict):
            theta = dict(zip(self.names, theta))
        d = self.driver
        # Default to the first (or only) calibration window, so a single-element
        # ``decades:`` list and the ``decade_start`` / ``decade_end`` shorthand
        # score the same span through one mechanism (:attr:`windows`).
        w0 = self.windows[0]
        return self.model.score(
            modules=self.modules, routing_N=d['routing_N'],
            spin_up_cycles=d['spin_up_cycles'],
            metric=metric if metric is not None else d['metric'],
            start=start if start is not None else w0.get('start'),
            end=end if end is not None else w0.get('end'),
            **self.run_kwargs(theta))

    def score_windows(self, theta, windows=None, metric=None):
        """Score ``theta`` on each calibration window; a list of CalibResult.

        ``windows`` is a list of ``{'start', 'end'}`` dicts and defaults to
        :attr:`windows` (the driver's ``decades:`` list, or the single
        ``decade_start`` / ``decade_end`` window). Aggregating the per-window
        results — a mean score for optimisation, concatenated residuals for a
        likelihood — is the caller's choice, so this stays sampler-agnostic.
        """
        return [self.score(theta, start=w.get('start'), end=w.get('end'),
                            metric=metric)
                for w in (windows if windows is not None else self.windows)]
