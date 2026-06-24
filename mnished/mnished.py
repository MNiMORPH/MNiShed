#! /usr/bin/python3

########################################
# Then, methought, the air grew denser #
#         - Edgar Allan Poe            #
#              THE RAVEN               #
########################################

# Started by A. Wickert
# 25 July 2019
# Updated slightly by J. Jones
# 08 Oct 2019
# Significant Update by A. Wickert
# October 2022
# CLI added by A. Wickert
# November 2023

import argparse
import math
import sys
import warnings

import matplotlib.dates as mdates
import numpy as np
import pandas as pd
import yaml
from matplotlib import pyplot as plt

try:
    import numba as _numba
    _numba_available = True
    _numba_import_error = None
except ModuleNotFoundError:
    # numba is not installed: pure Python is the expected default (the JIT is
    # the optional `jit` extra), so this case stays silent.
    _numba_available = False
    _numba_import_error = None
except ImportError as _exc:
    # numba is installed but failed to import — commonly the NumPy/Numba
    # version mismatch the `jit` extra guards against (numpy<2.3). The user
    # opted into the JIT but would silently get pure Python; flagged once at
    # the first run() (see _notify_jit_unavailable).
    _numba_available = False
    _numba_import_error = _exc

# c_p / L_f: water's specific heat divided by the latent heat of fusion.
# c_p = 4186 J kg⁻¹ °C⁻¹, L_f = 334 000 J kg⁻¹  →  ≈ 0.01253 °C⁻¹
# Gives mm SWE melted per mm of rain per °C of rain temperature.
_CP_LF = 4186.0 / 334000.0

# Set True once the first run() falls back to the pure-Python loop for a reason
# worth surfacing, so the notice fires at most once per process.
_jit_unavailable_notified = False


def _notify_jit_unavailable(reason):
    """
    Emit a one-time UserWarning that the pure-Python time loop is in use.

    Called from Buckets.run() when the Numba JIT is *not* used for a reason
    worth surfacing — numba installed-but-broken, or a configuration the JIT
    does not support. Fires at most once per process. The plain "numba not
    installed" default does not call this (pure Python is the expected
    behaviour without the ``jit`` extra), and ``reason=None`` is a no-op.
    """
    global _jit_unavailable_notified
    if _jit_unavailable_notified or reason is None:
        return
    _jit_unavailable_notified = True
    warnings.warn(
        "MNiShed is using the pure-Python time loop (~100x slower than the "
        f"Numba JIT): {reason}.",
        UserWarning, stacklevel=3)


if _numba_available:
    @_numba.jit(nopython=True, cache=True)
    def _jit_run(P_arr, ET_arr, T_arr, T_min_arr, T_max_arr,
                 sc_start_idx, sc_end_idx, sc_area_fractions,
                 H_init, H_snow_init, fgi_init, H_deficit_carry_init,
                 tau_arr, b_arr, H_ref_arr, f_dis_in, junction_arr,
                 leakance_R_arr, H_threshold_arr, Hmax_arr,
                 f_tile_arr, tau_tile_arr, H_tile_init,
                 multipath_thr_arr, multipath_tau_arr, pdm_H0_arr,
                 melt_factor, snow_insulation_k, fgi_decay_coeff, fdd_threshold,
                 direct_runoff_fraction, wp_soil, wp_soil_sigma, et_alpha, dt,
                 has_snowpack, use_fgi, use_rain_on_snow, use_et_reservoir_draw,
                 use_direct_runoff, use_et_water_stress, has_trange):
        """JIT-compiled daily time loop replacing Buckets.run()/update().

        Replicates update()/_advance_sub_catchment()/_compute_snowpack()/
        _update_fgi()/_draw_et_from_reservoirs() exactly.  PDM saturation-excess
        (``pdm_H0_arr`` > 0) and et_water_stress (``use_et_water_stress``) are
        both supported; ``pdm_H0_arr[r] <= 0`` disables PDM for reservoir *r*.

        Parallel sub-catchments are flattened into the per-reservoir arrays;
        ``sc_start_idx[k]``/``sc_end_idx[k]`` give the half-open reservoir-index
        range of sub-catchment *k*, and ``sc_area_fractions[k]`` its basin-area
        weight. Snowpack SWE, FGI, and the carried deficit are per-sub-catchment
        (``H_snow_init``/``fgi_init``/``H_deficit_carry_init`` are length-n_sub
        arrays). Basin discharge, SWE, and storage are the area-weighted means
        over sub-catchments. A single sub-catchment (area 1.0) reproduces the
        original behaviour exactly.

        Two tile-drain mechanisms are supported (independently or together):
          - Constant-fraction bypass (``f_tile_arr``/``tau_tile_arr``):
            a fraction of recession outflow is routed through a downstream
            fast linear sub-reservoir, regardless of storage state.
          - Multipath threshold-activated parallel drain
            (``multipath_thr_arr``/``multipath_tau_arr``): a second linear
            outflow path operates from the same storage but only when the
            reservoir depth exceeds the threshold. ``multipath_tau_arr[r] <= 0``
            disables the multipath drain for reservoir *r*.

        junction_arr encoding: 0 = fraction, 1 = leakance, 2 = threshold.
        """
        n_steps = len(P_arr)
        n_res   = len(tau_arr)
        n_sub   = len(sc_start_idx)

        Q_out     = np.empty(n_steps)
        SWE_out   = np.empty(n_steps)
        H_sub_out = np.empty(n_steps)
        H_res_out = np.empty((n_steps, n_res))

        H_res  = H_init.copy()
        H_snow = H_snow_init.copy()             # per-sub-catchment SWE [mm]
        fgi    = fgi_init.copy()                 # per-sub-catchment FGI [°C·day]
        H_dc   = H_deficit_carry_init.copy()     # per-sub-catchment carried deficit

        # Local copy of f_to_discharge — frozen-ground overrides the top
        # reservoir of a sub-catchment temporarily.
        f_dis = f_dis_in.copy()

        # Per-reservoir cascade working arrays (flat across all sub-catchments)
        H_to_next = np.zeros(n_res)
        H_deficit_res = np.zeros(n_res)
        H_tile = H_tile_init.copy()

        for step in range(n_steps):
            P_t  = P_arr[step]
            ET_t = ET_arr[step]
            T_t  = 0.0

            # Read temperature when needed by snowpack or FGI
            if has_snowpack or use_fgi:
                T_t = T_arr[step]

            # Skip step if any required forcing is NaN.
            # T NaN only triggers skip when snowpack is active (matches Python).
            skip = math.isnan(P_t) or math.isnan(ET_t)
            if has_snowpack and not skip:
                skip = math.isnan(T_t)

            if skip:
                Q_out[step] = np.nan
                swe_tot = 0.0
                H_sub_total = 0.0
                for sc in range(n_sub):
                    area = sc_area_fractions[sc]
                    swe_tot += area * H_snow[sc]
                    for r in range(sc_start_idx[sc], sc_end_idx[sc]):
                        H_sub_total += area * (H_res[r] + H_tile[r])
                for r in range(n_res):
                    H_res_out[step, r] = H_res[r]
                SWE_out[step]   = swe_tot
                H_sub_out[step] = H_sub_total
                continue

            qi_basin    = 0.0
            swe_tot     = 0.0
            H_sub_total = 0.0

            for sc in range(n_sub):
                r_start = sc_start_idx[sc]
                r_end   = sc_end_idx[sc]
                area    = sc_area_fractions[sc]
                H_snow_cur = H_snow[sc]
                fgi_cur    = fgi[sc]
                H_dc_sc    = H_dc[sc]

                # et_water_stress: scale ET demand by shallow-reservoir wetness
                # (PDM saturation CDF of the sub-catchment's top reservoir).
                if use_et_water_stress and pdm_H0_arr[r_start] > 0.0:
                    _h_shallow = H_res[r_start] if H_res[r_start] > 0.0 else 0.0
                    _et_stress = 1.0 - math.exp(-_h_shallow / pdm_H0_arr[r_start])
                else:
                    _et_stress = 1.0

                # ── Snowpack ──────────────────────────────────────────────
                sp_infiltrated = 0.0
                sp_deficit     = 0.0
                excess_dd      = 0.0

                if has_snowpack:
                    if use_et_reservoir_draw:
                        sp_recharge = P_t + H_dc_sc
                    else:
                        sp_recharge = P_t - ET_t * _et_stress + H_dc_sc

                    # Snowpack.recharge(sp_recharge)
                    if sp_recharge >= 0.0:
                        if T_t <= 0.0:
                            H_snow_cur += sp_recharge
                            # sp_infiltrated stays 0
                        else:
                            sp_infiltrated = sp_recharge
                    else:
                        # Sublimation; deficit if snow insufficient
                        if H_snow_cur > -sp_recharge:
                            H_snow_cur += sp_recharge
                        else:
                            sp_deficit = sp_recharge + H_snow_cur
                            H_snow_cur = 0.0
                        # sp_infiltrated stays 0

                    # Snowpack.melt(dt, P)
                    if T_t > 0.0:
                        pdd     = melt_factor * T_t * dt
                        ros_P   = P_t if use_rain_on_snow else 0.0
                        ros     = _CP_LF * T_t * ros_P
                        total_m = pdd + ros
                        if total_m <= H_snow_cur:
                            actual_melt = total_m
                        else:
                            actual_melt = H_snow_cur
                            if melt_factor > 0.0:
                                excess_dd = (total_m - actual_melt) / melt_factor
                        sp_infiltrated += actual_melt
                        H_snow_cur     -= actual_melt

                    # Temporarily store snowpack deficit as carry for res recharge
                    H_deficit_carry_for_res0 = sp_deficit
                else:
                    H_deficit_carry_for_res0 = H_dc_sc

                # ── FGI update ────────────────────────────────────────────
                # save calibrated value before possible override
                f0_frozen = f_dis[r_start]

                if use_fgi and not math.isinf(fdd_threshold):
                    T_eff = T_t * math.exp(-snow_insulation_k * H_snow_cur)

                    if has_trange:
                        T_max_t = T_max_arr[step]
                        T_min_t = T_min_arr[step]
                        DTR = T_max_t - T_min_t
                        # NaN comparisons → False → A_t = 1.0 when data absent
                        if DTR > 0.0 and T_max_t > 0.0:
                            f_above = T_max_t / DTR
                            if f_above > 1.0:
                                f_above = 1.0
                            A_t = 1.0 - (1.0 - fgi_decay_coeff) * f_above
                        else:
                            A_t = 1.0
                    else:
                        A_t = fgi_decay_coeff

                    new_fgi = A_t * fgi_cur - T_eff - excess_dd
                    if new_fgi < 0.0 or math.isnan(new_fgi):
                        new_fgi = 0.0
                    fgi_cur = new_fgi
                    if fgi_cur > fdd_threshold:
                        f_dis[r_start] = 1.0   # frozen ground: all → direct runoff

                # ── Reservoir cascade ─────────────────────────────────────
                qi = 0.0

                for r in range(r_start, r_end):
                    # --- Recharge ---
                    H_res_excess  = 0.0
                    H_res_deficit = 0.0

                    if r == r_start:
                        if has_snowpack:
                            _rech = sp_infiltrated + H_deficit_carry_for_res0
                        else:
                            if use_et_reservoir_draw:
                                _rech = P_t + H_deficit_carry_for_res0
                            else:
                                _rech = (P_t - ET_t * _et_stress
                                         + H_deficit_carry_for_res0)

                        if use_direct_runoff and _rech > 0.0:
                            q_direct = _rech * direct_runoff_fraction
                            qi      += q_direct
                            _rech   -= q_direct
                    else:
                        _rech = H_to_next[r - 1] + H_deficit_res[r - 1]

                    # PDM: shed the saturated fraction of positive recharge as
                    # saturation-excess overland flow before storage (exp CDF of
                    # storage capacities; mutually exclusive with a finite Hmax).
                    if pdm_H0_arr[r] > 0.0 and _rech > 0.0:
                        _f_sat = 1.0 - math.exp(-H_res[r] / pdm_H0_arr[r])
                        H_res_excess = _f_sat * _rech
                        _rech        = (1.0 - _f_sat) * _rech

                    new_H = H_res[r] + _rech
                    if new_H < 0.0:
                        H_res_deficit = new_H   # negative
                        H_res[r] = 0.0
                    elif (not math.isinf(Hmax_arr[r])) and new_H > Hmax_arr[r]:
                        H_res_excess = new_H - Hmax_arr[r]
                        H_res[r]     = Hmax_arr[r]
                    else:
                        H_res[r] = new_H

                    H_deficit_res[r] = H_res_deficit

                    # --- Discharge ---
                    b_r  = b_arr[r]
                    H0_r = H_res[r]

                    # Threshold junction: recession only above H_threshold
                    if junction_arr[r] == 2:
                        H_eff = H0_r - H_threshold_arr[r]
                        if H_eff < 0.0:
                            H_eff = 0.0
                    else:
                        H_eff = H0_r

                    if b_r == 1.0 or H_eff <= 0.0:
                        dH = H_eff * (1.0 - math.exp(-dt / tau_arr[r]))
                    else:
                        # Exact integration: H(t+dt) = [H^(1-b) + (b-1)dt/tau_eff]^(1/(1-b))
                        tau_eff = tau_arr[r] * (H_ref_arr[r] ** (b_r - 1.0))
                        arg     = H_eff ** (1.0 - b_r) + (b_r - 1.0) * dt / tau_eff
                        H_new_r = arg ** (1.0 / (1.0 - b_r)) if arg > 0.0 else 0.0
                        dH      = H_eff - H_new_r
                        if dH < 0.0:
                            dH = 0.0

                    # Partition dH: leakance applies to non-bottom reservoirs
                    is_bottom = (r == r_end - 1)
                    if junction_arr[r] == 1 and not is_bottom:
                        # Q_leak = min(dH, max(0, H_this - H_next) / R)
                        diff    = H0_r - H_res[r + 1]   # H_res[r+1] pre-recharge/discharge
                        if diff < 0.0:
                            diff = 0.0
                        Q_leak  = diff / leakance_R_arr[r]
                        if Q_leak > dH:
                            Q_leak = dH
                        H_to_next[r] = Q_leak
                        H_exfiltrated    = dH - Q_leak
                    else:
                        H_exfiltrated    = dH * f_dis[r]
                        H_to_next[r] = dH * (1.0 - f_dis[r])

                    H_discharge = H_res_excess + H_exfiltrated
                    H_res[r]   -= dH
                    qi         += H_discharge

                    # Tile drain: intercept f_tile fraction of H_to_next,
                    # route through a linear (b=1) sub-reservoir, to stream.
                    if f_tile_arr[r] > 0.0:
                        tile_in       = f_tile_arr[r] * H_to_next[r]
                        H_to_next[r] -= tile_in
                        H_tile[r]    += tile_in
                        tile_dH       = H_tile[r] * (1.0 - math.exp(-dt / tau_tile_arr[r]))
                        H_tile[r]    -= tile_dH
                        qi           += tile_dH

                    # Multipath drainage: a parallel fast path from this
                    # reservoir direct to stream, active when storage exceeds
                    # the threshold. Applied after primary recession via
                    # operator splitting (first-order accurate; acceptable at
                    # daily dt). Bypasses the junction/f_to_discharge partition.
                    if multipath_tau_arr[r] > 0.0:
                        H_above_mp = H_res[r] - multipath_thr_arr[r]
                        if H_above_mp > 0.0:
                            dH_mp        = H_above_mp * (
                                1.0 - math.exp(-dt / multipath_tau_arr[r]))
                            H_res[r]    -= dH_mp
                            qi          += dH_mp

                # Restore calibrated f_to_discharge (undo frozen-ground override)
                f_dis[r_start] = f0_frozen

                # Carry deficit from this sub-catchment's bottom reservoir
                H_dc_sc = H_deficit_res[r_end - 1]

                # ── ET draw from reservoirs ───────────────────────────────
                if use_et_reservoir_draw:
                    demand0 = et_alpha * ET_t
                    avail0  = H_res[r_start] if H_res[r_start] > 0.0 else 0.0
                    actual0 = demand0 if demand0 <= avail0 else avail0
                    H_res[r_start] -= actual0
                    # Condensation (negative ET) overtopping the cap runs off.
                    if H_res[r_start] > Hmax_arr[r_start]:
                        qi             += H_res[r_start] - Hmax_arr[r_start]
                        H_res[r_start]  = Hmax_arr[r_start]

                    if (r_end - r_start) >= 2:
                        r1 = r_start + 1
                        demand1 = (1.0 - et_alpha) * ET_t
                        if wp_soil > 0.0:
                            if wp_soil_sigma > 0.0:
                                # Gaussian CDF: demand scaled; draw from full H
                                f_avail = 0.5 * (1.0 + math.erf(
                                    (H_res[r1] - wp_soil)
                                    / (wp_soil_sigma * math.sqrt(2.0))))
                                demand1 = demand1 * f_avail
                                avail1  = H_res[r1] if H_res[r1] > 0.0 else 0.0
                            else:
                                # Hard threshold: no draw below wp_soil
                                if H_res[r1] <= wp_soil:
                                    demand1 = 0.0
                                    avail1  = 0.0
                                else:
                                    avail1 = H_res[r1] - wp_soil
                        else:
                            avail1 = H_res[r1] if H_res[r1] > 0.0 else 0.0
                        actual1 = demand1 if demand1 <= avail1 else avail1
                        H_res[r1] -= actual1
                        if H_res[r1] > Hmax_arr[r1]:
                            qi        += H_res[r1] - Hmax_arr[r1]
                            H_res[r1]  = Hmax_arr[r1]

                # Write back this sub-catchment's state and aggregate to basin
                H_snow[sc] = H_snow_cur
                fgi[sc]    = fgi_cur
                H_dc[sc]   = H_dc_sc

                qi_basin    += area * qi
                swe_tot     += area * H_snow_cur
                for r in range(r_start, r_end):
                    H_sub_total += area * (H_res[r] + H_tile[r])

            # ── Record basin outputs ──────────────────────────────────────
            Q_out[step]   = qi_basin
            SWE_out[step] = swe_tot
            for r in range(n_res):
                H_res_out[step, r] = H_res[r]
            H_sub_out[step] = H_sub_total

        return (Q_out, SWE_out, H_sub_out, H_res_out, H_res, H_tile,
                H_snow, fgi, H_dc)


class Reservoir(object):
    """
    Generic reservoir. Accepts new water (recharge), and sends it to other
    reservoirs and/or out of the system (discharge) at a rate that is
    proportional to the amount of water held in the reservoir.
    """

    def __init__(self, recession_coeff, f_to_discharge=1., Hmax=np.inf, pdm_H0=None,
                 H0=0., f_tile=0.0, tau_tile=None,
                 junction_type='fraction', leakance_R=None, H_threshold=0.0,
                 multipath_threshold=None, multipath_timescale=None):
        """
        Initialize a reservoir.

        Parameters
        ----------
        recession_coeff : float
            Recession coefficient. For a linear reservoir (recession
            exponent b = 1) this equals the true e-folding timescale, in
            days. For b > 1 it is not a timescale — its units are
            day·mm^(b-1) and the actual mean residence time depends on
            storage level; use :meth:`mean_residence_time` to obtain a
            physically comparable timescale at a given reference discharge.
        f_to_discharge : float, optional
            Fraction of water lost each time step that exits as river
            discharge. The remainder (1 - f_to_discharge) infiltrates to
            the next-deeper reservoir. Default 1.0 (all to discharge).
            Used only when junction_type is 'fraction' or 'threshold'.
        Hmax : float, optional
            Maximum water depth the reservoir can hold. Default np.inf.
        pdm_H0 : float or None, optional
            Characteristic storage depth [mm] for the probability-distributed
            model (PDM) of saturation-excess overland flow.  Storage capacity
            is assumed exponentially distributed across the catchment with
            mean pdm_H0; the saturated fraction when reservoir depth is H is
            f_sat = 1 - exp(-H / pdm_H0).  That fraction of each positive
            recharge pulse is shed immediately as overland flow.  Mutually
            exclusive with a finite Hmax.  Default None (PDM off).
        H0 : float, optional
            Initial water depth at the start of the simulation. Default 0.
        f_tile : float, optional
            Fraction of subsurface drainage (H_to_next) that is diverted
            to a fast tile-drain sub-reservoir instead of passing to the
            next-deeper reservoir.  The tile sub-reservoir drains directly to
            stream with e-folding time tau_tile.  This is a *fractional-bypass*
            tile representation: a constant fraction of recession outflow is
            routed through a downstream fast linear reservoir, regardless of
            current storage. Default 0 (no tile drainage).  Requires tau_tile
            when > 0.

            For a *threshold-activated parallel drainage* representation,
            in which a second outflow path from the same storage activates
            only above a water-table depth, use ``multipath_threshold`` /
            ``multipath_timescale`` instead.
        tau_tile : float or None, optional
            E-folding residence time [days] of the tile-drain sub-reservoir
            (used with ``f_tile``).  Typical values: 3–21 days for
            agricultural tile systems.  Required when f_tile > 0; ignored
            when f_tile == 0.  Default None.
        multipath_threshold : float or None, optional
            Storage depth [mm] above which a *parallel* fast drainage path
            from this reservoir activates.  The parallel path drains directly
            to stream with e-folding timescale ``multipath_timescale`` and is
            added to the discharge alongside the primary recession.  The
            outflow is

                Q_multipath = max(0, H - multipath_threshold) / multipath_timescale.

            Below the threshold this term is exactly zero, so the reservoir
            behaves identically to the no-multipath case.  Above it, a second
            linear drain operates in parallel with the primary recession,
            giving a two-timescale response: slow matrix drainage at low
            storage and faster combined drainage at high storage.

            Physical analogue: a tile-drain system installed at a fixed depth
            that activates only once the water table rises above the drain
            elevation.  Contrast with ``f_tile``/``tau_tile``, which is a
            constant-fraction fast-routing path through a separate downstream
            reservoir regardless of storage.

            Both ``multipath_threshold`` and ``multipath_timescale`` must be
            provided together, or both left as None (default = multipath
            disabled).
        multipath_timescale : float or None, optional
            E-folding timescale [days] of the parallel multipath drain.
            Typical values: 3–21 days for agricultural tile systems
            represented this way. Required together with
            ``multipath_threshold``; ignored when it is None.  Default None.
        junction_type : str, optional
            How drainage is partitioned between river discharge and the
            next-deeper reservoir.  Options:

            ``'fraction'`` (default)
                Fixed split: ``f_to_discharge`` fraction exits as stream
                discharge; remainder infiltrates.  Current default behaviour.

            ``'leakance'``
                Head-difference driven flow to the next reservoir, physically
                representing Darcy flow through a confining unit (e.g. a shale
                layer).  ``Q_leak = max(H_this - H_next, 0) / leakance_R``,
                capped at the total drainage ``dH``.  The remainder goes to
                stream.  Requires ``leakance_R``.

            ``'threshold'``
                Dead-storage threshold: the recession law is applied only to
                ``max(H - H_threshold, 0)``; water below ``H_threshold``
                never drains.  The above-threshold drainage splits by
                ``f_to_discharge`` as in the ``'fraction'`` case.  Models a
                stream-aquifer connection that activates only when the water
                table exceeds the streambed elevation.
        leakance_R : float or None, optional
            Leakance resistance [days].  ``Q_leak = ΔH / leakance_R``.
            Required when junction_type is ``'leakance'``.  Default None.
        H_threshold : float, optional
            Dead-storage threshold depth [mm].  Recession is applied only
            to ``max(H - H_threshold, 0)``.  Used when junction_type is
            ``'threshold'``.  Default 0.0 (no threshold; full storage drains).

        Raises
        ------
        ValueError
            If recession_coeff <= 0, f_to_discharge < 0 or > 1, Hmax < 0,
            pdm_H0 <= 0, f_tile < 0 or > 1, f_tile > 0 with no tau_tile,
            junction_type is unrecognised, junction_type is 'leakance'
            with no leakance_R, multipath_threshold < 0,
            multipath_timescale <= 0, or only one of the two multipath
            parameters is provided.
        """
        self.Hwater = H0
        self.Hmax = Hmax
        self.pdm_H0 = pdm_H0
        self.recession_coeff = recession_coeff
        self.f_to_discharge = f_to_discharge
        self.junction_type = junction_type
        self.leakance_R    = leakance_R
        self.H_threshold   = H_threshold

        # Initialized here so all instance attributes exist before
        # recharge() and discharge() are first called
        self.H_excess = 0.
        self.H_deficit = 0.
        self.H_exfiltrated = 0.
        self.H_to_next = 0.
        self.H_discharge = 0.
        # Fast-path components of H_discharge, recorded per step for
        # diagnostics / BMI flux partition (already included in H_discharge).
        self.H_tile = 0.
        self.H_multipath = 0.

        # Check values and note whether they are reasonable
        if recession_coeff <= 0:
            raise ValueError("recession_coeff must be > 0.")
        if f_to_discharge < 0:
            raise ValueError("Negative f_to_discharge not possible.")
        elif f_to_discharge > 1:
            raise ValueError("f_to_discharge: Cannot discharge >100% of water.")
        elif f_to_discharge == 0:
            warnings.warn("All water infiltrates when f_to_discharge is 0:"+
                          " you may have created a\n"+
                          "redundant pass-through water-storage layer",
                          stacklevel=2)
        if Hmax < 0:
            raise ValueError("Hmax must be >= 0 (and >0 makes more sense)")
        if pdm_H0 is not None and pdm_H0 <= 0:
            raise ValueError("pdm_H0 must be > 0")
        if f_tile < 0 or f_tile > 1:
            raise ValueError("f_tile must be in [0, 1]")
        if f_tile > 0 and tau_tile is None:
            raise ValueError("tau_tile must be provided when f_tile > 0")
        _valid_junctions = ('fraction', 'leakance', 'threshold')
        if junction_type not in _valid_junctions:
            raise ValueError(f"junction_type must be one of {_valid_junctions}")
        if junction_type == 'leakance' and leakance_R is None:
            raise ValueError("leakance_R must be provided when junction_type='leakance'")

        self.f_tile = f_tile
        if f_tile > 0.0 and tau_tile is not None:
            self.tile_res = Reservoir(tau_tile, f_to_discharge=1.0)
        else:
            self.tile_res = None

        # Multipath: threshold-activated parallel drain (distinct from the
        # constant-fraction f_tile/tau_tile bypass; see __init__ docstring).
        # Both must be set together or both None.
        if (multipath_threshold is None) ^ (multipath_timescale is None):
            raise ValueError(
                "multipath_threshold and multipath_timescale must be set "
                "together (both numbers) or both left as None.")
        if multipath_threshold is not None:
            if multipath_threshold < 0.0:
                raise ValueError("multipath_threshold must be >= 0.")
            if multipath_timescale <= 0.0:
                raise ValueError("multipath_timescale must be > 0.")
        self.multipath_threshold = multipath_threshold
        self.multipath_timescale = multipath_timescale

        # Power-law recession: Q = (H/τ) * (H/recession_H_ref)^(recession_exponent-1).
        # recession_exponent=1 (default) recovers the linear reservoir exactly.
        # recession_H_ref is the storage at which τ has its usual linear meaning [mm].
        self.recession_exponent = 1.0
        self.recession_H_ref    = 1.0

    @property
    def has_multipath(self):
        """True when both multipath parameters are configured (drain active)."""
        return (self.multipath_threshold is not None
                and self.multipath_timescale is not None)

    def recharge(self, H):
        """
        Add or remove water from the reservoir.

        Recharge H can be positive (net water input, e.g. P > ET) or
        negative (net deficit, e.g. ET > P). Sets H_excess if the
        reservoir overflows Hmax, or H_deficit if more water is removed
        than the reservoir holds.

        Parameters
        ----------
        H : float
            Depth of water added (positive) or removed (negative).

        Raises
        ------
        ValueError
            If Hwater is already negative before recharge is applied.
        """
        # Extra water above a maximum cap
        self.H_excess = 0.
        # Water that this layer cannot hold and cannot be passed to a deeper layer
        self.H_deficit = 0.

        # ERROR if water is less than 0 -- may be able to remove
        # this check later
        if self.Hwater < 0:
            raise ValueError("Hwater in reservoir < 0; non-physical")

        # PDM: exponential distribution of storage capacities.
        # Saturated fraction of catchment sheds positive recharge immediately
        # as saturation-excess overland flow before the remainder enters storage.
        if self.pdm_H0 is not None and H > 0:
            f_sat = 1.0 - np.exp(-self.Hwater / self.pdm_H0)
            self.H_excess = f_sat * H
            H = (1.0 - f_sat) * H

        # What if more water is lost during "recharge" than exists in reservoir?
        # Create a deficit and bring Hwater to 0
        if self.Hwater + H < 0:
            self.H_deficit += self.Hwater + H
            self.Hwater = 0.
        # What if more water is added than maximum reservoir capacity?
        # Mark excess (straight to runoff) and bring Hwater to Hmax
        elif self.Hwater + H > self.Hmax:
            self.H_excess += self.Hwater + H - self.Hmax
            self.Hwater = self.Hmax
        # Otherwise, we're in a range in which 0 <= H <= Hmax
        # Yay! Things are easier!
        else:
            self.Hwater += H

    def discharge(self, dt, H_next=None):
        """
        Discharge water from the reservoir over one time step.

        Computes water lost by the recession law, partitions it between
        river discharge (H_exfiltrated) and infiltration to the next-deeper
        reservoir (H_to_next) according to junction_type, and adds
        overflow from recharge() (H_excess) to H_discharge.

        Parameters
        ----------
        dt : float
            Time step duration (same units as recession_coeff; typically days).
        H_next : float or None, optional
            Current water depth of the next-deeper reservoir [mm].  Used
            only when junction_type is ``'leakance'`` to compute the
            head-difference driven flux.  Ignored for other junction types.
        """
        b   = self.recession_exponent
        H0  = self.Hwater

        # For threshold junction: recession applies only above H_threshold.
        # Water below H_threshold is permanent dead storage.
        H_eff = max(0.0, H0 - self.H_threshold) if self.junction_type == 'threshold' else H0

        if b == 1.0 or H_eff <= 0.0:
            dH = H_eff * (1 - np.exp(-dt / self.recession_coeff))
        else:
            # Exact integration of dH/dt = -(H/recession_coeff)·(H/H_ref)^(b-1)
            #   = -H^b / (recession_coeff · H_ref^(b-1))
            # Substituting u = H^(1-b):  du/dt = (b-1)/tau_eff
            # => H(t+dt) = [H0^(1-b) + (b-1)·dt/tau_eff]^(1/(1-b))
            tau_eff = self.recession_coeff * self.recession_H_ref ** (b - 1.0)
            H_new   = (H_eff ** (1.0 - b) + (b - 1.0) * dt / tau_eff) ** (1.0 / (1.0 - b))
            dH      = H_eff - max(0.0, H_new)

        # Partition dH between stream discharge and infiltration to next reservoir.
        if self.junction_type == 'leakance' and H_next is not None:
            # Head-difference driven flux through confining unit (e.g., shale layer).
            # Q_leak = max(H_this - H_next, 0) / R, capped at available drainage dH.
            Q_leak          = min(dH, max(0.0, H0 - H_next) / self.leakance_R)
            self.H_to_next  = Q_leak
            self.H_exfiltrated = dH - Q_leak
        else:
            # 'fraction' and 'threshold': fixed f_to_discharge split.
            self.H_exfiltrated = dH * self.f_to_discharge
            self.H_to_next     = dH * (1 - self.f_to_discharge)

        self.H_discharge = self.H_excess + self.H_exfiltrated
        self.Hwater -= dH

        # Reset the per-step fast-path components (recorded for the BMI flux
        # partition; they are also folded into H_discharge below).
        self.H_tile = 0.0
        self.H_multipath = 0.0

        # Tile drainage: divert f_tile of H_to_next to the fast sub-reservoir.
        # The remainder continues to the next-deeper reservoir as normal.
        if self.tile_res is not None:
            tile_in = self.f_tile * self.H_to_next
            self.H_to_next -= tile_in
            self.tile_res.recharge(tile_in)
            self.tile_res.discharge(dt)
            self.H_tile = self.tile_res.H_discharge
            self.H_discharge += self.H_tile

        # Multipath drainage: a parallel fast path direct to stream, active
        # only when storage exceeds multipath_threshold. Applied after the
        # primary recession via operator splitting (first-order accurate in dt;
        # acceptable for daily timesteps with τ_matrix >> dt). Bypasses the
        # junction/f_to_discharge partition — goes straight to discharge.
        if self.has_multipath:
            H_above = self.Hwater - self.multipath_threshold
            if H_above > 0.0:
                dH_mp = H_above * (1.0 - np.exp(-dt / self.multipath_timescale))
                self.Hwater     -= dH_mp
                self.H_multipath = dH_mp
                self.H_discharge += dH_mp

    def mean_residence_time(self, Q_ref):
        """
        Mean residence time [days] at a reference steady-state discharge.

        For a nonlinear reservoir governed by
        :math:`Q = (H/\\tau)\\cdot(H/H_{\\mathrm{ref}})^{b-1}
        = H^{b}/\\tau_{\\mathrm{eff}}`, with effective constant
        :math:`\\tau_{\\mathrm{eff}} = \\tau\\,H_{\\mathrm{ref}}^{\\,b-1}`,
        the steady-state storage at discharge *Q_ref* is
        :math:`H_{ss} = (Q_{\\mathrm{ref}} \\cdot \\tau_{\\mathrm{eff}})^{1/b}`,
        giving

        .. math::

            \\mathrm{MRT} = \\frac{H_{ss}}{Q_{\\mathrm{ref}}} =
                \\frac{\\tau_{\\mathrm{eff}}^{1/b}}{Q_{\\mathrm{ref}}^{\\,1 - 1/b}}

        With the default :math:`H_{\\mathrm{ref}} = 1` this is
        :math:`\\tau^{1/b} / Q_{\\mathrm{ref}}^{\\,1-1/b}`.

        For a linear reservoir (*b* = 1) this reduces exactly to
        :math:`\\tau`.  For *b* > 1, MRT is smaller than :math:`\\tau`
        whenever :math:`Q_{\\mathrm{ref}} > 1` mm/day, reflecting the
        faster drainage at realistic operating storage depths.

        Unlike :attr:`recession_coeff`, which equals a timescale only when b=1
        (linear reservoir), MRT is a physically comparable timescale across
        reservoirs with different recession exponents.

        Parameters
        ----------
        Q_ref : float
            Representative steady-state discharge from this reservoir
            [mm/day]. Use the long-term mean flux attributed to this
            layer (e.g. mean annual discharge partitioned by reservoir).

        Returns
        -------
        float
            Mean residence time [days].

        Raises
        ------
        ValueError
            If Q_ref <= 0.
        """
        if Q_ref <= 0:
            raise ValueError("Q_ref must be > 0.")
        b = self.recession_exponent
        if b == 1.0:
            return self.recession_coeff
        tau_eff = self.recession_coeff * self.recession_H_ref ** (b - 1.0)
        return tau_eff ** (1.0 / b) / Q_ref ** (1.0 - 1.0 / b)


class Snowpack(object):
    """
    Snowpack reservoir driven by temperature.

    Accumulates precipitation as snow when mean temperature is at or below
    0 °C. Melts at a positive-degree-day rate when temperature is above 0 °C.
    All melt is routed to the top subsurface reservoir as infiltration; direct
    discharge to the river is not modeled.

    Should precede the subsurface reservoir list in a watershed model.
    The melt factor is a positive-degree-day factor [mm/°C/day].
    """

    def __init__(self, melt_factor=None):
        """
        Initialize an empty snowpack.

        Parameters
        ----------
        melt_factor : float, optional
            Positive-degree-day melt factor (mm SWE °C⁻¹ day⁻¹).
            Can be set or updated later via set_melt_factor().
        """
        self.Hwater = 0.  # SWE
        self.melt_factor = melt_factor
        self.T = 0.
        self.H_infiltrated = 0.
        self.H_deficit = 0.

    def set_melt_factor(self, melt_factor):
        """
        Set or update the positive-degree-day melt factor.

        Parameters
        ----------
        melt_factor : float
            Melt rate per positive degree-day (mm SWE °C⁻¹ day⁻¹).
        """
        self.melt_factor = melt_factor

    def set_temperature(self, T):
        """
        Set the mean air temperature for the current time step.

        Parameters
        ----------
        T : float
            Mean air temperature (°C).
        """
        self.T = T

    def recharge(self, H):
        """
        Apply net water input or deficit to the snowpack.

        If T <= 0, positive H accumulates as snow (SWE). If T > 0,
        positive H bypasses the snowpack and is passed directly to the
        top subsurface reservoir via H_infiltrated. Negative H (ET > P)
        is removed from the snowpack as sublimation; any remainder that
        exceeds available SWE becomes H_deficit.

        Parameters
        ----------
        H : float
            Net water depth for this time step (mm). Positive = input
            (P - ET > 0); negative = deficit (ET - P > 0).
        """

        self.H_deficit = 0.  # Water deficit with neg ET; just this time step
        # If positive recharge
        if H >= 0:
            if self.T <= 0:
                self.Hwater += H
                self.H_infiltrated = 0.
            else:
                # Incoming precip component; melt sums with this
                # This is then directly passed to the first layer of the
                # set of hydrological reservoirs
                self.H_infiltrated = H
        # If negative recharge: remove water from snowpack via sublimation.
        # Any deficit beyond available SWE is passed down as H_deficit.
        else:
            # Sublimation (effectively) if snow present;
            # Otherwise pass water deficit
            if self.Hwater > -H:
                self.Hwater += H
            else:
                self.H_deficit += H + self.Hwater
                self.Hwater = 0
            self.H_infiltrated = 0.

    def melt(self, dt, P=0.0):
        """
        Compute positive-degree-day and rain-on-snow melt; update state.

        Both terms are routed to H_infiltrated (→ top soil reservoir). If
        total available energy exceeds the SWE present, the leftover is
        returned as equivalent degree-days so the caller can credit it
        toward frozen-soil thawing (FGI reduction) rather than losing it.

        Rain-on-snow sensible heat: water arriving at T_mean > 0 °C
        carries (c_p / L_f) · T · P mm SWE of latent-heat capacity.
        Spring snowpacks are near-isothermal at 0 °C, so cold-content
        corrections are negligible and the latent-heat term dominates.

        References
        ----------
        McCabe et al. (2007) doi:10.1175/BAMS-88-3-319
        Würzer et al. (2016) doi:10.1175/JHM-D-15-0181.1

        Parameters
        ----------
        dt : float
            Timestep [days].
        P : float, optional
            Raw liquid precipitation [mm/day]. Used to compute rain-on-snow
            sensible-heat melt. Default 0 (PDD only).

        Returns
        -------
        excess_dd : float
            Leftover melt energy after the snowpack is exhausted, expressed
            as degree-day equivalent [°C·day] = leftover mm SWE / melt_factor.
            Zero when SWE is not fully depleted.

            The melt factor (mm SWE °C⁻¹ day⁻¹) serves as the bridge
            between the PDD snowmelt representation and the frozen ground
            index (FGI): dividing excess melt depth (mm SWE) by melt_factor
            recovers the equivalent thermal forcing in °C·day, which is the
            currency the FGI uses. See Buckets._update_fgi().
        """
        if self.T <= 0:
            return 0.0

        pdd_avail   = self.melt_factor * self.T * dt    # [mm SWE]
        ros_avail   = _CP_LF * self.T * P               # [mm SWE] rain-on-snow
        total_avail = pdd_avail + ros_avail

        if total_avail <= self.Hwater:
            actual_melt = total_avail
            excess_dd   = 0.0
        else:
            actual_melt = self.Hwater
            # Leftover energy → °C·day equivalent for soil-thaw credit
            excess_dd = (total_avail - actual_melt) / self.melt_factor

        self.H_infiltrated += actual_melt
        self.Hwater        -= actual_melt
        return excess_dd


class SubCatchment(object):
    """
    A parallel hydraulic compartment of the basin with its own reservoir
    cascade and per-compartment state.

    Sub-catchments represent spatially distinct zones (e.g. till uplands vs.
    lake-clay lowlands) that drain to the same river channel *in parallel*
    rather than in a vertical cascade. Each sub-catchment is internally an
    ordinary serial reservoir cascade; basin discharge is the area-weighted
    mean of the sub-catchments' per-unit-area discharges. A model with a
    single sub-catchment of ``area_fraction = 1.0`` reproduces the original
    (non-partitioned) MNiShed behaviour exactly.

    State that evolves independently per zone — snowpack, frozen-ground index,
    and the carried recharge deficit — lives on the sub-catchment so that
    per-sub-catchment forcing can be added later without changing the public
    API. In this version the forcing is basin-shared, so these states evolve
    identically across sub-catchments.

    Parameters
    ----------
    name : str
        Label for the sub-catchment. Must be unique within a basin.
    area_fraction : float
        Fraction of the basin area drained by this sub-catchment, in [0, 1].
        The fractions across all sub-catchments must sum to 1.
    reservoirs : list of Reservoir
        The vertical reservoir cascade for this sub-catchment, ordered
        shallowest (index 0) to deepest. Must contain at least one reservoir.
    snowpack : Snowpack or None, optional
        Snowpack state for this sub-catchment. None when the basin has no
        snowpack module enabled. Default None.
    fgi_init : float, optional
        Initial frozen-ground index [°C·day] for this sub-catchment.
        Default 0.0.
    H_deficit_carry_init : float, optional
        Initial recharge deficit [mm] carried into the first time step.
        Default 0.0.
    forcing : dict or None, optional
        Per-sub-catchment forcing specification. Reserved for a future
        release; supplying it raises NotImplementedError, because the current
        implementation shares basin-level forcing across all sub-catchments.
        Default None.
    """

    def __init__(self, name, area_fraction, reservoirs, *,
                 snowpack=None, fgi_init=0.0, H_deficit_carry_init=0.0,
                 forcing=None):
        if forcing is not None:
            raise NotImplementedError(
                "Per-sub-catchment forcing is not yet supported; basin-level "
                "forcing is shared across all sub-catchments. Remove the "
                f"'forcing' block from sub-catchment '{name}'.")
        if len(reservoirs) < 1:
            raise ValueError(
                f"Sub-catchment '{name}' has no reservoirs; each "
                "sub-catchment must have at least one reservoir.")
        self.name             = name
        self.area_fraction    = area_fraction
        self.reservoirs       = reservoirs
        self.snowpack         = snowpack
        self.fgi              = fgi_init
        self.H_deficit_carry  = H_deficit_carry_init
        self.forcing          = forcing


class Buckets(object):
    """
    Incorporates a list of reservoirs into a linear hierarchy that sends water
    either downwards or out to the surface. Reservoirs are ordered from top
    (nearest Earth's surface) to bottom (deepest groundwater); this order
    controls the direction of infiltration between layers.

    MNiShed is designed as a daily-timestep model. This is a deliberate
    design choice: the physical parameterisations — degree-day snowmelt,
    Thornthwaite ET, and linear reservoir drainage — are climatological
    approximations that are well-founded at daily resolution but lose physical
    meaning at finer scales. The daily timestep is enforced in initialize().
    """

    def __init__(self, T_monthly_normals=None):
        """
        Initialize the watershed model.

        If using the ThornthwaiteChang2019 ET method, pass
        T_monthly_normals here so that the thermal index I and exponent
        a are computed once from climatological normals and remain fixed
        throughout the simulation.

        Parameters
        ----------
        T_monthly_normals : array-like of length 12, optional
            Long-term mean monthly temperatures (°C) used to compute the
            Thornthwaite thermal index I and exponent a per Chang et al.
            (2019), https://doi.org/10.1002/ird.2309. Required when
            evapotranspiration_method is 'ThornthwaiteChang2019'.
        """
        # Thornthwaite thermal index and exponent, per Chang et al. (2019)
        # https://doi.org/10.1002/ird.2309
        # I is climatologically imposed by the local normal temperature regime
        # and must remain fixed during simulation (not recomputed each timestep).
        if T_monthly_normals is not None:
            self.Chang_I = self._compute_Chang_I(T_monthly_normals)
            self.Chang_a = self._compute_Chang_a(self.Chang_I)

        # Frozen ground index (Molnau & Bissell 1983).  Disabled by default
        # (threshold = inf); overridden by snowmelt.fdd_threshold in YAML or
        # by run_and_score(fdd_threshold=).  The FGI state itself lives on each
        # sub-catchment (SubCatchment.fgi), created in initialize().
        self.fdd_threshold = np.inf  # [°C·day]

    @property
    def reservoirs(self):
        """
        Flat list of every reservoir across all sub-catchments, ordered by
        sub-catchment and then shallowest-to-deepest within each.

        Backward-compatible view: for the common single-sub-catchment case
        this is exactly the reservoir cascade as before. Reservoirs may be
        mutated through this list (it returns the live objects), but the list
        itself is rebuilt on each access, so assign to ``sub_catchments`` to
        change the structure.
        """
        return [r for sc in self.sub_catchments for r in sc.reservoirs]

    @reservoirs.setter
    def reservoirs(self, value):
        # Backward-compatible assignment for programmatic construction:
        # ``buckets.reservoirs = [...]`` replaces the structure with a single
        # sub-catchment holding the given cascade (area_fraction = 1.0). Any
        # existing snowpack reference is preserved. To build multiple parallel
        # sub-catchments, assign to ``sub_catchments`` directly.
        snowpack = (self.sub_catchments[0].snowpack
                    if getattr(self, 'sub_catchments', None) else None)
        self.sub_catchments = [SubCatchment('basin', 1.0, list(value),
                                            snowpack=snowpack)]

    @property
    def n_sub_catchments(self):
        """Number of parallel sub-catchments (>= 1)."""
        return len(self.sub_catchments)

    @property
    def snowpack(self):
        """
        Snowpack of the first sub-catchment (the whole basin when there is a
        single sub-catchment). For per-sub-catchment snowpacks, iterate
        ``sub_catchments``. Read-only; snowpacks are created in initialize().
        """
        return self.sub_catchments[0].snowpack

    @property
    def _fgi(self):
        """
        Frozen-ground index of the first sub-catchment [°C·day]. For
        per-sub-catchment FGI, see ``sub_catchments``. The state is advanced
        per sub-catchment by _update_fgi(); the setter broadcasts a scalar to
        every sub-catchment (used to restore or reset single-sub-catchment
        calibration state).
        """
        return self.sub_catchments[0].fgi

    @_fgi.setter
    def _fgi(self, value):
        for sc in self.sub_catchments:
            sc.fgi = value

    @property
    def H_deficit_carry(self):
        """
        Basin recharge deficit carried into the next time step [mm], as the
        area-weighted sum over sub-catchments (so the basin mass balance is
        exact at any number of sub-catchments). The per-sub-catchment deficits
        are stored on each SubCatchment.

        The setter broadcasts a scalar to every sub-catchment; because the area
        fractions sum to 1, the area-weighted getter then returns that same
        value, so a set/get round-trips. Used to restore or reset state in
        single-sub-catchment calibration.
        """
        return sum(sc.area_fraction * sc.H_deficit_carry
                   for sc in self.sub_catchments)

    @H_deficit_carry.setter
    def H_deficit_carry(self, value):
        for sc in self.sub_catchments:
            sc.H_deficit_carry = value

    def _compute_Chang_I(self, T_monthly_normals):
        """
        Compute the Thornthwaite thermal index I from long-term monthly normal
        temperatures, per Chang et al. (2019), Eq. 1.
        https://doi.org/10.1002/ird.2309

        Parameters
        ----------
        T_monthly_normals : array-like, length 12
            Long-term mean monthly temperatures (°C). Negative values are
            treated as 0 per the Thornthwaite convention.

        Returns
        -------
        I : float
            Thermal index (dimensionless).
        """
        Tn = np.maximum(T_monthly_normals, 0)
        return np.sum((0.2 * Tn) ** 1.514)

    def _compute_Chang_a(self, I):
        """
        Compute the Thornthwaite exponent a from thermal index I, per
        Chang et al. (2019), Eq. 1.
        https://doi.org/10.1002/ird.2309

        Parameters
        ----------
        I : float
            Thermal index, as returned by _compute_Chang_I.

        Returns
        -------
        a : float
            Thornthwaite exponent (dimensionless).
        """
        return (6.75e-7 * I**3
                - 7.71e-5 * I**2
                + 1.7912e-2 * I
                + 0.49239)

    def export_Hlist(self):
        """
        Return the current water depths in all subsurface reservoirs.

        Useful for checkpointing reservoir state between a spin-up run
        and the main simulation, or for restarting a paused run.

        Returns
        -------
        list of float
            Water depth in each reservoir, ordered from shallowest
            (index 0) to deepest.
        """
        return [reservoir.Hwater for reservoir in self.reservoirs]

    def _depth_column_names(self):
        """
        Column names for per-reservoir stored depths (``store_depths=True``),
        aligned with the flat :attr:`reservoirs` order.

        A single sub-catchment keeps the legacy ``'H_reservoir_{i} (modeled)
        [mm]'`` names. With several sub-catchments the names are prefixed with
        the sub-catchment name and use the reservoir's index *within* that
        sub-catchment, e.g. ``'H_till_uplands_reservoir_0 (modeled) [mm]'``,
        so each depth column is unambiguously attributable to a zone.
        """
        if self.n_sub_catchments == 1:
            return [f'H_reservoir_{i} (modeled) [mm]'
                    for i in range(len(self.sub_catchments[0].reservoirs))]
        return [f'H_{sc.name}_reservoir_{j} (modeled) [mm]'
                for sc in self.sub_catchments
                for j in range(len(sc.reservoirs))]

    @staticmethod
    def _build_reservoir_cascade(res_cfg, H0_list):
        """
        Build a vertical reservoir cascade from a "reservoirs" config block.

        Parameters
        ----------
        res_cfg : dict
            A "reservoirs" config block: per-reservoir parameter lists keyed
            by name (recession_coefficients, exfiltration_fractions, ...). All
            lists must have the same length as H0_list.
        H0_list : list of float
            Initial water depth [mm] for each reservoir, ordered shallowest to
            deepest. Its length sets the number of reservoirs in the cascade.

        Returns
        -------
        list of Reservoir
            The cascade, ordered shallowest (index 0) to deepest.
        """
        n = len(H0_list)
        for _key, _vals in res_cfg.items():
            if len(_vals) != n:
                raise ValueError(
                    f'"{_key}" within "reservoirs" has {len(_vals)} entries, '
                    f'implying a different number of subsurface water '
                    f'reservoirs than the {n} initial water depths '
                    f'(water_reservoir_effective_depths__mm).')

        _pdm_H0        = res_cfg.get('pdm_H0__mm',               [None]  * n)
        _f_tile        = res_cfg.get('tile_fractions',           [0.0]   * n)
        _tau_tile      = res_cfg.get('tile_residence_times__days', [None] * n)
        _recession_exp = res_cfg.get('recession_exponents',     [1.0]   * n)
        _junction_types = res_cfg.get('junction_types',     ['fraction'] * n)
        _leakance_R    = res_cfg.get('leakance_R__days',        [None]  * n)
        _H_threshold   = res_cfg.get('H_threshold__mm',         [0.0]   * n)
        _mp_threshold  = res_cfg.get('multipath_thresholds__mm', [None] * n)
        _mp_timescale  = res_cfg.get('multipath_timescales__days', [None] * n)

        reservoirs = [
            Reservoir(
                recession_coeff = res_cfg['recession_coefficients'][i],
                f_to_discharge  = res_cfg['exfiltration_fractions'][i],
                Hmax            = res_cfg['maximum_effective_depths__mm'][i],
                pdm_H0          = _pdm_H0[i],
                H0              = H0_list[i],
                f_tile          = _f_tile[i],
                tau_tile        = _tau_tile[i],
                junction_type   = _junction_types[i],
                leakance_R      = _leakance_R[i],
                H_threshold     = _H_threshold[i] if _H_threshold[i] is not None else 0.0,
                multipath_threshold = _mp_threshold[i],
                multipath_timescale = _mp_timescale[i],
            )
            for i in range(n)]
        for i, b_exp in enumerate(_recession_exp):
            reservoirs[i].recession_exponent = float(b_exp)
        return reservoirs

    def _build_sub_catchments(self, sc_cfgs):
        """
        Build the list of parallel SubCatchments from a "sub_catchments" config.

        Parameters
        ----------
        sc_cfgs : list of dict
            One entry per sub-catchment. Each requires a unique ``name``, an
            ``area_fraction``, and a ``reservoirs`` block (same schema as the
            top-level legacy block). An optional ``initial_conditions`` block
            supplies per-reservoir initial depths
            (water_reservoir_effective_depths__mm); absent depths default to 0.
            A ``forcing`` block is reserved for a future release and raises
            NotImplementedError if present.

        Returns
        -------
        list of SubCatchment

        Raises
        ------
        ValueError
            If names are not unique, or the area fractions do not sum to 1.
        """
        if not isinstance(sc_cfgs, list) or len(sc_cfgs) < 1:
            raise ValueError(
                '"sub_catchments" must be a non-empty list of sub-catchment '
                'definitions.')

        _names = [sc['name'] for sc in sc_cfgs]
        if len(set(_names)) != len(_names):
            raise ValueError(
                f'Sub-catchment names must be unique; got {_names}.')

        _area_sum = sum(sc['area_fraction'] for sc in sc_cfgs)
        if abs(_area_sum - 1.0) > 1e-6:
            raise ValueError(
                f'Sub-catchment area_fraction values must sum to 1.0; '
                f'got {_area_sum:.6f}.')

        sub_catchments = []
        for sc in sc_cfgs:
            _res_cfg = sc['reservoirs']
            _ic = sc.get('initial_conditions', {})
            _n = len(_res_cfg['recession_coefficients'])
            _H0 = _ic.get('water_reservoir_effective_depths__mm', [0.0] * _n)
            _reservoirs = self._build_reservoir_cascade(_res_cfg, _H0)
            sub_catchments.append(
                SubCatchment(sc['name'], sc['area_fraction'], _reservoirs,
                             forcing=sc.get('forcing')))
        return sub_catchments

    def initialize(self, config_file=None, enforce_water_balance=None):
        """
        Set up the model from a YAML configuration file.

        Reads the configuration file, loads the input time series, builds
        the reservoir stack, instantiates snowpack if temperature data are
        present, optionally computes the water-year ET multiplier, and runs
        any requested spin-up cycles. Part of the CSDMS Basic Model Interface.

        Parameters
        ----------
        config_file : str, optional
            Path to the YAML configuration file. If None, all required
            values must be set on the object directly before calling
            update().
        enforce_water_balance : str or None, optional
            How to scale ET to close the water balance. When None (default),
            the value is read from ``general: enforce_water_balance`` in the
            YAML config, which itself defaults to ``'water-year'`` if absent.
            Accepted values:

            * ``'water-year'`` — Scale ET by a per-water-year multiplier so
              that P - Q - ET = 0 over each water year.
            * ``'global'`` — Scale ET by a single multiplier computed from
              sum(P - Q_obs) / sum(ET_raw) over the full record. No
              per-year overfitting; does not add hidden degrees of freedom
              to AIC comparisons.
            * ``'none'`` — Use raw ET without correction. Appropriate only
              when supplying trusted measured ET (e.g. eddy covariance).
              Using ``'none'`` with ThornthwaiteChang2019 will raise a
              warning because Thornthwaite ET carries large systematic
              biases.
        """
        if config_file is None:
            warnings.warn("No configuration file provided; all values needed "+
                          "for a model run therefore must be set independently.",
                          stacklevel=2)

        # Parse YAML configuration file
        # And assign variables except for optimization bounds and plotting
        if config_file is not None:
            try:
                with open(config_file, "r") as yamlfile:
                    self.cfg = yaml.load(yamlfile, Loader=yaml.FullLoader)
            except FileNotFoundError as e:
                # Raise rather than sys.exit so library / BMI / notebook callers
                # can handle a bad config; the CLI (main()) translates this into
                # a clean exit code.
                raise FileNotFoundError(
                    f"Config file not found: {config_file}") from e
            except yaml.YAMLError as e:
                raise yaml.YAMLError(
                    f"Could not parse config file {config_file}:\n{e}") from e

        # Read input time series from the CSV path specified in the config
        self.hydrodata = pd.read_csv(
            self.cfg['timeseries']['datafile'],
            parse_dates=['Date'])

        # Build the sub-catchment structure. A config may either give a single
        # basin-wide reservoir cascade (legacy: top-level "reservoirs" +
        # "initial_conditions" blocks) or partition the basin into parallel
        # "sub_catchments", each with its own cascade. The legacy form is
        # promoted to one sub-catchment of area_fraction 1.0, which reproduces
        # the original behaviour exactly (self.reservoirs then exposes the same
        # flat cascade as before).
        if 'sub_catchments' in self.cfg:
            self.sub_catchments = self._build_sub_catchments(
                self.cfg['sub_catchments'])
        else:
            _basin_reservoirs = self._build_reservoir_cascade(
                self.cfg['reservoirs'],
                self.cfg['initial_conditions'][
                    'water_reservoir_effective_depths__mm'])
            self.sub_catchments = [
                SubCatchment('basin', 1.0, _basin_reservoirs)]
        # Total reservoir count across all sub-catchments (informational;
        # equals the legacy single-cascade length when there is one sub-catchment).
        self.n_reservoirs = len(self.reservoirs)

        # Check if bottom reservoir discharges all to river: conserve mass.
        # But allow through with a warning in case the user wants a
        # deep and non-discharging reservoir (although this could be set up
        # explicitly too).
        if self.reservoirs[-1].f_to_discharge < 1:
            warnings.warn("f_to_discharge of bottom water-storage layer < 1.\n"+
                          "You are not conserving mass.", stacklevel=2)

        # Set scalar variables based on yaml
        self.melt_factor         = self.cfg['snowmelt']['PDD_melt_factor']
        self.snow_insulation_k   = self.cfg['snowmelt'].get('snow_insulation_k',   0.0)
        self.fgi_decay_coeff     = self.cfg['snowmelt'].get('fgi_decay_coeff',     0.97)
        self.fdd_threshold       = self.cfg['snowmelt'].get('fdd_threshold',       np.inf)
        self.et_method = self.cfg['catchment']['evapotranspiration_method']
        if self.et_method == 'ThornthwaiteChang2019' and not hasattr(self, 'Chang_I'):
            n_years = len(self.hydrodata) / 365.25
            if n_years < 20:
                warnings.warn(
                    f"ThornthwaiteChang2019: monthly temperature normals were not provided "
                    f"and will be computed from the {n_years:.1f}-year input record. "
                    f"For short records this may not represent long-term climatology. "
                    f"Pass T_monthly_normals to Buckets() for reliable results.",
                    UserWarning, stacklevel=2,
                )
            T_normals = (self.hydrodata['Mean Temperature [C]']
                             .groupby(pd.DatetimeIndex(self.hydrodata['Date']).month)
                             .mean()
                             .reindex(range(1, 13))
                             .values)
            self.Chang_I = self._compute_Chang_I(T_normals)
            self.Chang_a = self._compute_Chang_a(self.Chang_I)
        self.water_year_start_month = self.cfg['catchment']['water_year_start_month']
        self.drainage_basin_area__km2 = self.cfg['catchment']['drainage_basin_area__km2']
        self.baseflow_Q = self.cfg['catchment'].get('baseflow_Q', 0.0)
        # Per-step flux-partition components, populated by update() and read
        # by the BMI wrapper (mm/day).  baseflow_Q itself is an output-layer
        # term applied by run_and_score / the BMI, not by the cascade.
        self._flux_direct_runoff = 0.0
        self._flux_tile = 0.0
        self._flux_multipath = 0.0

        # Module enable/disable flags — read from config, default to on
        # (except direct_runoff, which defaults to off).
        _modules = self.cfg.get('modules', {})
        self.use_snowpack        = _modules.get('snowpack',        True)
        self.use_frozen_ground   = _modules.get('frozen_ground',   True)
        self.use_rain_on_snow    = _modules.get('rain_on_snow',    True)
        self.use_direct_runoff   = _modules.get('direct_runoff',   False)
        self.use_dtr_fgi_decay   = _modules.get('dtr_fgi_decay',   True)
        self.use_et_water_stress  = _modules.get('et_water_stress',  False)
        self.use_et_reservoir_draw = _modules.get('et_reservoir_draw', False)
        if self.use_et_water_stress and self.use_et_reservoir_draw:
            warnings.warn(
                "et_water_stress and et_reservoir_draw are mutually exclusive. "
                "et_reservoir_draw will be used; et_water_stress ignored.",
                UserWarning, stacklevel=2,
            )
            self.use_et_water_stress = False
        self.et_scale = 1.0   # universal ET multiplier; default 1.0 (no-op); override via
        #                       run_and_score(et_scale=) or as a free calibration parameter
        # Fraction of ET_pot drawn from reservoir 0 (shallow); 1-et_alpha from reservoir 1.
        # Read from general: et_alpha in config YAML; override via run_and_score(et_alpha=).
        self.et_alpha = self.cfg['general'].get('et_alpha', 1.0)
        # Wilting-point threshold for soil reservoir ET draw [mm].
        # wp_soil_sigma > 0 enables the spatially variable (Gaussian CDF) form.
        self.wp_soil = 0.0
        self.wp_soil_sigma = 0.0

        # Check if there is a mean temperature column for snowpack.
        # If not, note that no snowpack processes will be included
        self.has_snowpack = (self.use_snowpack and
                             'Mean Temperature [C]' in self.hydrodata.columns)

        # Detect optional daily T_min / T_max for DTR-based FGI decay.
        self._has_trange = (self.use_dtr_fgi_decay and
                            'Minimum Temperature [C]' in self.hydrodata.columns and
                            'Maximum Temperature [C]' in self.hydrodata.columns)
        if self.has_snowpack:
            # Each sub-catchment gets its own snowpack so that per-sub-catchment
            # forcing can drive them independently later. With basin-shared
            # forcing they evolve identically. self.snowpack aliases the first.
            for sc in self.sub_catchments:
                sc.snowpack = Snowpack(self.melt_factor)
        elif 'Mean Temperature [C]' in self.hydrodata.columns and not self.use_snowpack:
            pass  # snowpack deliberately disabled via modules config
        else:
            warnings.warn('"Mean Temperature [C]" has not been set. '
                          'No snowpack processes will be simulated.', stacklevel=2)

        # How many times to loop the full time series for the spin-up
        # Maybe I should permit a more sophisticated spin-up at some point!
        self.n_spin_up_cycles = self.cfg['general']['spin_up_cycles']

        # Resolve enforce_water_balance: keyword argument takes precedence over YAML,
        # which defaults to 'water-year' if the key is absent.
        if enforce_water_balance is None:
            enforce_water_balance = self.cfg['general'].get('enforce_water_balance',
                                                            'water-year')
        # Normalize legacy boolean values from old YAML configs or API calls.
        if enforce_water_balance is True:
            enforce_water_balance = 'water-year'
        elif enforce_water_balance is False:
            enforce_water_balance = 'none'
        self.enforce_water_balance = enforce_water_balance

        # Fraction of positive daily recharge that bypasses the reservoir
        # cascade and exits directly as runoff.  Conceptually inspired by
        # Hortonian (infiltration-excess) overland flow, but at a daily
        # timestep rainfall intensity is unavailable, so the fraction cannot
        # be a rigorous physical representation -- except in extreme events
        # where intense rainfall dominates the daily total.  In practice it
        # acts as a calibrated fast-bypass fraction, off by default.
        self.direct_runoff_fraction = self.cfg['general'].get(
            'direct_runoff_fraction', 0.0)

        # Initial conditions if resuming from prior run. Each sub-catchment
        # takes its initial SWE from its own initial_conditions block; the
        # legacy single-cascade config uses the top-level initial_conditions.
        if self.has_snowpack:
            if 'sub_catchments' in self.cfg:
                for sc, sc_cfg in zip(self.sub_catchments,
                                      self.cfg['sub_catchments']):
                    sc.snowpack.Hwater = sc_cfg.get('initial_conditions', {}).get(
                        'snowpack__mm_SWE', 0.0)
            else:
                self.sub_catchments[0].snowpack.Hwater = (
                    self.cfg['initial_conditions']['snowpack__mm_SWE'])
        # Reservoir H0 values are set during sub-catchment construction.

        # Enforce the daily timestep. MNiShed is a daily model by design:
        # degree-day snowmelt, Thornthwaite ET, and linear reservoir drainage
        # are all daily-scale parameterisations.
        if (self.hydrodata['Date'].diff()[1:] == pd.Timedelta('1 day')).all():
            self.dt = 1.
        else:
            raise ValueError(
                "MNiShed requires a continuous daily time series "
                "(exactly 1-day intervals throughout). Sub-daily or "
                "irregular timesteps are not supported; the physical "
                "parameterisations (degree-day snowmelt, Thornthwaite ET, "
                "linear reservoir drainage) are daily-scale approximations."
            )

        # Compute specific discharge from data
        self.hydrodata['Specific Discharge [mm/day]'] = (
            self.hydrodata['Discharge [m^3/s]']
            / (self.drainage_basin_area__km2*1E3) * 86400)

        # Create columns for model output
        self.hydrodata['Specific Discharge (modeled) [mm/day]'] = pd.NA
        self.hydrodata['Snowpack (modeled) [mm SWE]'] = pd.NA
        self.hydrodata['Subsurface storage (modeled total) [mm]'] = pd.NA
        self._store_depths = False

        # Start out at first timestep
        # Could modify this to pick up a run in the middle
        # Or start at the beginning of a water year
        # for example
        self._timestep_i = self.hydrodata.index[0]

        # Carry-over of any water deficit from the previous timestep that the
        # deepest reservoir could not satisfy (ET > P + all storage).  This is
        # the unpaid debt passed forward one step; distinct from
        # Reservoir.H_deficit and Snowpack.H_deficit, which are per-timestep
        # only. Stored per sub-catchment; self.H_deficit_carry reads the
        # area-weighted sum.
        for sc in self.sub_catchments:
            sc.H_deficit_carry = 0.

        # Compute the water years based on the input month for
        # water-year rollover
        self.compute_water_year()

        # Compute ET, optionally scaling to close the water balance.
        self.global_et_multiplier = 1.0   # default; overwritten below if global mode
        if self.enforce_water_balance == 'global':
            self.compute_global_ET_multiplier()
        elif self.enforce_water_balance == 'water-year':
            self.compute_ET_multiplier()
        elif self.et_method == 'ThornthwaiteChang2019':
            warnings.warn(
                "enforce_water_balance='none' with ThornthwaiteChang2019: Thornthwaite ET "
                "will not be rescaled to close the water balance. "
                "Thornthwaite ET carries large systematic biases; omitting "
                "the correction is likely to produce significant mass-balance "
                "errors. Consider enforce_water_balance='water-year' or 'global', "
                "or supply measured ET via evapotranspiration_method: datafile.",
                stacklevel=2,
            )
        self.compute_ET()

        # Model spin-up, if requested
        for _ in range(self.n_spin_up_cycles):
            self.run()  # Spin-up; run() resets _timestep_i each call

    def compute_water_year(self):
        """
        Assign a water-year label to each row in self.hydrodata.

        Adds a 'Water Year' column. A water year begins on
        water_year_start_month and is labelled by the calendar year in
        which it ends. For example, with a start month of October (USGS
        convention), 1 Oct 2020 – 30 Sep 2021 is water year 2021.

        When water_year_start_month is 1 (January), the water year equals
        the calendar year and no offset is applied.
        """
        self.hydrodata['Water Year'] = pd.DatetimeIndex(self.hydrodata['Date']).year
        if self.water_year_start_month > 1:
            self.hydrodata['Water Year'] += (
                pd.DatetimeIndex(self.hydrodata['Date']).month
                >= self.water_year_start_month
            )

    def compute_global_ET_multiplier(self, start=None, end=None):
        """
        Compute a single ET scale factor to close the long-term water balance.

        Computes scale = sum(P - Q_obs) / sum(ET_raw) over all days where
        P, ET, and Q_obs are all finite, and stores it as self.global_et_multiplier.
        Unlike compute_ET_multiplier(), which fits one multiplier per water
        year, this uses a single ratio and does not overfit interannual
        variability. Appropriate when enforce_water_balance is set to 'global'.

        Parameters
        ----------
        start : str or datetime-like, optional
            Start of the window used to compute the multiplier (inclusive).
            If None, uses the full record.
        end : str or datetime-like, optional
            End of the window used to compute the multiplier (inclusive).
            If None, uses the full record.
        """
        if self.et_method == 'datafile':
            _raw_ET = np.asarray(self.hydrodata['Evapotranspiration [mm/day]'])
        else:
            _raw_ET = np.asarray(self.evapotranspiration_Chang2019())

        # Mask to days where all three water-balance terms are finite so that
        # P_sum, ET_sum, and Q_sum are computed over the same day-set.
        mask = (self.hydrodata['Specific Discharge [mm/day]'].notna()
                & self.hydrodata['Precipitation [mm/day]'].notna()
                & np.isfinite(_raw_ET))
        if start is not None:
            mask &= self.hydrodata['Date'] >= pd.Timestamp(start)
        if end is not None:
            mask &= self.hydrodata['Date'] <= pd.Timestamp(end)
        P_sum  = float(self.hydrodata.loc[mask, 'Precipitation [mm/day]'].sum())
        Q_sum  = float(self.hydrodata.loc[mask, 'Specific Discharge [mm/day]'].sum())
        ET_sum = float(_raw_ET[mask].sum())

        if ET_sum <= 0:
            warnings.warn("Global ET sum is zero or negative; multiplier set to 1.0.",
                          UserWarning, stacklevel=2)
            self.global_et_multiplier = 1.0
        else:
            self.global_et_multiplier = (P_sum - Q_sum) / ET_sum
            if self.global_et_multiplier <= 0:
                warnings.warn(
                    f"Global ET multiplier = {self.global_et_multiplier:.3f} (≤ 0); "
                    "P − Q ≤ 0 over the water-balance window. Setting to 1.0.",
                    UserWarning, stacklevel=2)
                self.global_et_multiplier = 1.0

    def compute_ET_multiplier(self):
        """
        Compute per-water-year ET scaling factors to enforce water balance.

        For each water year, computes the ratio of required ET (P - Q) to
        measured or computed ET, and stores this as 'ET multiplier' in
        self.hydrodata_WY_means. This multiplier is later applied in
        compute_ET() to scale ET so that P - Q - ET = 0 over each water year.
        """
        # Originally used "sum", but then used "mean" so the headers would
        # still be sensible
        self.hydrodata_WY_means = self.hydrodata.groupby(
            self.hydrodata['Water Year']).mean(numeric_only=True)
        # Not needed, but no real harm in calculating
        self.hydrodata_WY_means['Runoff ratio'] = (
            self.hydrodata_WY_means['Specific Discharge [mm/day]']
            / self.hydrodata_WY_means['Precipitation [mm/day]'])
        _ET_required = -(self.hydrodata_WY_means['Specific Discharge [mm/day]'] -
                         self.hydrodata_WY_means['Precipitation [mm/day]'])
        self.hydrodata_WY_means['ET multiplier'] = (
            _ET_required / self.hydrodata_WY_means['Evapotranspiration [mm/day]'])

        _bad_wy = self.hydrodata_WY_means.index[
            self.hydrodata_WY_means['ET multiplier'] <= 0]
        if len(_bad_wy) > 0:
            warnings.warn(
                f"ET multiplier <= 0 in water year(s) {list(_bad_wy)}. "
                "Annual discharge exceeds precipitation for those years; "
                "scaled ET will be zero or negative (water-generating). "
                "Check gauge data or consider removing those years.",
                stacklevel=2,
            )

    def compute_ET(self):
        """
        Build the ET time series used in the model.

        Obtains raw daily ET from the input data file or the Thornthwaite–Chang
        2019 equation (see evapotranspiration_Chang2019()). Five modes:

        et_scale (default 1.0) is applied as a universal final multiplier in
        every mode.  Four modes determine the base ET before et_scale:

        1. et_water_stress=True (or et_reservoir_draw + 'none'): base = raw ET.
           et_scale is then the sole WB-closure mechanism.
        2. enforce_water_balance='water-year' (stress off): base = raw ET ×
           per-water-year multiplier so P - Q - ET ≈ 0 each year (before et_scale).
        3. enforce_water_balance='global' (stress off): base = raw ET × single
           full-record multiplier from sum(P - Q_obs) / sum(ET_raw).
        4. enforce_water_balance='none' (stress off): base = raw ET.

        In modes 2–4, et_scale=1.0 (the default) leaves existing behaviour
        unchanged.  Set et_scale ≠ 1 to add a calibrated offset from the
        water-balance multiplier (e.g. to capture decade-specific land-cover
        or climate sensitivity); this relaxes exact WB closure.

        The result is stored as 'ET for model [mm/day]' in self.hydrodata.
        """
        if self.et_method == 'datafile':
            _raw_ET = self.hydrodata['Evapotranspiration [mm/day]']
        elif self.et_method == 'ThornthwaiteChang2019':
            _raw_ET = self.evapotranspiration_Chang2019()
        else:
            raise ValueError('evapotranspiration_method must be "datafile" or '+
                             '"ThornthwaiteChang2019".')

        if self.use_et_water_stress or (self.use_et_reservoir_draw and
                                         self.enforce_water_balance == 'none'):
            # In stress modes, et_scale is the sole WB-closure mechanism;
            # no mode-specific multiplier is applied before it.
            _et_mode = np.asarray(_raw_ET)
        elif self.enforce_water_balance == 'global':
            _et_mode = np.asarray(_raw_ET) * self.global_et_multiplier
        elif self.enforce_water_balance == 'water-year':
            # Merge per-water-year multiplier into hydrodata, then apply.
            # Drop any previous 'ET multiplier' column first so that calling
            # compute_ET() more than once (e.g. after a module flag override)
            # does not produce duplicate _x/_y suffix columns from pandas merge.
            # Use .to_numpy() to multiply by position rather than pandas index
            # so that any index reset from the merge cannot silently misalign rows.
            if 'ET multiplier' in self.hydrodata.columns:
                self.hydrodata = self.hydrodata.drop(columns=['ET multiplier'])
            self.hydrodata = self.hydrodata.merge(
                self.hydrodata_WY_means['ET multiplier'],
                on='Water Year')
            nan_wy = (self.hydrodata_WY_means['ET multiplier']
                          .index[self.hydrodata_WY_means['ET multiplier'].isna()]
                          .tolist())
            if nan_wy:
                warnings.warn(
                    f"ET multiplier is NaN for water year(s) {nan_wy} "
                    f"(no discharge observations). Raw ET used for those years "
                    f"(enforce_water_balance ineffective).",
                    UserWarning, stacklevel=2,
                )
                self.hydrodata['ET multiplier'] = (
                    self.hydrodata['ET multiplier'].fillna(1.0))
            _et_mode = np.asarray(_raw_ET) * self.hydrodata['ET multiplier'].to_numpy()
        else:
            _et_mode = np.asarray(_raw_ET)

        # et_scale (default 1.0) is a universal calibration multiplier applied
        # after the mode-specific water-balance correction.  Under stress modes
        # it is the sole correction (mode factor = raw ET); under other modes it
        # provides an additional degree of freedom — e.g. for decade-specific
        # land-cover or climate sensitivity — at the cost of exact WB closure.
        self.hydrodata['ET for model [mm/day]'] = _et_mode * self.et_scale

    def _et_stress_factor(self, sub_catchment):
        """
        Water-availability multiplier applied to potential ET each time step.

        Returns 1 - exp(-H_shallow / H0), where H_shallow is the current water
        depth in the sub-catchment's shallowest reservoir and H0 is its PDM
        characteristic storage depth (pdm_H0).  The multiplier is zero when the
        reservoir is empty and approaches 1 as it fills, so actual ET =
        potential ET * factor.

        Returns 1.0 (no stress) when et_water_stress is disabled or when the
        shallow reservoir has no pdm_H0 set.

        Parameters
        ----------
        sub_catchment : SubCatchment
            The sub-catchment whose shallow reservoir governs the stress.
        """
        if not self.use_et_water_stress:
            return 1.0
        H0 = sub_catchment.reservoirs[0].pdm_H0
        if H0 is None:
            return 1.0
        return 1.0 - np.exp(-max(sub_catchment.reservoirs[0].Hwater, 0.0) / H0)

    def _draw_et_from_reservoirs(self, ET_pot, sub_catchment=None):
        """
        Remove actual ET from a sub-catchment's reservoirs after drainage.

        ET_pot is partitioned between reservoir 0 (shallow) and reservoir 1
        (soil) by et_alpha.  Each draw is capped at available storage so
        Hwater never goes negative; unmet demand is implicitly lost as
        water-stress reduction of actual ET.  A negative ET_pot represents
        condensation/dew (a water input); if that input would push a
        reservoir above its Hmax cap, the surplus is shed directly to
        runoff rather than stored above the cap.

        Parameters
        ----------
        ET_pot : float
            Potential ET for this time step [mm/day], already scaled by
            et_scale or the global water-balance multiplier.
        sub_catchment : SubCatchment, optional
            The sub-catchment to draw from. Defaults to the first
            sub-catchment (the whole basin in the single-sub-catchment case).

        Returns
        -------
        float
            Condensation surplus [mm] shed to runoff because it overtopped
            a reservoir's Hmax.  Zero in the usual (ET removal) case; the
            caller adds it to the day's discharge.
        """
        if sub_catchment is None:
            sub_catchment = self.sub_catchments[0]
        reservoirs = sub_catchment.reservoirs
        fractions = [self.et_alpha, 1.0 - self.et_alpha]
        excess = 0.0
        for i, frac in enumerate(fractions):
            if i >= len(reservoirs):
                break
            demand = frac * ET_pot
            H = reservoirs[i].Hwater
            if i == 1 and self.wp_soil > 0.0:
                # Soil reservoir: scale demand by fraction of catchment above WP.
                if self.wp_soil_sigma > 0.0:
                    # Gaussian CDF: spatially variable WP — σ→0 recovers hard threshold.
                    f_avail = 0.5 * (1.0 + math.erf((H - self.wp_soil)
                                                     / (self.wp_soil_sigma * math.sqrt(2.0))))
                    demand = demand * f_avail
                else:
                    # Hard threshold: no ET extraction below wp_soil.
                    if H <= self.wp_soil:
                        continue
                    H = H - self.wp_soil   # only water above WP is available
            actual = min(demand, max(0.0, H))
            reservoirs[i].Hwater -= actual
            # Condensation (negative ET) that overtops the cap runs off.
            if reservoirs[i].Hwater > reservoirs[i].Hmax:
                excess += reservoirs[i].Hwater - reservoirs[i].Hmax
                reservoirs[i].Hwater = reservoirs[i].Hmax
        return excess

    def _compute_snowpack(self, time_step, sub_catchment):
        """
        Update a sub-catchment's snowpack for one timestep; return excess melt.

        Sets temperature, applies net water input, then calls melt() with
        the raw precipitation so rain-on-snow sensible heat is included.
        Updates the sub-catchment's H_deficit_carry from the snowpack before
        returning.

        Parameters
        ----------
        time_step : int
            Current row index in self.hydrodata.
        sub_catchment : SubCatchment
            The sub-catchment whose snowpack is advanced.

        Returns
        -------
        excess_dd : float
            Leftover melt energy [°C·day] after SWE is fully depleted.
            Pass to _update_fgi() to credit toward frozen-soil thawing.
        """
        snowpack = sub_catchment.snowpack
        T = self.hydrodata['Mean Temperature [C]'][time_step]
        P = self.hydrodata['Precipitation [mm/day]'][time_step]
        snowpack.set_temperature(T)
        if self.use_et_reservoir_draw:
            _sp_recharge = P + sub_catchment.H_deficit_carry
        else:
            _sp_recharge = (P
                            - self.hydrodata['ET for model [mm/day]'][time_step]
                            * self._et_stress_factor(sub_catchment)
                            + sub_catchment.H_deficit_carry)
        snowpack.recharge(_sp_recharge)
        excess_dd = snowpack.melt(self.dt,
                                  P=(P if self.use_rain_on_snow else 0.0))
        sub_catchment.H_deficit_carry = snowpack.H_deficit
        return excess_dd

    def _update_fgi(self, time_step, sub_catchment, excess_dd=0.0):
        """
        Update the frozen ground index; flag top reservoir as frozen if needed.

        FGI(t) = max(0, A_t · FGI(t-1) - T_eff - excess_dd)
          A_t   = daily FGI decay coefficient (climate-dependent; see below)
          T_eff = T_mean · exp(-snow_insulation_k · SWE)
          T_mean < 0  → FGI rises  (freezing degree-days accumulate)
          T_mean > 0  → FGI falls  (warm air thaws)
          excess_dd   → additional thaw credited from leftover snowmelt
                        energy [°C·day] = leftover mm SWE / melt_factor

        DTR-based decay (when Minimum/Maximum Temperature columns present):
          On days entirely below freezing (T_max ≤ 0), A_t = 1.0 — no
          sub-daily thawing, FGI accumulates unimpeded. When the diurnal
          cycle straddles 0°C (T_min < 0 < T_max), a fraction of the day
          is above freezing and drives partial thawing:

              f_above = T_max / (T_max - T_min)   [linear diurnal cycle]
              A_t     = 1 - (1 - fgi_decay_coeff) · f_above

          When f_above → 1, A_t → fgi_decay_coeff (M&B 0.97 maximum decay).
          When f_above = 0, A_t = 1.0 (no passive decay).
          This naturally gives near-unity A for continental climates (where
          T_max rarely crosses 0 during cold spells) and M&B-like decay for
          maritime climates with frequent freeze-thaw oscillations.

        Fallback (T_min / T_max absent): A_t = fgi_decay_coeff (constant,
          original M&B behaviour).

        When FGI exceeds fdd_threshold, the top reservoir's f_to_discharge
        is set to 1.0 so all drainage becomes direct runoff, simulating
        frozen-soil blockage of deep infiltration.

        References
        ----------
        Molnau & Bissell (1983) https://westernsnowconference.org/bibliography/1983Molnau.pdf
            (Western Snow Conference proceedings — original source for the
            FGI formulation and the exponential snow insulation factor)
        Shanley & Chalmers (1999) doi:10.1002/(SICI)1099-1085(199909)13:12/13
            <1843::AID-HYP879>3.0.CO;2-G
        Dunne & Black (1971) doi:10.1029/WR007i005p01160
        Snow insulation parameterisation (exponential form, original: M&B 1983):
            LISFLOOD: van der Knijff et al. (2010) doi:10.1080/02626660902852568
            GSSHA: Downer & Ogden (2004) doi:10.1061/(ASCE)1084-0699(2004)9:3(254)

        Parameters
        ----------
        time_step : int
            Current row index in self.hydrodata.
        sub_catchment : SubCatchment
            The sub-catchment whose FGI and top reservoir are updated.
        excess_dd : float, optional
            Degree-day equivalent of leftover melt energy from
            _compute_snowpack() [°C·day]. Reduces FGI alongside air
            temperature. Default 0 (temperature-only FGI, per Molnau &
            Bissell). Computed as leftover_mm_SWE / melt_factor, where
            melt_factor (mm SWE °C⁻¹ day⁻¹) converts the residual melt
            depth back to the degree-day units that the FGI operates in.
            See Snowpack.melt(). Not scaled by snow insulation because
            meltwater delivers heat directly to the soil surface.

        Returns
        -------
        f0 : float
            Calibrated f_to_discharge of the top reservoir, saved before any
            frozen-ground override. Restore it after the discharge loop.
        """
        f0 = sub_catchment.reservoirs[0].f_to_discharge
        if not self.use_frozen_ground or np.isinf(self.fdd_threshold):
            return f0

        if 'Mean Temperature [C]' not in self.hydrodata.columns:
            raise ValueError(
                "fdd_threshold is set but 'Mean Temperature [C]' is missing "
                "from the input data. FGI requires temperature data."
            )
        _swe  = (sub_catchment.snowpack.Hwater
                 if sub_catchment.snowpack is not None else 0.0)
        T     = self.hydrodata['Mean Temperature [C]'][time_step]
        T_eff = T * np.exp(-self.snow_insulation_k * _swe)

        if self._has_trange:
            T_max = self.hydrodata['Maximum Temperature [C]'][time_step]
            T_min = self.hydrodata['Minimum Temperature [C]'][time_step]
            DTR   = T_max - T_min
            if DTR > 0 and T_max > 0:
                f_above = min(1.0, T_max / DTR)
                A_t = 1.0 - (1.0 - self.fgi_decay_coeff) * f_above
            else:
                A_t = 1.0  # entirely below freezing; no sub-daily thawing
        else:
            A_t = self.fgi_decay_coeff

        sub_catchment.fgi = max(0.0, A_t * sub_catchment.fgi - T_eff - excess_dd)
        if sub_catchment.fgi > self.fdd_threshold:
            sub_catchment.reservoirs[0].f_to_discharge = 1.0
        return f0

    @staticmethod
    def _sub_catchment_storage(sub_catchment):
        """
        Total subsurface storage [mm] in a sub-catchment: reservoir water plus
        any tile-drain sub-reservoir water.
        """
        return (sum(r.Hwater for r in sub_catchment.reservoirs)
                + sum(r.tile_res.Hwater for r in sub_catchment.reservoirs
                      if r.tile_res is not None))

    def _advance_sub_catchment(self, time_step, sub_catchment):
        """
        Advance one sub-catchment by a single time step.

        Runs the snowpack, frozen-ground index, the vertical reservoir cascade,
        and reservoir-draw ET for this sub-catchment, mutating its state in
        place. Forcing is read from the basin-level columns of self.hydrodata.

        Parameters
        ----------
        time_step : int
            Row index in self.hydrodata for this step.
        sub_catchment : SubCatchment
            The sub-catchment to advance.

        Returns
        -------
        qi : float
            Per-unit-area specific discharge from this sub-catchment [mm/day],
            unweighted by area fraction.
        flux_direct, flux_tile, flux_multipath : float
            Per-unit-area flux-partition components [mm/day]: direct runoff,
            tile drainage, and multipath drainage.
        """
        reservoirs = sub_catchment.reservoirs
        excess_dd = (self._compute_snowpack(time_step, sub_catchment)
                     if self.has_snowpack else 0.0)
        f0 = self._update_fgi(time_step, sub_catchment, excess_dd)

        qi = 0.0
        for i in range(len(reservoirs)):
            if i == 0:
                if self.has_snowpack:
                    _recharge = (sub_catchment.snowpack.H_infiltrated
                                 + sub_catchment.H_deficit_carry)
                else:
                    if self.use_et_reservoir_draw:
                        _recharge = (
                            self.hydrodata['Precipitation [mm/day]'][time_step] +
                            sub_catchment.H_deficit_carry)
                    else:
                        _recharge = (
                            self.hydrodata['Precipitation [mm/day]'][time_step] -
                            self.hydrodata['ET for model [mm/day]'][time_step]
                            * self._et_stress_factor(sub_catchment) +
                            sub_catchment.H_deficit_carry)
                # Hortonian-inspired bypass: fraction exits without entering reservoirs.
                _q_direct = (max(0.0, _recharge) * self.direct_runoff_fraction
                             if self.use_direct_runoff else 0.0)
                qi += _q_direct
                reservoirs[i].recharge(_recharge - _q_direct)
            else:
                # Let water infiltrate to lower layers effectively
                # instantaneously; this isn't quite realistic, but
                # should be a simpler approach for parameter calibration
                # (Plus, this is just the water that did exit that above
                # container, which is already free to discharge, so this
                # seems more self-consistent.)
                # The amount of infiltrated water from above could be
                # negative; this represents ET in excess of what the
                # unsaturated zone ("soil zone"; top reservoir) holds.
                # Deeper loss of water could be due to plants tapping into
                # groundwater, direct lake evaporation, etc. -- or related
                # to this model not being physical or distributed, so just
                # needing to balance mass.
                reservoirs[i].recharge(
                    reservoirs[i-1].H_to_next
                    + reservoirs[i-1].H_deficit)
            H_next = (reservoirs[i + 1].Hwater
                      if i + 1 < len(reservoirs) else None)
            reservoirs[i].discharge(self.dt, H_next=H_next)
            qi += reservoirs[i].H_discharge

        reservoirs[0].f_to_discharge = f0
        sub_catchment.H_deficit_carry = reservoirs[-1].H_deficit

        flux_direct    = _q_direct
        flux_tile      = sum(r.H_tile for r in reservoirs)
        flux_multipath = sum(r.H_multipath for r in reservoirs)

        if self.use_et_reservoir_draw:
            qi += self._draw_et_from_reservoirs(
                self.hydrodata['ET for model [mm/day]'][time_step],
                sub_catchment)

        return qi, flux_direct, flux_tile, flux_multipath

    def update(self, time_step=None):
        """
        Advance the model by one time step.

        Advances each parallel sub-catchment (snowpack, frozen ground, and its
        reservoir cascade from shallowest to deepest) and aggregates them to
        basin means, weighting per-unit-area quantities by sub-catchment area
        fraction. Stores modeled specific discharge, snowpack SWE, and total
        subsurface storage in self.hydrodata for the current time step. Part of
        the CSDMS Basic Model Interface.

        Parameters
        ----------
        time_step : int, optional
            Index into self.hydrodata for the time step to update. If None,
            uses and then increments the internal counter self._timestep_i.
        """

        if time_step is None:
            time_step = self._timestep_i
            # Advance internal variable if external time step is not selected.
            # This should be a different variable and therefore not
            # modify the value of "time_step" by reference.
            self._timestep_i += 1

        # Skip timesteps with missing forcing: leave reservoir states unchanged
        # and record NaN output so the scoring mask excludes these days.
        _P  = self.hydrodata['Precipitation [mm/day]'][time_step]
        _ET = self.hydrodata['ET for model [mm/day]'][time_step]
        _skip = not np.isfinite(_P) or not np.isfinite(_ET)
        if not _skip and self.has_snowpack:
            _skip = not np.isfinite(
                self.hydrodata['Mean Temperature [C]'][time_step])
        if _skip:
            self.hydrodata.at[time_step,
                              'Specific Discharge (modeled) [mm/day]'] = np.nan
            if self.has_snowpack:
                self.hydrodata.at[time_step, 'Snowpack (modeled) [mm SWE]'] = sum(
                    sc.area_fraction * sc.snowpack.Hwater
                    for sc in self.sub_catchments)
            self.hydrodata.at[time_step,
                              'Subsurface storage (modeled total) [mm]'] = sum(
                sc.area_fraction * self._sub_catchment_storage(sc)
                for sc in self.sub_catchments)
            self._flux_direct_runoff = np.nan
            self._flux_tile = np.nan
            self._flux_multipath = np.nan
            return

        qi             = 0.0
        swe            = 0.0
        storage        = 0.0
        flux_direct    = 0.0
        flux_tile      = 0.0
        flux_multipath = 0.0
        for sc in self.sub_catchments:
            a = sc.area_fraction
            qi_sc, dr_sc, tile_sc, mp_sc = self._advance_sub_catchment(
                time_step, sc)
            qi             += a * qi_sc
            if self.has_snowpack:
                swe        += a * sc.snowpack.Hwater
            storage        += a * self._sub_catchment_storage(sc)
            flux_direct    += a * dr_sc
            flux_tile      += a * tile_sc
            flux_multipath += a * mp_sc

        # Record per-step flux-partition components (mm/day) for diagnostics
        # and BMI coupling.  These are already reflected in the total
        # specific discharge above; recording them changes no result.
        self._flux_direct_runoff = flux_direct
        self._flux_tile          = flux_tile
        self._flux_multipath     = flux_multipath

        self.hydrodata.at[time_step, 'Specific Discharge (modeled) [mm/day]'] = qi
        if self.has_snowpack:
            self.hydrodata.at[time_step, 'Snowpack (modeled) [mm SWE]'] = swe
        self.hydrodata.at[time_step,
                          'Subsurface storage (modeled total) [mm]'] = storage
        if self._store_depths:
            for _res, _col in zip(self.reservoirs, self._depth_column_names()):
                self.hydrodata.at[time_step, _col] = _res.Hwater

    def evapotranspiration_Chang2019(self, Tmax=None, Tmin=None, photoperiod=None,
                                    k=0.69):
        """
        Modified daily Thornthwaite ET₀ equation.

        Chang et al. (2019), Eq. 1–4. https://doi.org/10.1002/ird.2309

        Parameters
        ----------
        Tmax : array-like
            Daily maximum temperature (°C).
        Tmin : array-like
            Daily minimum temperature (°C).
        photoperiod : array-like
            Photoperiod N (hours), computed from latitude and Julian day
            per Allen et al. (1998), Eqs. 2–4 of Chang et al. (2019).
        k : float
            Calibration coefficient in the T_ef formula. Default 0.69,
            recommended by Pereira & Pruitt (2004) for daily ET₀
            (https://doi.org/10.1016/j.agrformet.2004.01.005).
            Use 0.72 for monthly ET₀ per Camargo et al. (1999).

        Returns
        -------
        ET0 : array-like
            Daily reference evapotranspiration (mm day⁻¹).
        """
        if Tmax is None:
            Tmax = self.hydrodata['Maximum Temperature [C]']
        if Tmin is None:
            Tmin = self.hydrodata['Minimum Temperature [C]']
        if photoperiod is None:
            photoperiod = self.hydrodata['Photoperiod [hr]']

        Tef = 0.5 * k * (3 * Tmax - Tmin)
        C = photoperiod / 360.

        quadratic  = C * (-415.85 + 32.24 * Tef - 0.43 * Tef**2)
        power_law  = 16. * C * (10. * Tef / self.Chang_I) ** self.Chang_a

        ET0 = np.where(np.isnan(Tef), np.nan,
               np.where(Tef > 26,   quadratic,
               np.where(Tef > 0,    power_law,
                                    0.)))
        return ET0

    def run(self, start=None, end=None, store_depths=False):
        """
        Advance the model through time steps in self.hydrodata.

        Resets the internal time counter before iterating, so run() is safe
        to call more than once (e.g. spin-up then main run). Captures storage
        at the start of the run for check_mass_balance(). Part of the CSDMS
        Basic Model Interface.

        Parameters
        ----------
        start : str or datetime-like, optional
            First date to simulate (inclusive). If None, starts from the
            beginning of self.hydrodata.
        end : str or datetime-like, optional
            Last date to simulate (inclusive). If None, runs to the end of
            self.hydrodata.

        Notes
        -----
        When start or end is provided, only the matching rows are stepped
        through; reservoir states at the boundaries carry in/out naturally,
        making windowed runs chainable (e.g. spin-up on pre-decade data
        followed by a scored run on the decade itself).
        """
        self._timestep_i = self.hydrodata.index[0]
        # Area-weighted basin storage at the start of this run (reservoir water
        # plus snowpack; tiles start empty), for check_mass_balance().
        self._run_initial_storage = sum(
            sc.area_fraction * (
                sum(r.Hwater for r in sc.reservoirs)
                + (sc.snowpack.Hwater
                   if (self.has_snowpack and sc.snowpack is not None) else 0.0))
            for sc in self.sub_catchments)
        if start is not None or end is not None:
            date_mask = pd.Series(True, index=self.hydrodata.index)
            if start is not None:
                date_mask &= self.hydrodata['Date'] >= pd.Timestamp(start)
            if end is not None:
                date_mask &= self.hydrodata['Date'] <= pd.Timestamp(end)
            _run_idx = self.hydrodata.index[date_mask]
        else:
            _run_idx = None

        self._store_depths = store_depths
        if store_depths:
            for _col in self._depth_column_names():
                self.hydrodata[_col] = pd.NA

        # JIT path: available whenever numba imported successfully. PDM
        # (pdm_H0) and et_water_stress are both supported by the compiled loop.
        _can_jit = _numba_available

        # Surface a slow pure-Python fallback once, so it is not silent — but
        # only for an installed-but-broken numba; a plain "numba not installed"
        # stays quiet (pure Python is the expected default without the extra).
        if not _can_jit and _numba_import_error is not None:
            _notify_jit_unavailable(
                f"numba is installed but failed to import "
                f"({_numba_import_error}), usually a NumPy/Numba version "
                "mismatch (the mnished[jit] extra pins numpy<2.3)")

        if _can_jit:
            _idx = _run_idx if _run_idx is not None else self.hydrodata.index
            _hd  = self.hydrodata.loc[_idx]

            _needs_T = self.has_snowpack or (
                self.use_frozen_ground and not np.isinf(self.fdd_threshold))
            _T    = (_hd['Mean Temperature [C]'].to_numpy(dtype=np.float64)
                     if _needs_T
                     else np.zeros(len(_idx), dtype=np.float64))
            _Tmin = (_hd['Minimum Temperature [C]'].to_numpy(dtype=np.float64)
                     if self._has_trange
                     else np.full(len(_idx), np.nan))
            _Tmax = (_hd['Maximum Temperature [C]'].to_numpy(dtype=np.float64)
                     if self._has_trange
                     else np.full(len(_idx), np.nan))

            # Flatten the sub-catchments into reservoir-index ranges (matching
            # the order of self.reservoirs) plus their area weights, and gather
            # the per-sub-catchment snowpack / FGI / carried-deficit state.
            _sc_start = np.empty(self.n_sub_catchments, dtype=np.int64)
            _sc_end   = np.empty(self.n_sub_catchments, dtype=np.int64)
            _sc_area  = np.empty(self.n_sub_catchments, dtype=np.float64)
            _off = 0
            for _k, _sc in enumerate(self.sub_catchments):
                _sc_start[_k] = _off
                _off         += len(_sc.reservoirs)
                _sc_end[_k]   = _off
                _sc_area[_k]  = _sc.area_fraction
            _snow_init = np.array(
                [sc.snowpack.Hwater
                 if (self.has_snowpack and sc.snowpack is not None) else 0.0
                 for sc in self.sub_catchments], dtype=np.float64)
            _fgi_init = np.array([sc.fgi for sc in self.sub_catchments],
                                 dtype=np.float64)
            _dc_init  = np.array([sc.H_deficit_carry
                                  for sc in self.sub_catchments], dtype=np.float64)

            _jmap = {'fraction': 0, 'leakance': 1, 'threshold': 2}
            (_Q, _SWE, _Hsub, _Hres_out, _finalH, _finalHTile, _finalSnow,
             _finalFgi, _finalDC) = _jit_run(
                _hd['Precipitation [mm/day]'].to_numpy(dtype=np.float64),
                _hd['ET for model [mm/day]'].to_numpy(dtype=np.float64),
                _T, _Tmin, _Tmax,
                _sc_start, _sc_end, _sc_area,
                np.array([r.Hwater          for r in self.reservoirs], dtype=np.float64),
                _snow_init,
                _fgi_init,
                _dc_init,
                np.array([r.recession_coeff      for r in self.reservoirs], dtype=np.float64),
                np.array([r.recession_exponent  for r in self.reservoirs], dtype=np.float64),
                np.array([r.recession_H_ref     for r in self.reservoirs], dtype=np.float64),
                np.array([r.f_to_discharge      for r in self.reservoirs], dtype=np.float64),
                np.array([_jmap[r.junction_type] for r in self.reservoirs], dtype=np.int64),
                np.array([r.leakance_R if r.leakance_R is not None else np.inf
                          for r in self.reservoirs], dtype=np.float64),
                np.array([r.H_threshold          for r in self.reservoirs], dtype=np.float64),
                np.array([r.Hmax                 for r in self.reservoirs], dtype=np.float64),
                np.array([r.f_tile               for r in self.reservoirs], dtype=np.float64),
                np.array([r.tile_res.recession_coeff if r.tile_res is not None else 1.0
                          for r in self.reservoirs], dtype=np.float64),
                np.array([r.tile_res.Hwater if r.tile_res is not None else 0.0
                          for r in self.reservoirs], dtype=np.float64),
                np.array([r.multipath_threshold if r.has_multipath else 0.0
                          for r in self.reservoirs], dtype=np.float64),
                np.array([r.multipath_timescale if r.has_multipath else 0.0
                          for r in self.reservoirs], dtype=np.float64),
                np.array([r.pdm_H0 if r.pdm_H0 is not None else 0.0
                          for r in self.reservoirs], dtype=np.float64),
                self.melt_factor if self.has_snowpack else 1.0,
                self.snow_insulation_k,
                self.fgi_decay_coeff,
                self.fdd_threshold,
                self.direct_runoff_fraction,
                self.wp_soil,
                self.wp_soil_sigma,
                self.et_alpha,
                self.dt,
                self.has_snowpack,
                self.use_frozen_ground,
                self.use_rain_on_snow,
                self.use_et_reservoir_draw,
                self.use_direct_runoff,
                self.use_et_water_stress,
                self._has_trange,
            )
            self.hydrodata.loc[_idx, 'Specific Discharge (modeled) [mm/day]'] = _Q
            if self.has_snowpack:
                self.hydrodata.loc[_idx, 'Snowpack (modeled) [mm SWE]'] = _SWE
            self.hydrodata.loc[_idx, 'Subsurface storage (modeled total) [mm]'] = _Hsub
            if store_depths:
                for _i, _col in enumerate(self._depth_column_names()):
                    self.hydrodata.loc[_idx, _col] = _Hres_out[:, _i]
            for _i, _h in enumerate(_finalH):
                self.reservoirs[_i].Hwater = float(_h)
            for _i, _r in enumerate(self.reservoirs):
                if _r.tile_res is not None:
                    _r.tile_res.Hwater = float(_finalHTile[_i])
            # Write each sub-catchment's final snowpack / FGI / carried-deficit
            # state back from the per-sub-catchment result arrays.
            for _k, _sc in enumerate(self.sub_catchments):
                if self.has_snowpack and _sc.snowpack is not None:
                    _sc.snowpack.Hwater = float(_finalSnow[_k])
                _sc.fgi             = float(_finalFgi[_k])
                _sc.H_deficit_carry = float(_finalDC[_k])

        elif _run_idx is not None:
            for i in _run_idx:
                self.update(time_step=i)
        else:
            for _ in self.hydrodata.index:
                self.update()

    def finalize(self):
        """
        Report model skill and display output plots.

        Calls compute_NSE(verbose=True) to print the Nash–Sutcliffe Efficiency
        to stdout, then calls plot() to display a time-series comparison of
        observed and modeled specific discharge. Part of the CSDMS Basic Model
        Interface.
        """
        # Goodness of fit
        # Add options to print and/or save values later
        self.compute_NSE(verbose=True)
        # Plot
        # Add flag for plotting (or not) later
        self.plot()

    def plot(self):
        """
        Display a time-series comparison of precipitation and specific discharge.

        Produces a dual-axis figure: precipitation as a bar chart on the left
        axis and both observed and modeled specific discharge as line plots on
        the right axis.
        """
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%y-%m-%d'))
        plt.xlabel('Date', fontsize=14)
        plt.xticks(rotation=45, horizontalalignment='right')
        plt.ylabel('Precipitation [mm/day]', fontsize=14, color='C0')
        plt.bar(self.hydrodata['Date'].values,
                height=self.hydrodata['Precipitation [mm/day]'].values/self.dt,
                width=1., align='center', label='Precipitation [mm/day]',
                linewidth=0, color='C0', alpha=0.5)  # C0 is the default bar-plot color
        plt.twinx()
        plt.plot(self.hydrodata['Date'].values,
                 self.hydrodata['Specific Discharge [mm/day]'].values,
                 'royalblue', label='Data', linewidth=2, alpha=0.8)
        plt.plot(self.hydrodata['Date'].values,
                 self.hydrodata['Specific Discharge (modeled) [mm/day]'].values,
                 'k', label='Model', linewidth=2, alpha=0.8)
        plt.ylim(0, plt.ylim()[-1])
        plt.legend(title='Specific Discharge', fontsize=11,
                   title_fontsize=11, labelcolor='linecolor')
        plt.ylabel('Specific Discharge [mm/day]', fontsize=14, color='0.3')
        plt.tight_layout()
        plt.show()

    def check_mass_balance(self, time_step=None):
        """
        Compute the mass-balance discrepancy at a given time step.

        Compares cumulative inputs (P - ET) from the start of the record
        through time_step with cumulative outputs (discharge) plus current
        storage (snowpack + subsurface reservoirs) and any carried-over
        deficit. Returns the excess mass still in the model; a value near
        zero indicates good mass conservation.

        Parameters
        ----------
        time_step : int, optional
            Row index in self.hydrodata at which to evaluate the balance.
            Defaults to the last row.

        Returns
        -------
        excess_mass_in_model : float
            Excess water remaining in the model budget (mm). Should be ~0
            for a mass-conserving run.
        """
        if time_step is None:
            time_step = self.hydrodata.index[-1]
        # Additions equals discharge out; set up this way, and can check.
        total_additions = \
            self.hydrodata['Precipitation [mm/day]'][:time_step+1].sum() \
            - self.hydrodata['ET for model [mm/day]'][:time_step+1].sum()
        # Storage reservoirs; snowpack is 0 when not simulated
        snow_storage = (self.hydrodata['Snowpack (modeled) [mm SWE]'][time_step]
                        if self.has_snowpack else 0.)
        subsurface_storage = self.hydrodata['Subsurface storage (modeled total) [mm]'][time_step]
        # Mass removal
        outlet_discharge = self.hydrodata[
            'Specific Discharge (modeled) [mm/day]'][:time_step+1].sum()
        # Unpaid water deficit carried forward from the last timestep
        deficit = self.H_deficit_carry

        # Initial storage at the start of the last run() call (not at initialize()
        # time, since spin-up changes storage before the scored run begins).
        initial_storage = getattr(self, '_run_initial_storage', 0.0)

        # Discrepancy: inputs = outputs + ΔS, so excess ≈ 0 when mass is conserved.
        excess_mass_in_model = (outlet_discharge + subsurface_storage
                                + snow_storage - total_additions + deficit
                                - initial_storage)

        return excess_mass_in_model

    def compute_NSE(self, return_nse=True, verbose=False):
        """
        Compute the Nash–Sutcliffe Efficiency of the discharge simulation.

        Compares modeled and observed specific discharge for all time steps
        where both values are non-missing. Stores the result as self.NSE.

        Parameters
        ----------
        return_nse : bool, optional
            If True (default), return the NSE value.
        verbose : bool, optional
            If True, print the NSE value to stdout.

        Returns
        -------
        NSE : float or None
            Nash–Sutcliffe Efficiency coefficient. Returns None if
            return_nse is False. A value of 1 indicates perfect agreement;
            values below 0 indicate the model performs worse than the
            observed-mean predictor.
        """

        q_data  = self.hydrodata['Specific Discharge [mm/day]']
        q_model = self.hydrodata['Specific Discharge (modeled) [mm/day]']

        # Calculate NSE
        _realvalue = ~q_model.isna() & ~q_data.isna()
        NSE_num = np.sum((q_model[_realvalue] - q_data[_realvalue])**2)
        NSE_denom = np.sum((q_data[_realvalue] -
                            np.mean(q_data[_realvalue]))**2)
        if verbose and np.sum(~_realvalue):
            print("Excluded", np.sum(~_realvalue), "no-data points from NSE calculation")

        self.NSE = 1 - NSE_num / NSE_denom

        if verbose:
            print("NSE:", self.NSE)

        if return_nse:
            return self.NSE


def main():
    """Command-line entry point: parse ``-y``/``--configfile`` and run the model."""
    parser = argparse.ArgumentParser(
        description='Pass the configuration file path to run MNiShed.')
    parser.add_argument('-y', '--configfile', type=str,
                        help='YAML file from which all inputs are read.')

    # Parse args if anything is passed.
    # If nothing is passed, then print help and exit.
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    b = Buckets()
    try:
        b.initialize(args.configfile)
    except (FileNotFoundError, yaml.YAMLError) as e:
        # initialize() raises on a bad config; at the CLI boundary, translate
        # that into a clean stderr message and exit code 2 (no traceback).
        print(f"\n{e}\n", file=sys.stderr)
        sys.exit(2)
    b.run()
    b.finalize()


if __name__ == "__main__":
    main()
