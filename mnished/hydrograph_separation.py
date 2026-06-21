"""
mnished.hydrograph_separation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Estimate reservoir timescales and initial storage depths from an observed
discharge time series, without running the hydrological model.

The pipeline has three stages:

1. **Spectral stage** (fallback / diagnostic) — fit a cascade of
   linear-reservoir transfer functions (product of Lorentzians) to the power
   spectral density of Q.  AIC selects the number of fast reservoirs.  This
   stage provides initial fast-component timescale estimates and serves as a
   fallback when the time-domain stage cannot resolve a timescale.

2. **Time-domain stage** (primary) — estimate the fast timescales directly
   from recession behaviour in the time domain (per-segment-minimum recession
   fitting).  These estimates replace the corresponding spectral values when
   available; the spectral values from stage 1 remain as fallbacks.

3. **Recession stage** — remove the fast components with a recursive digital
   filter, leaving a slowly varying residual.  Fit log(Q_residual) vs. time
   on dry-season recession segments to extract τ_karst and its initial storage.

The primary output is physically informed initial conditions and calibration
parameter priors for MNiShed, reducing spin-up requirements and improving
optimizer starting points.

References
----------
Maillet (1905) : multi-exponential recession decomposition of spring flows.
Young (1984)   : recursive estimation and modelling of nonstationary time series.
Jakeman & Hornberger (1993) : IHACRES transfer-function identification.
Kirchner (2009): catchments as simple dynamical systems — spectral approach.
Eckhardt (2005): recursive digital baseflow separation filter.
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import csd, welch
from scipy.stats import linregress


class HydrographSeparation:
    """
    Separate observed discharge into reservoir components and estimate
    timescales and initial storage depths.

    Parameters
    ----------
    Q : array-like
        Observed specific discharge [mm/day].
    dt : float
        Timestep [days]. Default 1.0.
    n_reservoirs : int or None
        Number of reservoirs. If None, selected by AIC on the spectral fit.
        The slowest reservoir (karst/deep) is always fitted by recession
        analysis regardless of this setting.
    max_reservoirs : int
        Maximum number of reservoirs to test when n_reservoirs is None.
        Default 4.
    precip : array-like or None
        Daily precipitation [mm/day], same length as Q.  Used to mask
        recession periods; if None, only discharge decline is used.
    recession_min_duration : int
        Minimum qualifying recession length [days]. Default 10.
    recession_precip_threshold : float
        Maximum daily precipitation within a recession segment [mm/day].
        Default 0.5.
    recession_antecedent_days : int
        Number of days immediately before a recession that must also be dry.
        Default 3.
    tau_deep : float or None
        Externally supplied timescale [days] for the deep (karst) reservoir.
        When provided, recession fitting is skipped for τ_karst and this
        value is used directly.  Recommended when τ_karst >> T_record, since
        a decay too slow to observe cannot be reliably estimated from the
        discharge record alone.
    tau_fast : array-like or None
        Externally supplied timescales [days] for the fast reservoirs (all
        except the deepest), e.g. from a prior calibration run.  When
        provided, spectral fitting is skipped.  Note: spectral corner periods
        are 2πτ, not τ — a 29-day reservoir has its -3 dB corner at 182 days,
        and a 466-day soil reservoir at ~8 years, which may lie near or
        beyond the record length and be unreliable from spectral fitting alone.

    Attributes (populated after fit())
    ------------------------------------
    tau : ndarray
        Reservoir timescales [days], ordered slowest to fastest.
    H0 : ndarray
        Initial storage depths [mm], same order as tau.
    n_reservoirs_fitted : int
        Number of reservoirs selected (spectral fast + 1 karst).
    aic_scores : list of float
        AIC values for k = 1 … max_reservoirs-1 spectral fits.
    tau_karst : float
        Deep-reservoir timescale [days] from recession fitting.
    """

    def __init__(self, Q, dt=1.0, n_reservoirs=None, max_reservoirs=4,
                 precip=None, recession_min_duration=10,
                 recession_precip_threshold=0.5, recession_antecedent_days=3,
                 tau_deep=None, tau_fast=None):
        """
        Parameters
        ----------
        Q : array-like
            Observed specific discharge [mm/day].
        dt : float, optional
            Timestep [days]. Default 1.0.
        n_reservoirs : int or None, optional
            Number of reservoirs. If None, selected by AIC. Default None.
        max_reservoirs : int, optional
            Maximum number of reservoirs to test when n_reservoirs is None.
            Default 4.
        precip : array-like or None, optional
            Daily precipitation [mm/day], same length as Q. Used to mask
            recession periods; if None, only discharge decline is used.
        recession_min_duration : int, optional
            Minimum qualifying recession length [days]. Default 10.
        recession_precip_threshold : float, optional
            Maximum daily precipitation within a recession segment [mm/day].
            Default 0.5.
        recession_antecedent_days : int, optional
            Number of days before a recession that must also be dry.
            Default 3.
        tau_deep : float or None, optional
            Externally supplied deep-reservoir timescale [days]; skips
            recession fitting for that component when provided.
        tau_fast : array-like or None, optional
            Externally supplied fast-reservoir timescales [days]; skips
            spectral fitting when provided.
        """
        self.Q    = np.asarray(Q, dtype=float)
        self.dt   = float(dt)
        self.n_reservoirs            = n_reservoirs
        self.max_reservoirs          = max_reservoirs
        self.precip                  = (np.asarray(precip, dtype=float)
                                        if precip is not None else None)
        self.recession_min_duration      = recession_min_duration
        self.recession_precip_threshold  = recession_precip_threshold
        self.recession_antecedent_days   = recession_antecedent_days
        self.tau_deep = float(tau_deep) if tau_deep is not None else None
        # Externally supplied fast timescales [days], shortest first.
        # When provided, skip spectral fitting and use these directly.
        self.tau_fast = (np.sort(np.asarray(tau_fast, dtype=float))
                         if tau_fast is not None else None)

        # Results
        self.tau                 = None
        self.H0                  = None
        self.n_reservoirs_fitted = None
        self.aic_scores          = None
        self.tau_karst           = None
        self.tau_karst_longterm  = None   # long-term absolute-time estimate (comparison)
        self._tau_spectral       = None   # fast timescales only (fastest first)
        self._Q_residual         = None   # Q after removing fast components
        self._Q_components_fast  = None   # list of fast-component time series
        self._psd_freqs          = None
        self._psd                = None
        self._tau_shallow_td     = None   # time-domain τ_shallow estimate
        self._tau_soil_td        = None   # time-domain τ_soil estimate

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def fit(self):
        """Run the full separation pipeline. Returns self."""
        self._compute_psd()
        self._fit_spectral_model()
        self._fit_timescales_td()
        self._separate_components()
        self._fit_recession()
        self._compute_initial_conditions()
        return self

    def get_initial_conditions(self):
        """
        Initial reservoir storage depths, ordered fastest to slowest
        (shallow first, karst last), matching the MNiShed reservoir
        index convention.

        Returns
        -------
        dict
            ``{'H0': [H0_shallow, H0_soil, ..., H0_karst]}`` in mm.
        """
        self._check_fitted()
        return {'H0': list(self.H0[::-1])}

    def get_parameter_priors(self):
        """
        Suggested initial values and bounds (in params.yml units) for the
        MNiShed log__recession_coeff_* calibration parameters.

        Bounds are ±0.5 in log10 around the estimated timescale — wide enough
        to allow the optimizer freedom, tight enough to stay physical.

        Returns
        -------
        dict
            Keyed by MNiShed parameter name.  Each value is a dict with
            ``'initial'``, ``'lower'``, ``'upper'`` in log10(days).
        """
        self._check_fitted()
        all_names = ['log__recession_coeff_shallow', 'log__recession_coeff_soil', 'log__recession_coeff_karst']
        tau_fast_to_slow = self.tau[::-1]   # shallow first, karst last
        n = len(tau_fast_to_slow)
        # Always end with karst; fill fast slots from the front.
        # For 2 reservoirs: ['shallow', 'karst'] — skip soil.
        names = all_names[:n - 1] + [all_names[-1]]
        priors = {}
        for name, tau_i in zip(names, tau_fast_to_slow):
            if not np.isfinite(tau_i):
                # For τ_soil, fall back to the time-domain rolling-minimum
                # estimate (not used for LP separation but good enough as a
                # prior hint — still biased low, so bounds are left wide).
                if name == 'log__recession_coeff_soil' and self._tau_soil_td is not None:
                    tau_i = self._tau_soil_td
                else:
                    priors[name] = None   # could not be estimated; keep params.yml defaults
                    continue
            log_tau = np.log10(tau_i)
            priors[name] = {
                'initial': round(float(log_tau),       3),
                'lower':   round(float(log_tau) - 0.5, 3),
                'upper':   round(float(log_tau) + 0.5, 3),
            }
        return priors

    def summary(self):
        """Print a brief summary of fitted timescales and initial conditions."""
        self._check_fitted()
        tau_f2s = self.tau[::-1]
        H0_f2s  = self.H0[::-1]
        labels  = ['shallow', 'soil', 'karst'] + [f'res{i}' for i in range(3, len(tau_f2s))]
        print('HydrographSeparation results')
        print(f'  n_reservoirs : {self.n_reservoirs_fitted}')
        print(f'  {"Reservoir":<10}  {"τ (days)":>12}  {"H₀ (mm)":>12}')
        print(f'  {"-"*10}  {"-"*12}  {"-"*12}')
        for label, tau_i, H0_i in zip(labels, tau_f2s, H0_f2s):
            print(f'  {label:<10}  {tau_i:>12.1f}  {H0_i:>12.1f}')
        if self.aic_scores is not None:
            print(f'\n  AIC by fast-reservoir order: '
                  f'{[round(a, 1) for a in self.aic_scores]}')
        print('\n  Time-domain estimates:')
        td_sh = (f'{self._tau_shallow_td:.1f} d'
                 if self._tau_shallow_td is not None else 'not estimated')
        print(f'    τ_shallow : {td_sh}')
        if self.n_reservoirs_fitted >= 3 or self._tau_soil_td is not None:
            td_so = (f'{self._tau_soil_td:.1f} d'
                     if self._tau_soil_td is not None
                     else 'not estimated (insufficient long dry-period recessions)')
            print(f'    τ_soil    : {td_so}')
        fallback_tau = 10.0 * len(self.Q) * self.dt
        is_fallback  = abs(self.tau_karst / fallback_tau - 1.0) < 0.01
        fb_note      = '  ← fallback: no declining trend detected' if is_fallback else ''
        print(f'\n  τ_karst  (segment minima):  {self.tau_karst:.1f} days{fb_note}')
        if self.tau_karst_longterm is not None and not is_fallback:
            agree = abs(np.log(self.tau_karst / self.tau_karst_longterm)) < 0.5
            print(f'  τ_karst  (long-term slope): {self.tau_karst_longterm:.1f} days'
                  f'  {"← agree" if agree else "← disagree — check interannual trend"}')

    # ------------------------------------------------------------------
    # Stage 1a: time-domain τ estimation (primary for fast reservoirs)
    # ------------------------------------------------------------------

    def _fit_recession_scale(self, Q_signal, min_dur, smooth_days, antecedent_days,
                              precip_threshold=None, min_tau=1.0, max_tau=None,
                              mode='within'):
        """
        Estimate a recession timescale from recession-segment analysis.

        Both current callers (τ_shallow and τ_soil in :meth:`_fit_timescales_td`)
        use ``mode='within'``.  τ_karst has its own dedicated pipeline in
        :meth:`_fit_recession`, which uses the same per-segment-minimum approach
        as ``mode='cross'`` but adds a long-term cross-check and reconciliation
        step not available here.

        Two modes:

        ``'within'`` (default) — for fast/intermediate timescales:
            Fit log(Q) vs relative time within each qualifying segment; take
            the median slope.  Captures the LOCAL recession rate, which is
            dominated by the fastest reservoir still active over that duration.
            Use short segments (5-20 d) for τ_shallow; longer segments (≥15 d)
            on the rolling-minimum envelope for τ_soil.

        ``'cross'`` — for inter-annual timescales:
            Record the minimum Q within each segment and fit log(Q_min) vs
            absolute time across all segments.  Captures the INTER-ANNUAL
            drainage trend that persists after fast reservoirs empty.
            Not called internally; available for external or diagnostic use.

        Parameters
        ----------
        Q_signal : ndarray
            Discharge signal to analyse [mm/day].
        min_dur : int
            Minimum qualifying segment length [days].
        smooth_days : float
            Moving-average window [days] for the declining test only.
        antecedent_days : int
            Days before segment start that must also satisfy the precipitation
            gate.
        precip_threshold : float or None
            Precipitation gate [mm/day].  Defaults to
            self.recession_precip_threshold.
        min_tau : float
            Reject fits with τ < min_tau [days].
        max_tau : float or None
            Reject fits with τ > max_tau [days].
        mode : {'within', 'cross'}
            Regression strategy (see above).

        Returns
        -------
        tau : float or None
            Estimated timescale [days], or None if fewer than three segments
            qualify or the slope is non-negative.
        n_segs : int
            Number of qualifying segments used.
        """
        n         = len(Q_signal)
        threshold = (precip_threshold if precip_threshold is not None
                     else self.recession_precip_threshold)

        # Smooth signal for declining test only
        w = max(1, int(round(smooth_days / self.dt)))
        if w > 1:
            kernel   = np.ones(w) / w
            Q_smooth = np.convolve(Q_signal, kernel, mode='same')
        else:
            Q_smooth = Q_signal.copy()
        Q_smooth = np.maximum(Q_smooth, 0.0)

        # Precipitation gate
        if self.precip is not None:
            dry = self.precip < threshold
            for lag in range(1, antecedent_days + 1):
                dry[lag:] &= dry[:-lag]
        else:
            dry = np.ones(n, dtype=bool)

        # Declining test on smoothed signal
        declining      = np.zeros(n, dtype=bool)
        declining[:-1] = Q_smooth[1:] < Q_smooth[:-1]
        candidate      = dry & declining

        # Collect qualifying segments
        mask = np.zeros(n, dtype=bool)
        i    = 0
        while i < n:
            if candidate[i]:
                j = i + 1
                while j < n and candidate[j]:
                    j += 1
                if j - i >= min_dur:
                    mask[i:j] = True
                i = j
            else:
                i += 1

        if not np.any(mask):
            return None, 0

        idx       = np.where(mask)[0]
        breaks    = np.where(np.diff(idx) > 1)[0] + 1
        segs      = np.split(idx, breaks)
        min_valid = max(3, min_dur // 3)

        if mode == 'within':
            slopes = []
            for seg in segs:
                Q_seg = Q_signal[seg]
                valid = Q_seg > 0
                if valid.sum() < min_valid:
                    continue
                t_rel = np.arange(valid.sum()) * self.dt
                s, _, _, _, _ = linregress(t_rel, np.log(Q_seg[valid]))
                if s < -1e-12:
                    slopes.append(s)
            n_segs = len(slopes)
            if n_segs < 3:
                return None, n_segs
            slope = float(np.median(slopes))

        else:  # 'cross': per-segment minimum vs absolute time
            t_min_list = []
            Q_min_list = []
            for seg in segs:
                Q_seg = Q_signal[seg]
                valid = Q_seg > 0
                if valid.sum() < min_valid:
                    continue
                min_pos = np.argmin(Q_seg[valid])
                abs_idx = seg[valid][min_pos]
                t_min_list.append(abs_idx * self.dt)
                Q_min_list.append(float(Q_seg[valid][min_pos]))
            n_segs = len(t_min_list)
            if n_segs < 3:
                return None, n_segs
            slope, _, _, _, _ = linregress(np.array(t_min_list),
                                            np.log(np.array(Q_min_list)))

        if slope >= -1e-12:
            return None, n_segs

        tau = -1.0 / slope
        if tau < min_tau:
            return None, n_segs
        if max_tau is not None and tau > max_tau:
            return None, n_segs

        return tau, n_segs

    def _fit_timescales_td(self):
        """
        Estimate fast-reservoir timescales via time-domain recession analysis
        and update self._tau_spectral with the results.

        Two stages of sequential peeling:

        1. **τ_shallow** — within-segment log(Q) slope on raw Q with short (≥5 day)
           segments and a 1-day smooth.  Accepts τ in [2, 200] days.

        2. **τ_soil** — per-segment minimum on Q after removing the shallow
           component (low-pass at τ_shallow).  Uses a moving-average window
           equal to τ_shallow to suppress residual fast noise, a relaxed
           precipitation threshold (2 mm/day — small rain barely recharges the
           soil), and an antecedent window ≈ τ_shallow to wait for the shallow
           reservoir to drain.  Only attempted when n_reservoirs ≥ 3.

        Time-domain estimates replace the corresponding spectral values in
        self._tau_spectral when available; spectral values serve as fallbacks
        only if they are well-separated from the shallow scale (> 5×τ_shallow).
        Diagnostic attributes _tau_shallow_td and _tau_soil_td are set
        regardless of whether they replace spectral values.
        """
        if self.tau_fast is not None:
            self._tau_shallow_td = None
            self._tau_soil_td    = None
            return

        n_fast_target = (self.n_reservoirs - 1 if self.n_reservoirs is not None
                         else self.max_reservoirs - 1)

        # ---- τ_shallow -------------------------------------------------------
        # Within-segment median slope on short (5-15 d) recession segments.
        # Short segments capture the fast shallow reservoir; deeper reservoirs
        # act as a near-constant baseline and barely affect the slope.
        tau_sh, _ = self._fit_recession_scale(
            self.Q,
            min_dur=5,
            smooth_days=1.0,
            antecedent_days=self.recession_antecedent_days,
            min_tau=2.0,
            max_tau=200.0,
            mode='within',
        )
        self._tau_shallow_td = tau_sh

        # ---- τ_soil (only when ≥ 3 reservoirs expected) ---------------------
        # Problem with LP filtering for τ_soil estimation: applying a low-pass
        # filter at τ_shallow creates a signal (Q_slow) whose FILTER IMPULSE
        # RESPONSE decays with timescale τ_shallow.  Recession segments in Q_slow
        # that start within ~3×τ_shallow days of a rain event are still in the
        # filter-decay regime, so within-segment slopes give τ_eff ≈ τ_shallow,
        # not τ_soil.  With Minnesota dry periods of only 20-35 days, most
        # segments start before the filter has fully settled.
        #
        # Fix: use a ROLLING MINIMUM over a 2×τ_shallow window.  The rolling
        # minimum extracts the lower envelope of Q (slow baseflow), cleanly
        # suppressing storm pulses without creating a lingering exponential tail.
        # The resulting Q_rmin ≈ Q_soil + Q_karst essentially everywhere it is
        # declining, with no multi-week filter-decay artefact.
        # Within-segment slopes on Q_rmin give
        #   τ_eff = (Q_soil/Q_rmin)/τ_soil + (Q_karst/Q_rmin)/τ_karst
        # Theoretically τ_eff > τ_soil (karst floor slows the apparent decay),
        # but in practice residual shallow contamination in early recession
        # segments biases τ_eff short.  The min_tau = 8×τ_shallow guard rejects
        # the worst-contaminated segments; the result is a rough lower bound on
        # τ_soil, useful as a prior but not a precise estimate.
        self._tau_soil_td = None
        if n_fast_target >= 2 and tau_sh is not None:
            n_sig  = len(self.Q)
            w_rmin = max(1, int(round(2.0 * tau_sh / self.dt)))
            Q_rmin = np.array([
                float(np.min(self.Q[max(0, t - w_rmin + 1):t + 1]))
                for t in range(n_sig)
            ])
            tau_so, _ = self._fit_recession_scale(
                Q_rmin,
                min_dur=15,
                smooth_days=tau_sh,
                antecedent_days=5,
                precip_threshold=2.0,
                min_tau=8.0 * tau_sh,
                mode='within',
            )
            self._tau_soil_td = tau_so

        # ---- update _tau_spectral --------------------------------------------
        # NOTE: _tau_soil_td is intentionally NOT added to _tau_spectral here.
        # The rolling-minimum estimate is still biased low (transition regime),
        # so including it in the LP separation would corrupt Q_residual and
        # degrade τ_karst. It is used only in get_parameter_priors().
        n_sp     = len(self._tau_spectral)
        td_taus  = []

        # Shallow position (fastest)
        if tau_sh is not None:
            td_taus.append(tau_sh)
        elif n_sp >= 1:
            td_taus.append(float(self._tau_spectral[0]))

        # Soil position (second fastest), only when ≥ 3 reservoirs requested
        if n_fast_target >= 2:
            if n_sp >= 2:
                # Accept spectral τ_soil only when it is well-separated from
                # τ_shallow — if spectral fitting collapsed both to the same
                # corner (common when seasonal peaks dominate the PSD), discard
                ref = td_taus[0] if td_taus else None
                sp_soil = float(self._tau_spectral[1])
                if ref is not None and sp_soil > 5.0 * ref:
                    td_taus.append(sp_soil)
                # else: drop to 2-reservoir initialisation

        if td_taus:
            self._tau_spectral       = np.sort(np.array(td_taus))  # shortest first
            self.n_reservoirs_fitted = len(td_taus) + 1            # +1 for karst

    # ------------------------------------------------------------------
    # Stage 1b: spectral fitting (fallback / diagnostic)
    # ------------------------------------------------------------------

    def _compute_psd(self):
        """
        Estimate the spectral target for Lorentzian fitting.

        When precipitation is provided, compute the transfer function magnitude
        |H(f)|² = |S_QP(f)|² / S_PP(f) (cross-spectrum / input PSD).  This
        cancels the colour of the precipitation spectrum and gives a clean
        estimate of the linear-reservoir cascade response regardless of whether
        the input is white noise, intermittent, or seasonal.

        Without precipitation, fall back to the output PSD S_QQ(f), which
        conflates input colour with reservoir structure but is still useful
        when timescale separation is large.
        """
        Q = self.Q.copy()
        Q[~np.isfinite(Q)] = np.nanmean(Q)
        nperseg = min(len(Q) // 4, 1024)
        kw = dict(fs=1.0 / self.dt, nperseg=nperseg,
                  window='hann', detrend='constant')

        if self.precip is not None:
            P = self.precip.copy()
            P[~np.isfinite(P)] = 0.0
            freqs, S_PP  = welch(P, **kw)
            freqs, S_QP  = csd(Q, P, **kw)
            # Transfer function magnitude squared; guard against near-zero S_PP
            safe_spp = np.where(S_PP > 0, S_PP, np.nan)
            H2 = np.abs(S_QP) ** 2 / safe_spp
            spectral = H2
        else:
            freqs, spectral = welch(Q, **kw)

        # The above uses nperseg = T_record/4, which may not resolve the
        # slow-reservoir corner.  Supplement with a full-record periodogram
        # (nperseg = N) for the low-frequency bins that are missing.  The
        # periodogram is noisier but gives the only spectral estimates below
        # 4/T_record, where the soil or karst corners may sit.
        n = len(Q)
        if nperseg < n:
            kw_full = dict(fs=1.0 / self.dt, nperseg=n,
                           window='hann', detrend='constant')
            if self.precip is not None:
                f2, S_PP2 = welch(P, **kw_full)
                f2, S_QP2 = csd(Q, P, **kw_full)
                safe2     = np.where(S_PP2 > 0, S_PP2, np.nan)
                spec2     = np.abs(S_QP2) ** 2 / safe2
            else:
                f2, spec2 = welch(Q, **kw_full)
            # Use full-record estimate for frequencies below the Welch resolution
            f_cutover = 4.0 / (n * self.dt)
            low_mask  = (f2 > 0) & (f2 < f_cutover) & np.isfinite(spec2) & (spec2 > 0)
            freqs   = np.concatenate([f2[low_mask],   freqs])
            spectral = np.concatenate([spec2[low_mask], spectral])

        # Drop DC and any non-finite bins
        mask = (freqs > 0) & np.isfinite(spectral) & (spectral > 0)
        self._psd_freqs = freqs[mask]
        self._psd       = spectral[mask]

    @staticmethod
    def _log_cascade_psd(log_freqs, log_A, *log_taus):
        """
        log PSD of a cascade of k linear reservoirs driven by white noise:
            log S(f) = log_A - Σ_i log(1 + (2π f τ_i)²)
        Parameters passed as (log_A, log_τ_1, ..., log_τ_k) so that
        scipy.optimize.curve_fit can treat them as a flat vector.
        """
        freqs = np.exp(log_freqs)
        log_S = np.full_like(freqs, log_A)
        for log_tau in log_taus:
            tau = np.exp(log_tau)
            log_S -= np.log(1.0 + (2.0 * np.pi * freqs * tau) ** 2)
        return log_S

    def _fit_k_reservoirs(self, k):
        """
        Fit k Lorentzians to log PSD.

        Returns
        -------
        taus : ndarray or None
            Fitted timescales [days], sorted longest first.
        aic : float
            AIC of the fit (lower is better).
        """
        log_f = np.log(self._psd_freqs)
        log_S = np.log(self._psd)

        # Spread initial τ guesses log-uniformly across the observable range
        T_record = len(self.Q) * self.dt
        tau_init  = np.logspace(0, np.log10(T_record / 4.0), k)
        p0 = [float(np.mean(log_S))] + list(np.log(tau_init))

        tau_lo = np.log(0.5)
        tau_hi = np.log(5.0 * T_record)
        bounds  = ([-np.inf] + [tau_lo] * k,
                   [ np.inf] + [tau_hi] * k)

        try:
            popt, _ = curve_fit(self._log_cascade_psd, log_f, log_S,
                                p0=p0, bounds=bounds, maxfev=20000)
        except (RuntimeError, ValueError):
            return None, np.inf

        taus   = np.sort(np.exp(popt[1:]))[::-1]   # longest first
        resid  = log_S - self._log_cascade_psd(log_f, *popt)
        n      = len(log_S)
        rss    = np.dot(resid, resid)
        aic    = n * np.log(rss / n) + 2.0 * (k + 1)   # k taus + amplitude
        return taus, aic

    def _fit_spectral_model(self):
        """
        Determine fast-reservoir timescales (all reservoirs except the deepest).

        If tau_fast was supplied at construction, use it directly and skip
        spectral fitting.  Spectral fitting works well when the system is a
        nearly pure cascade (small exfiltration fractions so most water
        reaches the deep reservoir) and when the corner periods 2πτ_i fall
        within the observable frequency band.  For real data with seasonal
        forcing and short records relative to τ_soil, externally supplied
        values from a prior calibration run are more reliable.
        """
        if self.tau_fast is not None:
            self._tau_spectral       = self.tau_fast          # already sorted
            self.n_reservoirs_fitted = len(self.tau_fast) + 1 # +1 for deep
            self.aic_scores          = None
            return

        if self.n_reservoirs is not None:
            n_fast_max = self.n_reservoirs - 1
        else:
            n_fast_max = self.max_reservoirs - 1

        aic_scores = []
        tau_by_k   = []
        for k in range(1, n_fast_max + 1):
            taus, aic = self._fit_k_reservoirs(k)
            aic_scores.append(aic)
            tau_by_k.append(taus)

        self.aic_scores = aic_scores

        if self.n_reservoirs is not None:
            best_k = self.n_reservoirs - 1
        else:
            best_k = int(np.argmin(aic_scores)) + 1

        taus_best = tau_by_k[best_k - 1]
        if taus_best is None:
            # curve_fit failed for the selected k; fall back to 1 fast reservoir
            for k_try in range(1, n_fast_max + 1):
                if tau_by_k[k_try - 1] is not None:
                    taus_best = tau_by_k[k_try - 1]
                    best_k    = k_try
                    break
            else:
                taus_best = np.array([1.0])
                best_k    = 1
        self._tau_spectral       = np.sort(taus_best)  # fastest first
        self.n_reservoirs_fitted = best_k + 1          # +1 for karst

    # ------------------------------------------------------------------
    # Stage 2: component separation
    # ------------------------------------------------------------------

    def _apply_lowpass(self, signal, tau):
        """
        Forward recursive low-pass filter with timescale τ [days].
        Returns the slowly varying component of signal.
        """
        alpha  = np.exp(-self.dt / tau)
        result = np.empty_like(signal)
        result[0] = signal[0]
        for t in range(1, len(signal)):
            result[t] = alpha * result[t - 1] + (1.0 - alpha) * signal[t]
        return result

    def _separate_components(self):
        """
        Peel fast components from Q one by one (fastest first).
        Each pass: apply low-pass filter at τ_i, subtract → fast component.
        Residual after all passes is the slow (karst) signal.
        """
        Q = self.Q.copy()
        Q[~np.isfinite(Q)] = np.nanmean(Q)

        components = []
        for tau_i in self._tau_spectral:        # shortest τ first
            Q_slow = self._apply_lowpass(Q, tau_i)
            components.append(Q - Q_slow)       # fast component at this τ
            Q = Q_slow

        self._Q_components_fast = components    # list: fastest component first
        self._Q_residual = np.maximum(Q, 0.0)  # slow residual for recession fit

    # ------------------------------------------------------------------
    # Stage 3: recession fitting for τ_karst
    # ------------------------------------------------------------------

    def _recession_mask(self):
        """
        Boolean mask of timesteps within qualifying recession segments:
        Q_residual monotonically declining, precipitation absent, minimum
        duration satisfied.
        """
        n = len(self._Q_residual)

        # Precipitation gate
        if self.precip is not None:
            dry = self.precip < self.recession_precip_threshold
            for lag in range(1, self.recession_antecedent_days + 1):
                dry[lag:] &= dry[:-lag]
        else:
            dry = np.ones(n, dtype=bool)

        # Discharge declining
        declining        = np.zeros(n, dtype=bool)
        declining[:-1]   = self._Q_residual[1:] < self._Q_residual[:-1]

        candidate = dry & declining

        # Enforce minimum run length
        mask = np.zeros(n, dtype=bool)
        i = 0
        while i < n:
            if candidate[i]:
                j = i + 1
                while j < n and candidate[j]:
                    j += 1
                if j - i >= self.recession_min_duration:
                    mask[i:j] = True
                i = j
            else:
                i += 1
        return mask

    def _fit_recession(self):
        """
        Estimate τ_karst using per-segment minimum-flow analysis, cross-checked
        against a long-term absolute-time regression on all recession timesteps.

        Primary method — per-segment minima:
            For each qualifying recession segment, record the minimum
            Q_residual value and its absolute time position.  The minimum
            represents the most-depleted state within that segment, after fast
            reservoirs have substantially drained.  OLS on log(Q_min) vs
            absolute time gives a slope -1/τ_karst.

        Cross-check — long-term slope:
            Fit log(Q_residual) vs absolute time across all recession timesteps.
            This averages over more data points and is more robust to soil
            contamination of the per-segment estimate (unremoved soil components
            in Q_residual shorten the per-segment estimate because per-segment
            minima track the seasonal soil cycle rather than the interannual
            karst trend).

        Reconciliation:
            When both estimates are available and they disagree by more than
            factor 2 (|log ratio| > log 2 ≈ 0.69), the LONGER τ is adopted as
            the final tau_karst.  Physically, soil contamination always biases
            the per-segment estimate SHORT, so the longer of the two is a
            better lower-bound on τ_karst.

        Both raw estimates are stored (tau_karst and tau_karst_longterm).

        If tau_deep was supplied at construction, skip fitting and use it.
        """
        if self.tau_deep is not None:
            self.tau_karst          = self.tau_deep
            self.tau_karst_longterm = self.tau_deep
            return

        mask         = self._recession_mask()
        fallback_tau = 10.0 * len(self.Q) * self.dt

        # ---- long-term slope (always computed for comparison) ----
        if np.any(mask):
            idx   = np.where(mask)[0]
            Q_dry = self._Q_residual[idx]
            pos   = Q_dry > 0
            if pos.sum() >= 5:
                slope_lt, _, _, _, _ = linregress(idx[pos] * self.dt,
                                                   np.log(Q_dry[pos]))
                self.tau_karst_longterm = (fallback_tau if slope_lt >= -1e-12
                                           else -1.0 / slope_lt)
            else:
                self.tau_karst_longterm = fallback_tau
        else:
            self.tau_karst_longterm = fallback_tau

        if not np.any(mask):
            self.tau_karst = fallback_tau
            return

        # ---- per-segment minima ----
        idx    = np.where(mask)[0]
        breaks = np.where(np.diff(idx) > 1)[0] + 1
        segs   = np.split(idx, breaks)

        t_min_list  = []
        Q_min_list  = []
        for seg in segs:
            Q_seg = self._Q_residual[seg]
            valid = Q_seg > 0
            if valid.sum() < self.recession_min_duration:
                continue
            min_pos = np.argmin(Q_seg[valid])
            abs_idx = seg[valid][min_pos]
            t_min_list.append(abs_idx * self.dt)
            Q_min_list.append(Q_seg[valid][min_pos])

        if len(t_min_list) < 3:
            self.tau_karst = self.tau_karst_longterm
            return

        t_min  = np.array(t_min_list)
        logQ_m = np.log(np.array(Q_min_list))
        slope, _, _, _, _ = linregress(t_min, logQ_m)

        tau_segmin = fallback_tau if slope >= -1e-12 else -1.0 / slope

        # Reconcile: prefer the LONGER estimate when there is strong disagreement.
        # Per-segment minimum is biased short when Q_residual still contains a
        # soil component whose seasonal drainage dominates the segment minima.
        tau_lt = self.tau_karst_longterm
        if (tau_lt < fallback_tau and
                abs(np.log(tau_segmin / tau_lt)) > np.log(2.0)):
            self.tau_karst = max(tau_segmin, tau_lt)
        else:
            self.tau_karst = tau_segmin

    # ------------------------------------------------------------------
    # Stage 4: initial conditions
    # ------------------------------------------------------------------

    def _compute_initial_conditions(self):
        """
        H₀_i = Q_i(t = 0) × τ_i for each reservoir.
        tau and H0 are stored slowest-first (karst, …, shallow).

        When fewer fast reservoirs were estimated than n_reservoirs specifies
        (e.g., τ_soil not identifiable from recession analysis), the missing
        intermediate reservoirs are padded with H₀ = 0 and τ = NaN so that
        get_initial_conditions() returns the expected number of entries.
        """
        H0_karst = float(self._Q_residual[0]) * self.tau_karst

        H0_fast = [max(float(Q_c[0]), 0.0) * tau_i
                   for Q_c, tau_i in zip(self._Q_components_fast,
                                         self._tau_spectral)]

        # Pad with NaN/0 for any fast reservoirs that could not be estimated
        n_fast_needed = (self.n_reservoirs - 1 if self.n_reservoirs is not None
                         else len(self._tau_spectral))
        n_missing = n_fast_needed - len(H0_fast)
        if n_missing > 0:
            # Insert unestimated reservoirs between the fastest and the karst.
            # tau = NaN signals "unknown"; H0 = 0 is a conservative start.
            H0_fast = [0.0] * n_missing + H0_fast
            tau_pad = [float('nan')] * n_missing + list(self._tau_spectral[::-1])
        else:
            tau_pad = list(self._tau_spectral[::-1])

        # Slowest first
        self.tau = np.array([self.tau_karst] + tau_pad)
        self.H0  = np.array([H0_karst]       + H0_fast)

    # ------------------------------------------------------------------
    # Utility
    # ------------------------------------------------------------------

    def _check_fitted(self):
        if self.tau is None:
            raise RuntimeError('Call fit() before accessing results.')
