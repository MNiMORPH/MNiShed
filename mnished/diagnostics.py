"""
Calibration diagnostics for MNiShed.

:class:`SeasonalMassBalance` turns "the seasonality is wrong" into "*which* flux,
source, or timing is wrong" — the decomposition that separates an ET-phasing
error from a routing, storage, or recession error. For a calibrated run it
tabulates, per season, the basin-mean water-balance terms (precipitation, ET,
storage change, and discharge) with **discharge split by source** (fast/event,
slow/baseflow, and lake outlet), against the observed discharge, plus a monthly
snowpack/ET climatology to locate melt timing.

Reading the table (the patterns that make it diagnostic):

* summer **ET > P** but modeled Q still over observed → the surplus is
  **slow-store release** (baseflow not receding), not an ET deficit;
* the freshet missing **and** spring ET ≈ P at melt → ET is **consuming the
  melt** (a temperature-index ET phasing problem; see the phenology coefficient);
* a per-source flow **flat across seasons** → that reservoir is not responding
  seasonally (an over-buffered lake, or a non-receding groundwater store).

The per-source split requires the run to have recorded it::

    result = run_and_score(cfg, ..., store_fluxes=True)
    smb = SeasonalMassBalance(result.buckets)
    print(smb.report())

``store_fluxes=True`` runs the pure-Python time loop (the JIT loop does not record
the partition); this is a post-calibration analysis, not a hot-loop concern.
"""

import numpy as np
import pandas as pd


class SeasonalMassBalance:
    """Per-season water-balance decomposition with discharge split by source.

    Parameters
    ----------
    buckets : Buckets
        A model run **with** ``store_fluxes=True`` (so the per-source discharge
        columns are present). Pass ``run_and_score(..., store_fluxes=True).buckets``
        or a ``Buckets`` advanced with ``run(store_fluxes=True)``.
    start, end : str or datetime-like, optional
        Restrict the analysis window (inclusive). Defaults to the full run.
    """

    SEASONS = (('DJF', (12, 1, 2)), ('MAM', (3, 4, 5)),
               ('JJA', (6, 7, 8)), ('SON', (9, 10, 11)))
    _FLUX = ('Discharge: fast [mm/day]', 'Discharge: slow [mm/day]',
             'Discharge: lake [mm/day]')

    def __init__(self, buckets, start=None, end=None):
        hd = buckets.hydrodata
        missing = [c for c in self._FLUX if c not in hd.columns]
        if missing:
            raise ValueError(
                "SeasonalMassBalance needs the per-source discharge columns; "
                "run the model with store_fluxes=True (missing: "
                f"{missing}).")

        def num(col, default=None):
            if col not in hd.columns:
                return pd.Series(default, index=hd.index)
            return pd.to_numeric(hd[col], errors='coerce')

        swe = num('Snowpack (modeled) [mm SWE]', 0.0)
        storage = num('Subsurface storage (modeled total) [mm]', np.nan) + swe
        # Basin-mean ET, area-weighted to include lake open-water evaporation:
        # land zones use the land ET demand, lakes the phenology-free open-water
        # column (= land ET when no phenology, so non-lake basins are unchanged).
        lake_frac = sum(sc.area_fraction for sc in buckets.sub_catchments
                        if getattr(sc, 'kind', 'land') == 'lake')
        land_et = num('ET for model [mm/day]')
        _ow = 'ET for model (open water) [mm/day]'
        lake_et = num(_ow) if _ow in hd.columns else land_et
        et = (1.0 - lake_frac) * land_et + lake_frac * lake_et
        df = pd.DataFrame({
            'date': pd.to_datetime(hd['Date']),
            'P':    num('Precipitation [mm/day]'),
            'ET':   et,
            'obs':  num('Specific Discharge [mm/day]'),
            'mod':  num('Specific Discharge (modeled) [mm/day]'),
            'fast': num(self._FLUX[0]),
            'slow': num(self._FLUX[1]),
            'lake': num(self._FLUX[2]),
            'SWE':  swe,
            'dS':   storage.diff(),       # daily total-storage tendency [mm/day]
        })
        if start is not None:
            df = df[df['date'] >= pd.Timestamp(start)]
        if end is not None:
            df = df[df['date'] <= pd.Timestamp(end)]
        df = df[df['mod'].notna()]        # scored days only
        df['month'] = df['date'].dt.month
        self.df = df.reset_index(drop=True)

    def seasonal_table(self):
        """Per-season means (mm/day): P, ET, dS, obs, mod, mod/obs, fast, slow,
        lake. ``dS`` is the mean daily total-storage (subsurface + snowpack)
        tendency, so ``P - ET - mod - dS`` is the model's residual closure."""
        cols = ['P', 'ET', 'dS', 'obs', 'mod', 'mod/obs', 'fast', 'slow', 'lake']
        rows = {}
        for name, months in self.SEASONS:
            s = self.df[self.df['month'].isin(months)]
            obs = s['obs'].mean()
            m = {c: s[c].mean() for c in
                 ('P', 'ET', 'dS', 'obs', 'mod', 'fast', 'slow', 'lake')}
            m['mod/obs'] = m['mod'] / obs if obs else np.nan
            rows[name] = m
        return pd.DataFrame(rows).T[cols]

    def monthly_table(self):
        """Per-month means: SWE, P, ET, obs, mod — to locate melt timing."""
        g = self.df.groupby('month')
        return pd.DataFrame({c: g[c].mean() for c in
                             ('SWE', 'P', 'ET', 'obs', 'mod')}).reindex(range(1, 13))

    def annual(self):
        """Annual-mean P, ET, obs, mod (mm/day)."""
        return {c: float(self.df[c].mean()) for c in ('P', 'ET', 'obs', 'mod')}

    def report(self):
        """A formatted multi-line text report (seasonal + monthly + annual)."""
        st = self.seasonal_table()
        lines = ["Seasonal water balance (mm/day) — discharge split by source",
                 f"{'seas':<5}{'P':>6}{'ET':>6}{'dS':>6}{'obs':>7}{'mod':>7}"
                 f"{'mod/obs':>8} |{'fast':>7}{'slow':>7}{'lake':>7}"]
        for name, r in st.iterrows():
            lines.append(
                f"{name:<5}{r['P']:>6.2f}{r['ET']:>6.2f}{r['dS']:>6.2f}"
                f"{r['obs']:>7.3f}{r['mod']:>7.3f}{r['mod/obs']:>8.2f} |"
                f"{r['fast']:>7.3f}{r['slow']:>7.3f}{r['lake']:>7.3f}")
        mt = self.monthly_table()
        lines.append("\nMonthly climatology — SWE [mm], P/ET/obs/mod [mm/day]")
        for mo, r in mt.iterrows():
            lines.append(
                f"  m{mo:>2}: SWE={r['SWE']:>6.1f}  P={r['P']:.2f}  "
                f"ET={r['ET']:.2f}  obs={r['obs']:.2f}  mod={r['mod']:.2f}")
        a = self.annual()
        lines.append(
            f"\nAnnual: P={a['P']:.2f}  ET={a['ET']:.2f}  obs={a['obs']:.3f}  "
            f"mod={a['mod']:.3f}  (P-ET={a['P'] - a['ET']:.2f})")
        return "\n".join(lines)

    def __repr__(self):
        return self.report()
