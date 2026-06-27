# Handoff: document Thornthwaite ET's phenology/seasonality limitation

*From the Wickert2026 hydroRaVENS decadal-optimization work (Crow Wing River),
2026-06-27. Self-contained — decide for yourself what (if anything) to do here.*

> **STATUS — RESOLVED 2026-06-27.** All actions below are actioned; kept as the
> origin record of the finding. Disposition:
> - (1) docs caveat → DONE: `.. warning::` in configuration.rst + model_description
>   echo + an accept-vs-update tradeoff section (pushed).
> - (2) NDVI/phenology-driven ET → issue #26 (kept separate from the land-cover
>   DDF prior #22).
> - (3) seasonal mass-balance diagnostic → issue #25 (scoped: docs + functionality).
> - Beyond the original list: a **simple GDD phenology Kc** was implemented
>   (`phenology:` config block; pushed), and building its water-year side exposed
>   and fixed a **pre-existing Thornthwaite + water-year water-balance non-closure**
>   bug (`_demand_ET`; local, staged in [Unreleased] for a bundled patch).

## The finding (general)

`ThornthwaiteChang2019` derives ET from **temperature alone**, so it cannot
represent seasonality when the actual ET cycle is **out of phase with
temperature**. The clearest case is **cold-region forests**: canopy leaf-out
lags spring warming by weeks, so Thornthwaite ramps ET up with the spring
temperature rise while little transpiration actually occurs before leaf-out.
Early-spring ET is therefore over-estimated. In a snowmelt basin this can
**consume the spring freshet** — evaporating meltwater that should appear as
streamflow — and, because calibration then lowers `et_scale` to recover the
annual water balance, the surplus **inflates flow in the remaining seasons**.

## How it surfaced (evidence)

Calibrating a Crow Wing (north-central MN, forested) frozen-ground / two-layer
model, a residual persisted after a good multi-decade fit: **spring under,
fall over** (MAM mod/obs 0.77, SON 1.21). We first suspected snowmelt routing
through the lake. A **seasonal water-balance decomposition** (P, ET, Q, ΔStorage
by season) reframed it:

- *Not* a summer-ET deficit — summer ET (3.07 mm/d) already exceeded P (2.85).
- The spring melt arrives in April, when Thornthwaite ET ≈ April P — the freshet
  was being **evaporated**, not mis-routed.

A monthly leaf-out factor on ET (April suppressed, **May ramp**, full canopy
summer, senescing fall — verified against regional phenology: north-central MN
canopy leaf-out is mid-to-late May, ~9 days later than the Twin Cities), with a
full re-calibration, **fixed spring** (MAM 0.77 → 1.01), kept summer balanced,
raised multi-decade `KGE_logKGE` 0.704 → 0.736, and made the melt factor
identifiable at a physical value. (It also over-cut *winter* ET and did not touch
the *fall* surplus — that's a separate groundwater-recession issue — but spring,
the target, was clearly an ET-phasing problem.)

The general takeaway is robust and not Crow-Wing-specific: **temperature-index ET
mis-phases wherever vegetation phenology is decoupled from temperature.**

## Suggested action (your call)

1. **A docs caveat** at the `evapotranspiration_method` decision point in
   `docs/source/configuration.rst` (right after the existing
   `ThornthwaiteChang2019` `.. note::`). Draft text:

   ```rst
   .. warning::

       **Temperature-index ET and vegetation phenology.** Because
       ``ThornthwaiteChang2019`` derives ET from temperature alone, it cannot
       represent seasonality when the actual ET cycle is *out of phase with
       temperature*. The clearest case is cold-region forests, where canopy
       leaf-out lags spring warming by weeks: Thornthwaite ramps ET up with the
       spring temperature rise, but little transpiration occurs before leaf-out,
       so early-spring ET is over-estimated. In snowmelt basins this can
       *consume the spring freshet* — evaporating meltwater that should appear
       as streamflow — and, because calibration then lowers ``et_scale`` to
       recover the annual water balance, it inflates flow in the remaining
       seasons. Where phenology and temperature are out of phase, prefer
       measured ET (``evapotranspiration_method: datafile``), ideally from a
       remote-sensing / NDVI-based product that follows the actual green-up.
   ```

   Could also be echoed in `model_description.rst`'s water-balance / ET section.

2. **(Optional) a phenology/NDVI-driven ET option** as future work — an ET factor
   that follows observed green-up rather than temperature. This is the *same*
   remote-sensing input as the land-cover DDF prior (issue #22), so the two are
   arguably one "forcing/priors from land cover" body of work.

3. **The seasonal mass-balance diagnostic — worth documenting / providing.**
   This is the test that actually cut through the confusion (we first blamed lake
   routing; it was ET). Recommended as a docs/methods section, and arguably a
   small reusable utility — it only needs the existing `hydrodata` columns plus a
   per-source discharge split. The method, for any calibrated run:

   For each season, tabulate the modeled basin-mean **P, ET, ΔStorage, and Q —
   with Q split by source** (fast soil/overland, slow groundwater, lake outlet) —
   against observed Q, plus monthly SWE/ET to locate melt timing. The patterns
   are diagnostic:

   - summer **ET > P** but Q still over → surplus is **slow-store release**
     (baseflow not receding), *not* an ET deficit;
   - freshet missing **and** spring ET ≈ P at melt → ET is **consuming the melt**
     (phasing — the case here);
   - a per-source flow that is **flat across seasons** → that reservoir is not
     responding seasonally (e.g. an over-buffered lake or a non-receding gw
     store).

   It converts "the seasonality is wrong" into "*which* flux / source / timing is
   wrong" — separating ET-phasing from routing, storage, and recession errors.
   Without it we'd have kept chasing the lake. Local implementation:
   `crow_wing_river/calib_frozen/diagnostics/diag_seasonal_mass_balance.py` (instruments
   `_advance_sub_catchment` / `_advance_lake` for the per-source split).

Full local write-up:
`~/dataanalysis/Wickert2026-hydroRaVENS-decadal-optimization/crow_wing_river/NOTES_seasonal_mismatch_ET_phenology.md`
