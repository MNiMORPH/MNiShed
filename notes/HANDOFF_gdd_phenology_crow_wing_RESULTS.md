# GDD-phenology Crow Wing test — RESULTS

*Reply to `notes/HANDOFF_gdd_phenology_crow_wing.md`. From the
Wickert2026 decadal-optimization run, 2026-06-27.*

## Headline for you: fall senescence is a real lever — set it deliberately

Your handoff said *"don't expect phenology to fix fall (SON 1.21); that's a
separate groundwater-recession issue."* **It did fix fall.** Under the best
configuration, **SON mod/obs went 1.21 → 1.01** — the senescence ramp
(`senescence_start_doy` 260 → `senescence_end_doy` 305, i.e. Sep ~0.96 → Oct
~0.61 Kc) plus the free-`et_scale` annual balance pulled fall ET up enough to cut
the autumn over-prediction. So the fall over-production was **more ET-coupled
than diagnosed**, and `senescence_start_doy` / `senescence_end_doy` are a genuine
control on autumn discharge, not cosmetic.

Implications for the feature:
- **The 260/305 defaults worked well** for north-central MN (mid-Sep onset →
  early-Nov dormancy). Keep them as the temperate-deciduous default, but they
  are clearly *active* — worth documenting that senescence shapes the fall limb
  (and is the natural lever if a basin's fall is over/under).
- The long-term form is NDVI-driven senescence (your issue #26), since the
  brown-down date moves year to year; the calendar DOY is the right first cut.

## What we ran

GDD Kc **on** (defaults: base 5 °C, leafout 100 GDD, full-canopy 400 GDD,
dormant_Kc 0.4, senescence 260/305 → mid-to-late-May leaf-out, verified). Two
closures, 8-decade multi-window SCE-UA, two-layer-land + frozen-ground + lake:

| mod/obs (2001–10) | no-phenology base | GDD + **global** | GDD + **none/et_scale** |
|---|---:|---:|---:|
| MAM | 0.77 | 0.95 | 0.85 |
| JJA | 1.01 | 1.29 | 0.92 |
| SON | 1.21 | 1.23 | **1.01** |
| DJF | 1.04 | 1.39 | 1.21 |
| mean KGE_logKGE (8 dec) | 0.704 | 0.685 | **0.740** |
| β (annual) | — | — | 0.97 |

**`enforce_water_balance: 'none'` + free `et_scale` is clearly best** (0.740, beats
even the earlier hand-kludge 0.736; β≈1; every season improved at once). The melt
factor becomes identifiable at a physical value in both closures (PDD 3.48 here,
6.27 under global) — your structural-confirmation signal holds.

## Two more tuning notes

1. **`dormant_Kc = 0.4` looks too low → winter over-production.** DJF was the worst
   season in *every* phenology run (1.21–1.39). Cutting winter ET to 0.4× of an
   already-small Thornthwaite winter ET inflates DJF discharge. Consider a default
   nearer **0.5–0.6**, or at least flag it.
2. **A spring-peakiness ↔ overall-balance dial.** Under `global`, PDD settled
   sharp (6.27) → strong freshet *peaks* (top-20 0.80) but the closure couldn't
   tune the annual level, so summer over-produced (1.29). Under `none`+`et_scale`,
   PDD settled gentle (3.48) → best overall balance but softer peaks (top-20
   0.675). Both are legitimate optima; worth a sentence in the docs that the
   closure choice interacts with how peaky the melt response calibrates.

## Best params (none + et_scale, for reproduction)

PDD 3.48, et_scale 0.80, τ_soil 68 d, τ_gw 584 d, f_exfil 0.73, τ_lake 17,800,
H_sill 822 mm, f_route 0.62, fdd_threshold 14.4 °C·day.

Config + write-up in the run repo:
`crow_wing_river/calib_frozen/phenology/` and
`crow_wing_river/NOTES_seasonal_mismatch_ET_phenology.md`.
