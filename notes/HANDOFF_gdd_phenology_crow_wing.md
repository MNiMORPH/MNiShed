# Testing the new GDD phenology Kc on Crow Wing

**What landed:** MNiShed now has an optional growing-degree-day vegetation
coefficient (Kc) that reshapes the Thornthwaite ET demand to follow thermal-time
leaf-out instead of temperature — suppressing early-spring ET so it stops
evaporating the snowmelt freshet. This is the in-package version of the prescribed
monthly leaf-out factor that fixed the spring residual on Crow Wing.

**Provenance:** MNiShed `master` @ `0088eb8` (version 3.1.0, unreleased). Commits
`c94e70a` impl / `4fa5ae7` tests / `0088eb8` docs. Origin finding:
`notes/HANDOFF_thornthwaite_ET_phenology.md`. Related: #25 (seasonal mass-balance
diagnostic), #26 (NDVI ET).

**No install step:** `mnished` is an editable install of `~/models/MNiShed` in the
`mnished-jit` env, so the feature is already live — `git -C ~/models/MNiShed pull`
if your checkout is behind.

## Enable it — add to your model config YAML

```yaml
phenology:
  enabled: true
  base_temperature__C:  5.0    # GDD base temperature (~deciduous)
  leafout_GDD:          100.0  # cumulative GDD (reset each Jan 1) at green-up onset
  full_canopy_GDD:      400.0  # cumulative GDD at full canopy
  dormant_Kc:           0.4    # ET coefficient outside the growing season
  full_Kc:              1.0    # ET coefficient at full canopy
  senescence_start_doy: 260    # fall ramp-down begins (~mid-Sep)
  senescence_end_doy:   305    # fully dormant (~early Nov)
```

Defaults give a mid-to-late-May leaf-out (north-central MN). Raise `leafout_GDD`
to delay green-up, lower to advance.

## Water-balance closure

`'global'`, `'water-year'`, and `et_scale`/stress closure all now normalise
against the phenology-adjusted demand, so the annual total is preserved exactly
whichever you use (the earlier "water-year isn't Kc-aware" gap is fixed — see the
`_demand_ET` water-balance fix, which also fixed a pre-existing Thornthwaite +
water-year non-closure unrelated to phenology). `'global'` is still a clean,
simple choice for this test.

**Note:** the fix landed in `~/models/MNiShed` after this test may have started —
a long-running Python process won't pick it up until restarted, so a fresh run
(`git -C ~/models/MNiShed pull` if needed) gets both the closure fix and phenology.

## What to expect / check

- **Recalibrate under global closure.** Expect spring MAM mod/obs to recover from
  ~0.77 toward ~1.0, summer to stay balanced, and — the strongest signal — the
  **melt factor to become identifiable at a physical value** (that's what told us
  the original fix was structural, not added flexibility).
- **Don't expect phenology to fix fall (SON 1.21).** That over-prediction was
  diagnosed as partly a *separate* groundwater-recession issue. If phenology *does*
  move fall, the calendar senescence window (`senescence_*_doy`) is the lever; if
  it doesn't, leave fall to the recession side.
- **Optionally sweep `leafout_GDD` / `full_canopy_GDD`** as a small sensitivity —
  but prefer fixing them from regional phenology/GDD literature (prior, not free
  parameter), consistent with the model's minimum-free-parameters ethos.

## Quick verification it's active

```python
import pandas as pd
from mnished import Buckets
b = Buckets(); b.initialize('your_config.yml', enforce_water_balance='global')
Kc = b.phenology_Kc()                       # daily array; scalar 1.0 if disabled
print(pd.Series(Kc, index=b.hydrodata['Date'].dt.month)
        .groupby(level=0).mean().round(2))  # monthly-mean Kc
# expect ~0.40 Jan-Apr, ~0.72 May, 1.00 Jun-Aug, declining Sep-Nov
```

## Known gaps (deliberate, for a first cut)

- Fall senescence is **day-of-year** (calendar); the spring/freshet limb is fully
  GDD-driven.
- Kc params are config-set, **not wired as calibration targets** (cache assumes
  fixed within a run) — calibrating them is a follow-up.

(The `'water-year'` closure gap noted in the first cut is now fixed.)
