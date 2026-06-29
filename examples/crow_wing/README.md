# Crow Wing River — lakes, sub-catchments, and phenology in one basin

A worked example that exercises three of MNiShed's optional process modules
**together** in a single real catchment, with a physical reason for each. Where
the Cannon example (`../cannon_inverse/`) is about calibration and
identifiability, this one is about the *process structure*: how a lake, a split
land surface, and a vegetation-phenology ET coefficient combine to reproduce a
northern-forest hydrograph across eight decade-long windows.

**Catchment:** Crow Wing River, north-central Minnesota (~2335 km²; sandy
outwash, mixed forest, many lakes and wetlands)
**Period:** 1905–2024 daily forcing; calibrated over eight decade windows
(1931–2020)
**Metric:** `KGE_logKGE` — equal-weight composite of KGE (peaks) and logKGE
(low-flow timing), averaged across the eight windows

## The three modules, and why each is here

| Module | What it does here |
|--------|-------------------|
| **Lakes** | A series outlet element (`kind: lake`): stores water above a sill and releases it as `Q = a·(H − H_sill)^(5/3)` (Manning wide-channel spillover), exchanging groundwater bidirectionally with its partner land zone. Buffers the part of the basin that drains through it. |
| **Sub-catchments** | Two parallel land zones — `direct_land` (drains straight to the gauge) and `lake_basin_land` (drains *through* the lake, fraction `f_route_lake`). Each is a two-reservoir cascade: a fast "soil" store over a slow "groundwater" store. |
| **Phenology** | A growing-degree-day vegetation coefficient (Kc) on ET. The one calibrated knob is `leafout_GDD`; the rest of the curve (dormant/full Kc, senescence DOYs) is fixed from regional phenology. |

**Frozen ground** (`frozen_ground: true`) supplies the snowmelt freshet: when
the frozen-ground index exceeds `fdd_threshold`, the soil's top reservoir is
forced to shed all snowmelt as direct runoff, then infiltration resumes on thaw.
For this sandy basin — which only sheds water when frozen — that is the dominant
overland-flow control. (The PDM saturation-excess store is disabled: it fired on
summer rain and was rejected by calibration.)

## Why phenology matters here — the freshet problem

Thornthwaite ET is temperature-only: it ramps up with spring warmth regardless
of whether the canopy has leafed out. In a northern mixed forest the canopy does
not transpire until **mid-to-late May**, so a temperature-only ET evaporates the
April snowmelt freshet that should be reaching the gauge — spring under-produces
and (via the low `et_scale` needed to close the annual balance) fall
over-produces.

The GDD Kc fixes this by holding ET near `dormant_Kc` until thermal-time
leaf-out (`leafout_GDD`), then ramping to `full_Kc` and back down through autumn
senescence. Calibrating `leafout_GDD` lets the data place green-up; for
north-central MN it settles near **~200 GDD = a late-May onset**, consistent
with regional phenology (common lilac flowers ~May 21 in Itasca County). That is
*more* physical than the model's generic ~100-GDD (mid-May, southern-MN) default
— a reminder that the GDD leaf-out prior is latitude-dependent.

## Running it

```bash
conda activate mnished-jit                  # numba JIT + spotpy
python calibrate.py sceua [reps]            # best-fit (SCE-UA), config-driven
python calibrate.py dream [reps] [iid|ar1]  # posterior / UQ (DREAM)
```

`calibrate.py` is the same generic, config-driven runner as the Cannon example —
no per-basin Python. It reads each parameter's `target:` in `params.yml` to know
where it maps in the model, builds the model once with `mnished.ScoringModel`,
and scores the eight-decade multi-window objective. Requires the `mnished-jit`
environment (`pip install mnished[jit] spotpy`; the Numba JIT needs
`numpy < 2.3`). See issue #20.

## Expected result

SCE-UA reaches a composite **`KGE_logKGE` ≈ 0.74** over the eight decades, with
`leafout_GDD` settling at a mid-to-late-May green-up and the melt factor at a
forest-physical value:

| parameter | value | note |
|-----------|------:|------|
| `leafout_GDD` | ~177 GDD | green-up ~mid-to-late May (physical for ~46–47°N) |
| `PDD_melt_factor` | ~2.5 | mm SWE/°C/day — forest-physical |
| `et_scale` | ~0.75 | free (no water-balance rescaling) |

Seasonal mod/obs for a representative decade (2001–2010):

| season | mod/obs |
|--------|--------:|
| DJF | ~1.08 |
| MAM | ~0.89 |
| JJA | ~0.95 |
| SON | ~1.10 |

It is the phenology `Kc` (suppressing early-spring ET until leaf-out) that
recovers the freshet; the `leafout_GDD` threshold itself is only weakly
constrained — calibrating it barely moves the composite score — so what the
calibration recovers is the *physical placement* of green-up, not a score gain.

## Honest open items

This example reproduces the *process behaviour* well; it is a research setup, not
a polished operational calibration. Two residuals remain visible:

- **A winter↔fall trade-off.** `KGE_logKGE` is fairly flat across the seasonal
  shape, so runs can swap a slightly high winter (DJF ~1.08) for a slightly high
  fall and back at nearly the same score. Calibrating with the
  `KGE_logKGE_seasonal` metric — which weights the four meteorological seasons
  equally rather than by volume — presses on this directly.
- **Groundwater recession.** The deep store's recession is gentle; part of the
  residual fall flow is a baseflow-shape issue, partly separate from ET phasing.

## Files

| File | Description |
|------|-------------|
| `crow_wing_config.yml` | Model config — lake + two land zones + frozen ground + GDD phenology |
| `crow_wing_forcing.csv` | Daily forcing and observed discharge, 1905–2024 |
| `params.yml` | Parameter bounds, `target:` mappings, module toggles, and the eight decade windows |
| `calibrate.py` | Generic config-driven in-process runner (SCE-UA / DREAM); identical to the Cannon example's |
