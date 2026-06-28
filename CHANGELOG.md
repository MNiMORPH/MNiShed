# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [3.2.0] - 2026-06-28

### Added

- Lake (open-water) sub-catchments: a `kind: lake` sub-catchment is a single
  open-water store with a threshold power-law stage–discharge outlet
  (`Q_out = a*(H - H_sill)^b`, `b = 5/3` by default), fed by direct
  precipitation minus open-water evaporation (the basin `et_scale`; no separate
  Penman lake-ET term). It is coupled to a land sub-catchment's deepest
  reservoir by a bidirectional, volume-conserving groundwater exchange `Q_gw`
  that reuses that reservoir's recession coefficient and exponent (no new
  calibrated parameter) and flips direction with the head difference, giving
  seasonal store-and-release buffering. Configure with a `lake:` block
  (`outflow_coefficient`, `sill_storage__mm`, `outflow_exponent`, optional
  `gw_partner`, `f_route_lake`); calibrate the outlet through the
  `sub_catchments` override of `run_and_score` (lake `recession_coeff = 1/a`,
  `H_threshold = H_sill`). Supported on both the pure-Python and Numba JIT time
  loops (verified identical). New `Buckets.has_lake`; basins without a lake are
  unchanged.
- Channelized routing through a lake (`f_route_lake` in `[0, 1]`): a fraction of
  the partner land zone's discharge is routed through the lake — buffered by its
  outlet and re-released — instead of reaching the gauge directly, so the land
  cascade no longer has to fake slow release by stretching its recession
  timescale. The transfer is instantaneous (no channel lag, buffered by lake
  storage); `update()` runs two passes (land, then lakes) so a lake sees its
  partner's current-step discharge. `f_route_lake` is data-derived (lake network
  position), never calibrated; `0` (default) keeps the lake disconnected
  (terminal/closed lake). Supported on the pure-Python and JIT loops (verified
  identical). Automatically deriving `f_route_lake` and lake network position
  from terrain remains planned with the drainage-density / hydraulic-conductivity
  work (MNiMORPH/MNiShed#19).

- In-process, config-driven calibration: `Calibrator` (with `target_kwargs`)
  defines the calibration problem declaratively — each parameter's `target:` in a
  `params.yml` maps it into the model (flat `name` / `name[i]`, or nested
  `sub_catchments[I].key[j]` for multi-zone and lake setups), with no per-basin
  Python. It is sampler-agnostic: point any optimizer or sampler at
  `Calibrator.score`, which is bit-identical to the equivalent `run_and_score`
  call.
- `ScoringModel`: a build-once / score-many wrapper that constructs the model and
  reads the forcing once, then reuses a fresh copy per evaluation. For thousands
  of evaluations this removes the per-evaluation rebuild that otherwise dominates
  the wall-clock; `Calibrator` uses it internally. Bit-identical to `run_and_score`.
- Multi-window calibration objectives: `Calibrator.score_windows` scores a
  parameter set on each of several calibration windows (a driver `decades:` list
  of `{start, end}`), returning one `CalibResult` per window and leaving
  aggregation — a mean score, or concatenated residuals for a likelihood — to the
  caller (the decadal-backbone objective). Single-window configs are unchanged.
- `log_flow_residual_terms`: exposes the scored log-flow residual terms (observed
  and modeled log-discharge over the scoring mask) used by a formal Gaussian
  likelihood, so an external Bayesian sampler can build its own objective on
  MNiShed's scored flows.
- Parameter-identifiability diagnostics (`mnished.identifiability`): post-fit tools
  that characterise how well the data constrain each parameter — per-parameter
  `profile` / `profile_all` likelihood profiles, a curvature `eigenspectrum` that
  names stiff vs. degenerate (sloppy) parameter combinations, 2-D `ridge` maps, and
  an `IdentifiabilityReport`, with plotting. Bound-pegged parameters are excluded
  from the curvature spectrum so a single pegged parameter no longer voids it.
- The Cannon example gains in-process SPOTPY calibration runners — SCE-UA
  (best-fit) and DREAM (posterior / UQ, with iid and AR(1) log-flow likelihoods)
  driven by `Calibrator` — plus a Dakota DREAM Bayesian path. The in-process,
  JIT-warm evaluation is roughly 100× faster per evaluation than the Dakota fork
  interface, which remains available for expensive or cluster runs.

- Optional growing-degree-day vegetation-phenology coefficient (`Kc`) on the
  `ThornthwaiteChang2019` ET demand, configured with a `phenology:` block (off by
  default). `Kc` ramps from `dormant_Kc` to `full_Kc` as accumulated GDD (base
  `base_temperature__C`, reset each calendar year) rise from `leafout_GDD` to
  `full_canopy_GDD`, holds through summer, then declines over a day-of-year
  senescence window. Because GDD enters nonlinearly it corrects seasonal *phasing*
  rather than being absorbed by the annual `et_scale`, so it stops temperature-
  index ET from evaporating the snowmelt freshet before leaf-out. Applied to the
  demand so the water-balance correction preserves the annual total; supported on
  the pure-Python and JIT loops (verified identical). Of the curve's parameters,
  only `leafout_GDD` (green-up timing) is exposed as a calibration target
  (`run_and_score(leafout_GDD=...)` / `target: leafout_GDD`) — the rest are fixed
  from priors to avoid degeneracy (`full_Kc`/magnitude with `et_scale`,
  `base_temperature__C` with the threshold). Validated on the Crow Wing River
  (best multi-decade KGE_logKGE 0.740, every season improved); the senescence
  window is an active autumn-discharge lever, and `dormant_Kc` defaults to 0.5.
  `Kc` applies to land ET only — lake open-water evaporation is phenology-free
  (open water has no leaf phenology).

- `SeasonalMassBalance` diagnostic (`mnished.diagnostics`): a per-season
  water-balance decomposition (P, ET, storage change, discharge) with discharge
  **split by source** — fast/event, slow/baseflow, and lake outlet — against
  observed discharge, plus a monthly snowpack/ET climatology for melt timing. It
  separates an ET-phasing error from a routing, storage, or recession error
  (the decomposition that, on the Crow Wing River, reframed a suspected
  lake-routing problem as a Thornthwaite spring-ET phasing issue). Backed by a new
  `run(store_fluxes=True)` / `run_and_score(..., store_fluxes=True)` recording
  mode whose three source columns sum exactly to the modeled discharge
  (MNiMORPH/MNiShed#25).

- Worked example `examples/crow_wing/` (Crow Wing River, north-central Minnesota):
  a *process-structure* example exercising lakes, parallel sub-catchments, and
  growing-degree-day phenology together in one real basin (plus frozen ground),
  driven by the generic config-driven runner, with a documentation page covering
  the freshet/phenology narrative, the latitude-dependent `leafout_GDD` finding,
  and the seasonal mass-balance diagnostic (MNiMORPH/MNiShed#27, #28).

### Fixed

- Water-year ET scaling now normalises against the actual ET demand that is
  applied (`Buckets._demand_ET()`), not the raw input `Evapotranspiration` column.
  Previously `enforce_water_balance: 'water-year'` with
  `evapotranspiration_method: ThornthwaiteChang2019` divided `P - Q` by the input
  ET column while applying the multiplier to the Thornthwaite demand, so the
  per-water-year balance did not close (off by `mean(Thornthwaite)/mean(column)`,
  ~2× on test data). `datafile` mode and `'global'` closure are unaffected
  (bit-identical); only Thornthwaite + water-year results change, now closing
  correctly. Phenology is water-balance-aware under every closure mode.

## [3.1.0] - 2026-06-24

### Added

- Parallel sub-catchments: a basin can be partitioned into spatially distinct
  zones that drain to the same channel in parallel, each with its own
  reservoir cascade and snowpack/frozen-ground state. Basin discharge and
  storage are the area-weighted means over sub-catchments. Configure with a
  `sub_catchments:` YAML block (each entry has a `name`, `area_fraction`, its
  own `reservoirs` block, and optional `initial_conditions`); calibrate by
  passing a `sub_catchments=[...]` argument to `run_and_score`. Supported on
  both the pure-Python and Numba JIT time loops. A single sub-catchment of
  area 1.0 reproduces the previous single-cascade behaviour exactly, so
  existing configurations and calls are unchanged. New `SubCatchment` class
  and `Buckets.sub_catchments` / `Buckets.n_sub_catchments`.
  `run_and_score` chains per-sub-catchment storage state (reservoir depths,
  snowpack, frozen-ground index, carried deficit) across decade windows:
  `final_states` / `initial_states` are nested per sub-catchment when there
  are several, and stay flat/scalar for a single sub-catchment.
- The Numba JIT now covers **PDM saturation-excess** (`pdm_H0`) and
  **`et_water_stress`**; those configurations previously fell back to the
  pure-Python loop. The JIT and pure-Python loops remain verified-identical, so
  the JIT is now used for every supported configuration whenever Numba is
  importable.

### Changed

- The pure-Python time-loop fallback is no longer silent. `Buckets.run()` emits
  a one-time `UserWarning` when Numba is installed but fails to import (usually
  a NumPy/Numba version mismatch), so an unexpected ~100× slowdown is visible.
  A plain "Numba not installed" stays quiet, since pure Python is the expected
  default without the `jit` extra.

### Deprecated

- The flat single-sub-catchment state shape (`{'reservoirs': [...],
  'snowpack': ..., 'fgi': ...}`) for `run_and_score`'s `initial_states` /
  `post_spinup_states`. Passing it now emits a `DeprecationWarning`; it will be
  removed in v4.0 in favour of the uniform per-sub-catchment form
  (`{'sub_catchments': [...]}`). See
  [#18](https://github.com/MNiMORPH/MNiShed/issues/18). The nested form is
  accepted at any number of sub-catchments (including one) and does not warn.

### Fixed

- `run_and_score` now validates chained `initial_states` / `post_spinup_states`
  and raises a clear `ValueError` if they contain non-finite (NaN/inf) values,
  naming the offending key/index. Previously a NaN state from a partial-data or
  failed decade propagated silently — every modelled flow became NaN and the
  score looked merely poor rather than broken. (`None` reservoir entries in
  `post_spinup_states`, meaning "keep the spin-up value", are still allowed.)
- The BMI `snowpack__liquid_equivalent_depth` and
  `land_surface__frozen_ground_index` outputs are now the area-weighted basin
  mean over sub-catchments, instead of reporting only the first sub-catchment.
  Exact for a single sub-catchment, so K=1 couplers are unchanged
  ([#16](https://github.com/MNiMORPH/MNiShed/issues/16)).

## [3.0.0] - 2026-06-23

Major release. The library was renamed from **hydroRaVENS** to **MNiShed**,
the numerical core moved to a Numba JIT time loop, and many new reservoir
mechanics, hydrologic processes, analysis tools, and a CSDMS BMI wrapper
were added. Several renames are breaking; there are no backwards-compat
shims (clean break). See `RELEASE_NOTES_v3.0.0_DRAFT.md` for the full notes
and migration guide.

### Changed (breaking)

- Renamed the package `hydroRaVENS` → `mnished`
  (`pip install mnished`, `from mnished import ...`).
- Renamed the BMI class `BmiHydroRaVENS` → `BmiMNiShed`.
- Renamed the `Reservoir` constructor argument `t_efold` / `t_recession`
  → `recession_coeff`.
- Renamed the YAML key `e_folding_residence_times__days`
  → `recession_coefficients`.
- Renamed the attribute `Reservoir.H_infiltrated` → `Reservoir.H_to_next`.
- Replaced the boolean `scale_et` with the string-valued
  `enforce_water_balance` (`'water-year'` / `'global'` / `'none'`), as both
  a `run_and_score()` argument and a YAML key.
- Corrected the misspelled `evapotranspiration_method` value
  `ThorntwaiteChang2019` → `ThornthwaiteChang2019` (config files must update
  the value).
- Updated the BMI input/output variable names to current CSDMS Standard
  Names conventions: air temperature → `atmosphere_bottom_air__temperature`;
  daily extremes → `…__time_min_of_temperature` / `…__time_max_of_temperature`;
  discharge → `channel_exit_water_x-section__volume_flow_rate`; ET forcing
  input → `land_surface_water__uncorrected_evapotranspiration_volume_flux`.
  Existing `BmiMNiShed` couplers must update these variable names.

### Added

- **Numba JIT-compiled time loop** with roughly two orders of magnitude
  speedup. Numba is an optional dependency, installable via the `jit` extra
  (`pip install mnished[jit]`); the pure-Python loop is the fallback and is
  used for the PDM and `et_water_stress` configurations.
- **Power-law (nonlinear) recession** (`recession_exponents`) with exact
  integration, plus `Reservoir.mean_residence_time()`.
- **Junction types between reservoirs**: `fraction` (default), `leakance`
  (head-difference flow through a confining unit, `leakance_R__days`), and
  `threshold` (dead-storage cutoff, `H_threshold__mm`).
- **Tile-drain sub-reservoir** (`tile_fractions` / `tile_residence_times__days`).
- **Multipath threshold-activated parallel drain**
  (`multipath_thresholds__mm` / `multipath_timescales__days`).
- **Exponential PDM saturation-excess overland flow** (`pdm_H0`).
- **Hydrologic processes:** snowpack insulation (`snow_insulation_k`),
  DTR-based FGI decay, constant regional baseflow (`baseflow_Q`),
  direct-runoff bypass (`direct_runoff_fraction`), and three ET module
  toggles (`et_reservoir_draw`, `et_water_stress`, `et_scale`).
- **Run controls:** top-level `modules:` block, automatic spin-up cycle
  count, analytical steady-state initial conditions, and `post_spinup_states`
  for decade-chained initial storage.
- **Analysis and diagnostics:** `BrutsaertNieber` recession analysis,
  `HydrographSeparation`, `Priors` / `suggest_priors`, composite KGE metrics
  (`KGE_logKGE`, `KGE_logKGE_logFDC`, `KGE_logKGE_logFDC_BFI`,
  `logKGE_logFDC_BFI`), and `kge_logfdc` on `CalibResult`.
- **CSDMS BMI wrapper** (`BmiMNiShed`): 10-reservoir cap (up from 3),
  `channel_exit_water_x-section__volume_flow_rate` output, and
  flux-partition, frozen-ground-index, and actual-ET outputs. Input/output
  names follow current CSDMS Standard Names conventions
  (`atmosphere_bottom_air__temperature`, `…__time_min_of_temperature` /
  `…__time_max_of_temperature`); the ET forcing input is
  `land_surface_water__uncorrected_evapotranspiration_volume_flux` (a
  model-specific name for the as-supplied ET, before water-balance
  correction), distinct from the corrected ET output.
- **Calibration API:** AIC free-parameter counting for the new junction,
  multipath, and recession-exponent parameters via `*_calibrated` arguments;
  single full-record ET multiplier in decade mode.
- **GitHub Actions CI** with NumPy 2.x in the test matrix.
- **Performance benchmarks** (`benchmarks/bench_jit.py`, `plot_jit.py`) and a
  **Performance** documentation page quantifying the Numba JIT speedup
  (≈100×–400×, scaling with record length; ~72 µs vs. ~0.2 µs per simulated
  day).

### Fixed

- `Buckets.initialize()` now **raises** `FileNotFoundError` /
  `yaml.YAMLError` on a missing or unparseable config file instead of
  printing and calling `sys.exit(2)`, so library, BMI, and notebook callers
  can handle the error; the `mnished` CLI catches it and exits cleanly.
- **Numba JIT correctness:** the JIT time loop dropped negative-ET
  (condensation) input to the soil reservoir, and omitted tile storage from
  the subsurface total on skipped (missing-forcing) days. Both now match the
  pure-Python loop. (Affected only runs where the JIT actually executed.)
- **`recession_H_ref` standardized to 1.0** and
  `Reservoir.mean_residence_time()` corrected — for nonlinear reservoirs it
  was off by the `H_ref^((b-1)/b)` gauge factor; `run_and_score` no longer
  sets a hidden `H_ref` anchor.
- ET-reservoir-draw condensation (negative ET) is now capped at `Hmax`, with
  the surplus shed to runoff instead of stored above the cap.
- AIC free-parameter counting counts only finite `Hmax` entries
  (`inf` = no cap, not a calibrated parameter).
- The `land_surface__frozen_ground_index` BMI output returns `NaN` before the
  first `update()`, consistent with the other outputs.
- `H_deficit_carry` is now reset before full-record spin-up, before
  pre-decade spin-up, and when initial / post-spin-up states are applied
  (affected multi-decade chained calibrations).
- `f_to_discharge` index alignment for sparse reservoir lists.
- NaN propagation through the ET multiplier, Nash-cascade routing, and
  Thornthwaite ET.
- `compute_water_year()` off-by-one for January start months.
- `compute_ET()` idempotency on repeat invocations.
- Configuration robustness and YAML key/unit fixes
  (`__m` → `__mm` in several keys); documentation API-example fixes.

### Removed

- `Snowpack.discharge()` and the dead `H_discharge` attribute; the obsolete
  `scalar_dt` (the daily timestep is documented as a design choice).
- No backwards-compatibility shims are provided for the renames above.

[Unreleased]: https://github.com/MNiMORPH/MNiShed/compare/v3.2.0...HEAD
[3.2.0]: https://github.com/MNiMORPH/MNiShed/compare/v3.1.0...v3.2.0
[3.1.0]: https://github.com/MNiMORPH/MNiShed/compare/v3.0.0...v3.1.0
[3.0.0]: https://github.com/MNiMORPH/MNiShed/compare/v2.3.0...v3.0.0
