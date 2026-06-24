# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[3.0.0]: https://github.com/MNiMORPH/MNiShed/compare/v2.3.0...v3.0.0
