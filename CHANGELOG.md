# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [3.0.0] - unreleased

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

### Added

- **Numba JIT-compiled time loop** with roughly two orders of magnitude
  speedup. Numba is an optional dependency; the pure-Python loop is the
  fallback and is used for the PDM and `et_water_stress` configurations.
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
  `channel_exit_water__volume_flow_rate` output, and flux-partition,
  frozen-ground-index, and actual-ET outputs.
- **Calibration API:** AIC free-parameter counting for the new junction,
  multipath, and recession-exponent parameters via `*_calibrated` arguments;
  single full-record ET multiplier in decade mode.
- **GitHub Actions CI** with NumPy 2.x in the test matrix.

### Fixed

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
