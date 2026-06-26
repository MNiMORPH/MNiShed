# Cannon River — parameter estimation (Dakota)

Calibrates MNiShed for the Cannon River catchment using efficient
global optimisation (EGO) followed by pattern search, via the
[Dakota](https://dakota.sandia.gov/) toolkit.

**Catchment:** Cannon River near Red Wing, MN (USGS 05355200; 3800 km²)
**Period:** 1992–1995 (daily)
**Metric:** KGE\_logKGE\_logFDC — equal-weight composite of KGE (peaks),
logKGE (low-flow timing), and KGE on log-FDC (flow-regime shape)

## Requirements

- `mnished` installed (`pip install mnished`)
- [Dakota](https://dakota.sandia.gov/) (tested with v6.x)
- `pyyaml` (`pip install pyyaml`)

If Dakota and mnished are in separate environments, activate the one
that has both before running.

## Workflow

**1. Configure parameters and modules**

Edit `params.yml` to set parameter bounds and enable/disable process
modules. The `modules` block controls both model behaviour and whether
the corresponding parameter is free or fixed in Dakota:

```yaml
modules:
  frozen_ground: true   # calibrates log__fdd_threshold
  direct_runoff: false  # fixes f_direct_runoff at 0 (bypass disabled)
```

**2. Generate the Dakota input file**

```bash
python generate_dakota_in.py
```

Re-run this whenever `params.yml` changes. The generated `dakota.in` is
overwritten and should not be edited by hand.

**3. Run calibration**

```bash
bash run.sh <short-description>
# e.g.: bash run.sh kge_3res_nogamma
```

Dakota runs EGO (global search) then pattern search (local refinement).
Results are archived to `runs/<timestamp>_<description>/`.

**4. Inspect results**

```bash
python plot_best.py --dat runs/<run>/evaluations.dat \
                    --save runs/<run>/best_fit.png
```

Prints logKGE, NSE, KGE, KGE\_logFDC, AIC, and BFI; saves a two-panel
diagnostic figure (hydrograph + flow duration curve).

## Bayesian calibration (DREAM)

The optimisation above returns a single best-fit point. To instead
characterise **how well the data constrains each parameter** — the full
posterior, including the parameter trade-offs (equifinality) a point estimate
hides — run Dakota's `bayes_calibration` with the DREAM sampler.

Where `driver.py` returns the scalar `1 − KGE`, the Bayesian driver returns the
**modelled log-flow vector** as Dakota `calibration_terms`; Dakota differences
it against the observed log-flows and forms a Gaussian likelihood, inferring
the error scale (`calibrate_error_multipliers`).

```bash
python make_bayes_data.py            # write the observed log-flow data (once)
python generate_dakota_in.py --bayes # generate dakota_bayes.in
dakota -i dakota_bayes.in -o dakota_bayes.out
```

The chain (`dakota_bayes.dat`) holds the posterior samples: per-parameter
marginal widths are the identifiability, and the posterior correlations are the
parameter trade-offs. For a fast, local complement that needs no MCMC, see
`mnished.identifiability` (parameter profiles + curvature eigenspectrum).

**Note:** DREAM needs many thousands of evaluations, so run with the Numba JIT
active (the pure-Python fallback is ~100× slower).

## Files

| File | Description |
|------|-------------|
| `cannon_cfg_template.yml` | Model config; physical parameters overridden by driver |
| `CannonTestInput.csv` | Daily forcing and observed discharge |
| `params.yml` | Parameter bounds, module toggles, and solver settings |
| `generate_dakota_in.py` | Generates `dakota.in` (or `dakota_bayes.in` with `--bayes`) |
| `dakota.in` | Dakota optimisation input (generated; do not edit by hand) |
| `driver.py` | Dakota evaluation driver (returns `1 − KGE`) |
| `run_driver.sh` | Shell wrapper Dakota calls per evaluation |
| `run.sh` | Runs Dakota, plots best fit, archives results |
| `archive_run.sh` | Copies outputs to `runs/<name>/` |
| `plot_best.py` | Re-runs and plots the best-fit parameter set |
| `make_bayes_data.py` | Writes the observed log-flow `calibration_data_file` |
| `driver_bayes.py` | Dakota DREAM driver (returns modelled log-flows) |
| `run_driver_bayes.sh` | Shell wrapper Dakota calls per Bayesian evaluation |
| `dakota_bayes.in` | Bayesian Dakota input (generated; do not edit by hand) |

## Forward run

For a single forward run with fixed parameters see `../cannon_forward/`.
