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

## In-process calibration (SPOTPY — no Dakota)

The Dakota workflows above fork a fresh Python process per evaluation (~1.5 s of
interpreter startup each); for a fast model on a workstation that overhead
dominates. **In-process SPOTPY** calls the model directly with the Numba JIT
warm — ~100× faster per evaluation — and offers both an optimiser (SCE-UA) and a
Bayesian sampler (DREAM) on the same setup.

`calibrate.py` is the generic, **config-driven** runner — it works for any basin
with no per-basin Python, because each parameter's `target` in `params.yml`
declares where it maps in the model (see the `parameters:` comment there).
`run_sceua.py` / `run_spotpy.py` are the explicit, hand-mapped equivalents.

```bash
conda activate mnished-jit                 # numba JIT + spotpy
python calibrate.py sceua [reps]           # best-fit (SCE-UA), config-driven
python calibrate.py dream [reps] [iid|ar1] # posterior / UQ (DREAM), config-driven
# explicit per-basin equivalents:
python run_sceua.py  [reps]
python run_spotpy.py [reps] [iid|ar1]
```

- `run_sceua.py` builds the model once with `mnished.ScoringModel` and reuses it
  every evaluation (no per-eval CSV re-read or reconstruction), reaching the same
  optimum as Dakota EGO + pattern-search in a small fraction of the wall-clock.
- `run_spotpy.py` runs DREAM with a formal log-flow likelihood (`iid`, or `ar1`
  to account for the strong day-to-day autocorrelation of streamflow residuals);
  the saved chain is the posterior.
- **Serial is fastest here** — each evaluation returns a long simulation vector,
  so multiprocessing loses more to inter-process pickling than it saves.

Requires the `mnished-jit` environment (`pip install mnished[jit] spotpy`; the
Numba JIT needs `numpy < 2.3`). See issue #20.

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
| `calibrate.py` | Generic config-driven in-process runner (SCE-UA / DREAM; reads `target:` mappings) |
| `run_sceua.py` | In-process SCE-UA optimiser (SPOTPY; builds once via `ScoringModel`) |
| `run_spotpy.py` | In-process DREAM sampler (SPOTPY; iid / AR(1) log-flow likelihood) |

## Forward run

For a single forward run with fixed parameters see `../cannon_forward/`.
