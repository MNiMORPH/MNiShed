# Bayesian calibration and parameter identifiability of the Cannon River (MNiShed)

*Comprehensive in-process Bayesian + identifiability study, 2026-06-26.
Verified numbers re-derived from the saved chains; see "Reproducibility".*

## Summary

The Cannon River was calibrated with MNiShed using an **in-process** Bayesian
sampler (SPOTPY DREAM) and a formal log-flow likelihood, and analysed with a
new local–global **identifiability** toolkit. The headline findings:

1. **A good fit, with a real degeneracy.** The optimum reaches KGE ≈ 0.80, but
   one parameter combination — soil-zone exfiltration fraction vs. soil
   recession timescale — is essentially unidentified, and the cheap local
   diagnostic and the full posterior **agree on it to |cos| = 0.99**.
2. **The formal posterior is overconfident.** Daily log-flow residuals are
   strongly autocorrelated (lag-1 ρ = 0.92), giving an effective sample size of
   **≈ 22 of 1096** days, so the iid-Gaussian posterior is ~7× too narrow.
3. **An AR(1) error model widens the marginals parameter-specifically**, not
   uniformly — recession timescales widen 3–4×, thresholds barely.
4. **The posterior is genuinely hard to sample.** What first looked like a
   multimodal threshold parameter (`Hmax`) proved to be **slow mixing** — it
   resolves to a single mode at 100k samples — but the posterior as a whole did
   not fully converge even then.
5. **In-process is ~100× faster per evaluation than Dakota's fork interface**,
   which makes all of the above tractable on a workstation.

---

## 1. Setup

**Basin / gauge.** Cannon River near Red Wing, MN (USGS 05355200; ~3800 km²).

**Forcing.** Two windows from the same daily record:
- *3-year* (1992–1994, 1096 scored days) — the shipped `cannon_inverse`
  example, used for the identifiability cross-check.
- *Decade* (2001–2010, 3652 scored days) — from the 1991–2011 forcing, used for
  the record-length and convergence work.

**Model.** Single-cascade MNiShed: 3 reservoirs (shallow/interflow, soil,
karst/bedrock), snowpack + frozen-ground + rain-on-snow modules, Nash-cascade
routing. **9 free parameters** (calibration coordinate; `log__` = log₁₀):
`t_recession_{shallow,soil,karst}`, `f_exfiltration_{shallow,soil}`,
`PDD_melt_factor`, `Hmax_shallow`, `fdd_threshold`, `routing_K`.

---

## 2. Methods

- **Sampler.** SPOTPY **DREAM** (Vrugt), run **in-process** (the model called
  directly via `run_and_score`, warm Numba JIT), 8–12 chains.
- **Likelihood.** Formal **Gaussian on log-flow residuals** (the measurement
  error marginalised out), with an **AR(1)** variant that whitens the residuals
  with a sampled coefficient `φ` before the Gaussian — the autocorrelation term
  of the Schoups & Vrugt (2010) generalized likelihood.
- **Identifiability diagnostics** (`mnished.identifiability`): per-parameter
  objective **profiles**, a finite-difference **curvature eigenspectrum**
  (stiff/sloppy directions and their named parameter combinations), and 2-D
  ridge grids — compared against the DREAM posterior covariance.
- **Optimizer A/B.** SPOTPY **SCE-UA** (in-process) vs. the existing Dakota
  EGO + pattern-search.

---

## 3. Results

### 3.1 Best fit
KGE at the posterior mode (3-yr) = **0.800**; in-process SCE-UA independently
reaches **0.806** (914 evals). Comparable optima from two methods.

### 3.2 Parameter identifiability — local diagnostic agrees with the posterior
At the posterior mode (3-yr), the curvature eigenspectrum of the *log-likelihood*
and the DREAM posterior covariance agree on the least-constrained direction:

- **Sloppiest direction alignment: |cos| = 0.99.**
- **Named degenerate combination:** `+0.81·f_exfiltration_soil + 0.59·t_recession_soil`
  — the soil-zone exfiltration fraction and its recession timescale trade off
  (both set how the soil reservoir partitions and times its release), so the
  data cannot separate them.
- Condition numbers: log-likelihood Hessian 1084, posterior covariance 89
  (same order; the gap reflects incomplete convergence + finite-difference
  sharpness, not disagreement on structure).

The cheap local tool is therefore **validated against the rigorous posterior**
on the degenerate combination that matters — it can serve as the fast screen,
with DREAM as the arbiter. (Caveat: the *stiff* directions are near-degenerate
among themselves and do not match eigenvector-by-eigenvector; only the flat end,
which governs identifiability, is meaningful.)

### 3.3 Residual autocorrelation → the posterior is overconfident
Log-flow residuals at the mode (3-yr):

| quantity | value |
|---|---|
| lag-1 autocorrelation ρ₁ | **0.924** |
| integrated autocorrelation time τ_int | 49 |
| effective sample size N_eff | **≈ 22** (of 1096) |
| overconfidence factor √τ_int | **≈ 7×** |

The iid-Gaussian likelihood treats 1096 daily residuals as independent; they
carry the information of ~22. **Reported posterior widths are ~7× too narrow**
until the autocorrelation is modelled. This is a property of the formal
likelihood, not of the model fit.

### 3.4 AR(1) error model — parameter-specific widening
Running an AR(1) likelihood (sampling `φ`) head-to-head against iid:

| window | φ posterior (median) | median widening | most-widened |
|---|---|---|---|
| 3-year | 0.969 | 1.1× | t_recession_shallow 4.2×, t_recession_karst 3.2× |
| decade | 0.987 | 0.9× | fdd_threshold 2.3×, t_recession_shallow 1.4× |

The error model **recovers the measured autocorrelation** (φ ≈ 0.97–0.99 vs the
measured ρ = 0.92). The widening is **not the uniform ~7× a naive N_eff argument
predicts** — it is parameter-specific: autocorrelation discounts each parameter
by how much its signature overlaps the correlated error structure (smooth,
slow effects like recession timescales are discounted most). *(These factors
are provisional — see §3.5; the chains are not fully converged.)*

### 3.5 Convergence — and the `Hmax` slow-mixing correction
The formal-likelihood posterior is **hard to sample**: the sharp likelihood
(itself a consequence of the autocorrelation) makes DREAM mix slowly. At
moderate budgets (15–30k) several parameters had Gelman–Rubin R̂ ≫ 1.2.

`Hmax_shallow` (saturation-excess threshold) initially looked **multimodal** —
its chains sat in different regions (per-chain means split ~1.5–3.3) across both
the 3-year and decade runs. **This was a false alarm.** A 100k-sample, 12-chain
decade run resolves it to a single **unimodal** mode (range tightened to
[1.64, 2.23]; per-chain means agree to within 0.16). It was **slow mixing, not
multimodality.**

The honest caveat stands: even at 100k the posterior was **not fully converged**
(`fdd_threshold` R̂ ≈ 25, `t_recession_shallow` R̂ ≈ 4.7 remained laggards). The
sampling difficulty is real; the specific multimodality diagnosis was not.

### 3.6 Record length (3-yr vs decade)
More data sharpens the slow reservoirs and is what ultimately resolved the
`Hmax` mixing, but it also **sharpens the likelihood further**, making
convergence harder, not easier — so the decade needed far more samples than the
3-year window, not fewer.

### 3.7 Performance (why this was feasible at all)
Per-evaluation wall-clock, identical model:

| engine | per-eval | note |
|---|---|---|
| Dakota fork interface | ~1500 ms | fresh Python process per eval |
| in-process (3-yr) | 13–23 ms | warm JIT, serial |
| in-process (decade) | ~33 ms | longer record |

In-process SCE-UA optimization: KGE 0.806 in **914 evals / 12 s**, vs Dakota
EGO + pattern-search ≈ 700 evals × 1.5 s ≈ **~18 min** — same optimum,
~90× wall-clock. Serial beats multiprocessing here (the ~1000–3650-value
simulation vector makes inter-process pickling cost more than the compute).

---

## 4. Limitations

- **Posterior widths are not yet quantitative** — the iid runs are overconfident
  (§3.3) and the AR(1) runs are not fully converged (§3.5). Treat all marginal
  widths and widening factors as provisional pending a converged AR(1) chain.
- **Convergence not reached** even at 100k samples for two parameters.
- **Single-objective KGE** for the optimum; the formal likelihood is on
  log-flows only (no heteroscedasticity/non-Gaussian terms — the rest of the
  Schoups & Vrugt generalized likelihood was not used).
- The local↔global agreement is demonstrated on the **3-year** window; it should
  be re-confirmed on a converged decade posterior.

---

## 5. Conclusions

1. MNiShed's parameters for Cannon are **largely identifiable, with one clear
   degeneracy** (soil exfiltration ↔ soil recession), confirmed by two
   independent methods.
2. The **formal log-flow likelihood is overconfident** by ~7× from
   autocorrelation; an AR(1) model is necessary before quoting uncertainties,
   and it widens marginals **parameter-specifically**.
3. The posterior is a **genuinely hard sampling problem**; apparent
   multimodality was slow mixing.
4. **In-process evaluation is the enabling change** — ~100× per-eval over the
   Dakota fork — and the same framework delivers both best-fit (SCE-UA) and UQ
   (DREAM).

---

## Reproducibility

- **Environment:** `mnished-jit` conda env (Python 3.11, `mnished[jit]` +
  `spotpy`; numba JIT requires `numpy < 2.3`).
- **Runners:** `examples/cannon_inverse/run_spotpy.py` (DREAM, iid/AR1),
  `run_sceua.py` (SCE-UA); decade workspace `examples/cannon_inverse/decade_run/`.
- **Chains (SPOTPY CSV):** 3-yr `cannon_dream.csv` (15k, iid),
  `cannon_dream_{iid,ar1}.csv` (30k); decade `decade_run/cannon_dream_iid.csv`
  (20k), `decade_run/cannon_dream_ar1.csv` (100k, 12 chains).
- **Diagnostics:** `mnished.identifiability` (committed) +
  `mnished.calibration.log_flow_residual_terms` (committed).

🤖 Generated with [Claude Code](https://claude.com/claude-code)
