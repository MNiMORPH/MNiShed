# Photoperiod senescence + leaf-out prior (issue #35) — Crow Wing handoff

*From the MNiShed side, 2026-06-28. Builds on
`HANDOFF_gdd_phenology_crow_wing(_RESULTS).md` and
`HANDOFF_phenology_refactor_crow_wing.md`. Self-contained; copy/paste into the
Crow Wing run.*

Issue #35 added the **minimum-extra-data** half of the phenology story (the NDVI
complement is #26). Two pieces, both pushed-pending on the MNiShed master branch:

## 1. Photoperiod-driven autumn senescence (the new fall lever)

`senescence_method:` now switches the autumn brown-down between the old fixed
calendar window and a day-length cue:

```yaml
phenology:
  enabled: true
  leafout_GDD: 185                     # see (2) below
  senescence_method: photoperiod       # default is still 'doy' (back-compatible)
  senescence_photoperiod__hr: 11.0     # critical day length
  senescence_photoperiod_span__hr: 1.8 # decline span
```

It reads the existing `Photoperiod [hr]` forcing column — no new input. The ramp
is gated to the post-solstice half-year so spring's short days don't trigger it.

### The one thing to know before you re-run

**The default critical day length is 11.0 hr, deliberately BELOW the ~12-hr
autumn equinox.** A photoperiod cue only transfers across latitude in the correct
direction (north senesces earlier) when it is crossed *after* the equinox; a
near- or above-equinox threshold transfers weakly or backwards. My first instinct
(12.5 hr, picked to reproduce your 260/305 calendar window) was wrong for exactly
this reason — it sat just above the equinox and transferred the wrong way.

**Consequence for Crow Wing:** at ~47°N the 11.0-hr default browns down **~2–3
weeks later** than your calibrated `doy` 260/305 window (midpoint ~day 300 vs
282). That moves the autumn ET draw-down later, so **your SON fit will change** —
the senescence lever that took SON 1.21 → 1.01 is now placed differently. Re-check
it. If photoperiod senescence pulls fall too late, you have two honest options:

- raise `senescence_photoperiod__hr` toward ~11.5 hr (earlier onset, still below
  the equinox — milder but correct transfer), or
- keep `senescence_method: doy` with your tuned 260/305 (zero-transfer but already
  validated on this basin).

Either is defensible; which one is the *useful result* to report. Don't assume the
photoperiod default beats your calendar window here — test it.

## 2. `leafout_GDD` from a regional leaf-out date (the prior, not a free knob)

`mnished.leafout_GDD_from_date` converts a regional green-up **date** into the
`leafout_GDD` prior by accumulating *your* forcing's GDD to that date — grounding
it in the basin's thermal climate instead of a fabricated latitude→GDD curve. The
date is what spring-index climatologies give (USA-NPN Extended Spring Index /
SI-x; Schwartz et al. 2013), and it carries the latitude dependence implicitly.

```python
import pandas as pd
from mnished import leafout_GDD_from_date
df = pd.read_csv("crow_wing_forcing.csv", parse_dates=["Date"])
leafout_GDD_from_date(df, 5, 20)     # north-central MN ~ late-May leaf-out -> ~185 GDD
```

On your forcing a May-20 leaf-out gives **~185 GDD** — independently within ~10%
of the ~200 GDD your SCE-UA calibration settled on (a nice corroboration, neither
informing the other). Use it as the *fixed* `leafout_GDD` (zero free parameters),
or as the centre of its calibration bounds.

## Provenance / availability

Commits `3c5f4f9..a7fe5cf` on `~/models/MNiShed` master (mechanism, helper, tests,
the below-equinox default fix, docs + figure). **Not yet pushed.** `mnished-jit`
is an editable install, so a fresh process picks it up after a
`git -C ~/models/MNiShed pull`; restart any long-running calibration. Full suite
260 passed, ruff clean.

## Suggested experiment (one comparison, then report whichever wins)

Re-run your best 8-decade config three ways — `doy` 260/305 (your current),
`photoperiod` 11.0/1.8 (new default), `photoperiod` 11.5/1.8 (earlier onset) —
holding everything else fixed. Report the seasonal mod/obs (esp. SON) and
`KGE_logKGE`. The honest outcome is the finding regardless of direction.

## RESULTS — 3-way re-validation (2026-06-29)

Ran the comparison on the in-repo `examples/crow_wing/` (8-decade multi-window
SCE-UA, ~3100 evals each to convergence, all params free except senescence held
at each setting):

| senescence | KGE_logKGE | leafout_GDD | DJF | MAM | JJA | SON |
|---|---|---|---|---|---|---|
| `doy` 260/305      | 0.711 | 149 | 1.05 | 0.87 | 0.95 | 1.13 |
| `photoperiod` 11.0 | 0.702 | 129 | 1.14 | 0.86 | 0.95 | 1.08 |
| `photoperiod` 11.5 | 0.717 | 206 | 0.93 | 0.83 | 1.01 | 1.16 |

**Finding: a wash.** The three span 0.702–0.717 — within single-run SCE-UA noise
on this (flat) objective; the very different `leafout_GDD`/PDD landing at
near-identical scores is the documented equifinality (#37). Photoperiod
senescence neither helps nor hurts the Crow Wing fit, and the seasonal shapes are
similar.

**Why this is expected / what it does and doesn't show.** A *single-basin* test
is structurally blind to photoperiod's actual benefit — **latitude-transferability**
(one setting placing brown-down correctly across many basins without re-tuning).
At one latitude, `doy` and `photoperiod` are just two equally-fittable
parameterizations of the same autumn shape. So this confirms the photoperiod
option is **sound and reasonable on real data** (doesn't degrade the fit;
physical params), but the transfer claim needs a **multi-basin latitude-gradient**
test (fix senescence, calibrate the rest, compare `doy` vs `photoperiod` transfer
across basins) — not yet run. Reinforces that photoperiod-vs-`doy` is a
refinement, not the substance; the demonstrated upgrade is the data-grounded
`leafout_GDD_from_date` prior (~185 from data vs ~200 calibrated).
