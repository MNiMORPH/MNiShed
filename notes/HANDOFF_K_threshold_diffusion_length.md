# Handoff — K threshold and diffusion length as structural drivers in MNiShed

*Written 2026-06-24. From a calibration session in `~/dataanalysis/`.
Documents a physical insight from cross-basin calibration work that has
potential structural significance in how MNiShed sets up calibrations,
not just how it runs them. This is a thinking-piece for the MNiShed
session to decide whether and how to act on; it is not an implementation
spec.*

## The insight

The hydraulic diffusivity D = K·H/S and the calibration window t set a
diffusion length L_D = √(D·t) for hydraulic signals in a saturated
substrate. If L_D over the calibration window is much smaller than the
distance between a typical recharge location and the channel network
(or between adjacent tile drains), then **the matrix is hydraulically
disconnected from the gauge at the timescales the streamflow record can
see**. Matrix τ is unidentifiable from streamflow no matter what
architecture is used. Adding free DoFs for matrix τ in such cases just
gives the optimizer somewhere to park slack — the data does not
constrain those values.

## Threshold values worked out from cross-basin data

For a representative agricultural-basin geometry (H ~ 5 m, S ~ 0.05,
t = 10-year calibration window) the diffusion length comes out to:

| K (m/s) | L_D over 10 years | Comparable to … |
|---:|---:|---|
| 10⁻¹⁰ | 1.8 m | smaller than tile spacing — matrix silent |
| 10⁻⁸ | 18 m | tile spacing — tile network captures upper zone |
| 10⁻⁶ | 180 m | hillslope length — matrix enters baseflow signal |
| 10⁻⁴ | 4.3 km | basin-scale integration — decadal matrix-τ trends detectable |

The practical identifiability threshold is around **K ≈ 10⁻⁶ m/s**.
Above it, matrix τ is constrainable; below it, it is a nuisance
parameter that should be fixed at literature value or bounded tightly.

## Empirical confirmation from the Wickert2026 cross-basin work

- **Pure till basins** (Blue Earth, Le Sueur, Cottonwood) — bulk matrix
  K ~ 10⁻⁸ m/s; L_D over decades is below or near tile spacing.
  Multipath captures the engineered-drainage timescale cleanly; matrix τ
  in v6 lands at long timescales with broad uncertainty.
- **Wild Rice** (lake clay + till) — K ~ 10⁻⁹ to 10⁻⁷ m/s; L_D over the
  calibration window is below tile spacing for the matrix. The v8
  parallel sub-catchments architecture is physically defensible but the
  calibration sends deep-reservoir τ to upper bounds (decades-long
  timescales) because the data cannot pin them down. The till and clay
  zones come out hydrologically equivalent in the joint fit because
  matrix flow contributes a negligible fraction of the streamflow in
  either substrate.
- **Cannon River** — multi-layer with productive carbonate (PdC
  dolomite) and sandstone (Wonewoc) aquifers, K ~ 10⁻⁴ m/s. L_D over
  decades is multi-kilometer; the basin is hydraulically integrated. A
  matrix-τ shortening trend of −0.082/decade in log(τ_soil) over
  1931–2020 is detectable from streamflow alone. The signal exists
  because the productive aquifers connect the basin within the
  calibration window.

So the cross-basin pattern of where matrix τ identifiability shows up
versus where it doesn't tracks substrate K exactly as the diffusion-length
heuristic predicts. This is not an artifact — it is structural.

## Why this matters for MNiShed (and how it might matter structurally)

MNiShed currently lets users calibrate matrix τ freely on every
reservoir. There is no built-in awareness that the streamflow record may
not contain information about matrix τ for low-K substrates. The
optimizer happily produces decades-long matrix τ values that hit upper
bounds without the user noticing the parameter is just absorbing slack.

Possible structural responses to this, in increasing order of intrusiveness:

### 1. Documentation only (low risk, high value)

Add a section to `model_description.rst` (or a new `architecture_guide.rst`)
that walks through the K threshold and diffusion-length heuristic
*before* a user picks reservoirs and writes a calibration YAML.
Pattern: "before you decide whether to calibrate matrix τ, compute
L_D over your calibration window and compare to tile spacing and
hillslope length." Reference the [Brutsaert-Nieber recession analysis]
guidance as the complementary check on bulk Q(H) shape.

### 2. A diagnostic class

Add a small class, analogous to `BrutsaertNieber`, that takes
K (or a substrate name), H, S, basin geometry, calibration window, and
reports L_D and a verdict ("matrix τ likely identifiable / borderline /
not identifiable"). Use case: a user types
`MatrixIdentifiability(K=1e-8, H=5, S=0.05, calib_yr=70).verdict()` and
gets a clear answer before running calibrations. Output could feed into
`Priors`.

### 3. Integration into Priors

Extend the existing Priors module so that, given basin geometry and
substrate K (or a literature-default substrate selector), it returns
suggested matrix τ priors (or, more usefully, suggests whether matrix
τ should be calibrated at all). This is a natural extension of the
existing pattern — Priors already suggests parameter ranges from
observed data; this adds physics-based bounds from substrate
characteristics.

### 4. A warning in `run_and_score`

When `recession_coeff` is calibrated for a reservoir whose substrate K
implies non-identifiability of matrix τ, raise a `UserWarning` at
calibration setup. Similar to the existing
`enforce_water_balance='none'` warning. Optional opt-in via a
`substrate_K=` kwarg or a `--check_identifiability` flag in the CLI.

### 5. A "lock matrix τ to literature" mode

Allow the user to declare `recession_coeff: literature` for a
reservoir, in which case MNiShed picks a value from a built-in lookup
table keyed on substrate type (clay-rich till, lake clay, dolomite,
sandstone, sand-gravel, etc.) and basin size. Removes the parameter
from the calibration entirely. This is the most invasive option but
also the most useful for users who do not yet have the physical
intuition to set bounds by hand.

## Connections to existing MNiShed features

- **Multipath drainage** (already in v3.0.0 scope) captures the
  engineered-drainage timescale that is the *actual* identifiable
  parameter in low-K basins. The K threshold tells you when multipath
  is essentially the only parameter you have meaningful information
  about; matrix τ is along for the ride.
- **BrutsaertNieber MRC analysis** gives you the bulk Q(H) shape from
  observed data; the K threshold tells you which reservoirs in the
  resulting architecture are constrainable. The two diagnostics are
  complementary: B-N is data-side, K is physics-side; together they
  bracket what the calibration can find.
- **Priors / suggest_priors** is the natural home for K-based bounds
  on matrix τ. Priors currently works from observed-data diagnostics;
  this would add a substrate-physics path.

## Open design questions for the MNiShed session

These need product-design judgment more than coding:

1. **What's the right level of intrusiveness?** Documentation, diagnostic
   class, opt-in warning, or default warning? Probably some combination,
   but which.
2. **Where does the substrate K come from?** User-provided, looked up
   from a built-in table keyed on substrate name, or both?
3. **Should this be a separate PR or rolled into the next big feature
   PR?** I'd argue separate — it's conceptually distinct from
   sub-catchments and from multipath. Could land as v3.2.0 after
   sub-catchments (v3.1.0).
4. **How to handle anisotropic K?** Vertical vs horizontal differ by
   factors of 10–100 in many substrates. The diffusion length using
   bulk K is a simplification.
5. **Sensitivity to S and H assumptions.** L_D depends on both; the
   verdict could flip with reasonable variation. Worth showing a
   range, not a single number.

## Suggested triage

If you take only one thing forward, take **(1) Documentation only** as
the first PR — it's stewardship-friendly, low-risk, and immediately
useful. Then **(2) Diagnostic class** as a follow-up if the user
community asks for it. Skip (4) and (5) unless there's clear demand —
they trade safety for prescriptiveness and that's a judgment call best
made with user feedback.

## Where to find the full reasoning

`/home/awickert/.claude/projects/-home-awickert-dataanalysis/memory/feedback_k_threshold_diffusion_length.md`

That memory entry contains the full derivation, the cross-basin
empirical evidence, the comparison table, and the suggested workflow
for applying the K-threshold check on a new basin during a calibration
session. The MNiShed-side memory (`project_k_threshold_handoff.md`)
points back to this handoff and that calibration-side memory entry.

## Constraints (per `~/.claude/CLAUDE.md`)

Same as the existing MNiShed handoffs: stewardship-mode work, explicit
per-step authorization required for commits, tags, pushes, releases.
This handoff is descriptive, not prescriptive — implement what makes
sense for MNiShed's user community, in whatever form makes sense, when
the time is right. After the sub-catchments PR and probably after
v3.1.0.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
