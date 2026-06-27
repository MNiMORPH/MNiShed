# Phenology refactor — `leafout_GDD` as the one knob (Crow Wing handoff)

*Follows `notes/HANDOFF_gdd_phenology_crow_wing.md` and its `_RESULTS`. From the
MNiShed side, 2026-06-27.*

## What changed

The phenology component is now collapsed to **one calibratable parameter** —
`leafout_GDD`, the thermal-time green-up threshold. The other six `phenology:`
values are fixed from priors (they are degenerate with `et_scale` or with the GDD
threshold, so calibrating them blows up identifiability). Three concrete changes
matter for your runs:

1. **`leafout_GDD` is now a calibration target** — usable as a `run_and_score`
   keyword, a declarative `target: leafout_GDD`, or a `ScoringModel` argument.
2. **`dormant_Kc` default is now 0.5** (was 0.4) — from your DJF over-production
   finding.
3. **The Thornthwaite ET is cached and only the GDD curve recomputes per eval**
   (~0.2 ms), so calibrating `leafout_GDD` does **not** slow the loop.

**Provenance / availability:** commits `80da2d3..25e7b01`, **pushed**
(`origin/master` @ `25e7b01`). The `mnished-jit` env is an editable install of
`~/models/MNiShed`, so it is already live — a *fresh* run picks it up; a
long-running process needs a restart (`git -C ~/models/MNiShed pull` if your
checkout is behind).

## 1. The one knob — fix it, or calibrate it

**Default (recommended): fix it** — zero free parameters. Set it in the config and
leave it out of the calibration. Your best run used the default `leafout_GDD: 100`
and it worked; if you have a regional leaf-out date (e.g. a USA-NPN spring-index
GDD for the basin), set it there.

```yaml
phenology:
  enabled: true
  leafout_GDD: 100        # fix from regional phenology if you have it
  dormant_Kc: 0.5         # now the default; set explicitly if you want another value
```

**Or calibrate it (one added parameter)** — add it as a target in `params.yml`:

```yaml
parameters:
  leafout_GDD:
    lower:   40           # ~ early-April leaf-out
    upper:   280          # ~ mid-June leaf-out
    initial: 100
    target:  leafout_GDD
```

…or directly: `run_and_score(cfg, ..., leafout_GDD=120)`.

**Do not calibrate the other phenology values** (`base_temperature__C`,
`full_canopy_GDD`, `dormant_Kc`, `full_Kc`, `senescence_*_doy`). They are fixed
from priors on purpose. Set them deliberately in the `phenology:` block — but the
senescence window is the one to *set with intent*, because…

## 2. Senescence is an active fall lever (your finding, now documented)

`senescence_start_doy` / `senescence_end_doy` control the autumn ET draw and
therefore fall discharge — they took your SON 1.21 → 1.01. The 260/305 default is
a temperate-deciduous starting point; adjust per basin if fall is over/under.

## 3. Closure: `none` + free `et_scale` is best (confirmed)

Your 0.740 config — `enforce_water_balance: 'none'` with free `et_scale` — is now
the documented recommendation: `et_scale` sets the annual level while
`leafout_GDD` sets the phase. A per-year/global multiplier can over-produce other
seasons when it can't tune the annual level. (Note: the closure trades against
melt-factor peakiness — sharper PDD → stronger freshet peaks but harder annual
balance.)

## Quick check it's live and calibratable

```python
from mnished import run_and_score
a = run_and_score(cfg, ..., enforce_water_balance='none', leafout_GDD=60).score
b = run_and_score(cfg, ..., enforce_water_balance='none', leafout_GDD=200).score
print(a, b)   # should differ — leafout_GDD now moves the fit
```

## Suggested experiment

You fixed `leafout_GDD = 100` and got 0.740. Worth one comparison: add it as a
target (bounds ~40–280) and see whether the freshet timing tightens further —
**but** watch for it trading against PDD / `et_scale`. If the improvement isn't
identifiable (KGE flat, parameter wandering), fix it from regional phenology and
keep the zero-parameter version. That call is itself the useful result.
