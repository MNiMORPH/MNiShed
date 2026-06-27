# Handoff — Enable `f_route_lake > 0` (channelized inflow through the lake)

> **STATUS — IMPLEMENTED on `master` (2026-06-25).** `f_route_lake > 0` now
> ships: schema/validation + two-pass `update()` + routed inflow + JIT mirror
> (commit `b63720a`), with docs (`367c062`). The `NotImplementedError` is gone;
> `f_route_lake` is accepted in `[0, 1]`, requires a `gw_partner` (the routing
> source) when `> 0`, and is data-derived (never calibrated). Verified: basin
> mass balance closes for `f_route_lake ∈ {0, 0.5, 1}`; `f_route_lake = 0` is a
> bit-identical no-op; JIT == pure-Python for `f_route_lake ∈ {0, 0.5, 0.9, 1}`.
> The empirical motivation and the calibration-side **validation plan below
> still stand** — they are the next calibration step on Crow Wing.

*Written 2026-06-25. From a calibration session in `~/dataanalysis/`. Documents
the empirical motivation for shipping the `f_route_lake` half of lake support
(now implemented; was a `NotImplementedError`) and what the calibration side
needed from MNiShed before lake-rich basins could be honestly fit. Companion to
`HANDOFF_lake_support.md` and `DESIGN_lakes.md`. This handoff is descriptive and
motivation-focused, not a code spec — the design and triage decisions are laid
out in `DESIGN_lakes.md`.*

## TL;DR

Today's Crow Wing tests confirmed the prediction from `HANDOFF_lake_support.md`
and `DESIGN_lakes.md`: the lake-support v1 (Q_gw only, `f_route_lake = 0`) is
necessary but not sufficient. With the lake's groundwater capacitor doing real
physical work, **the missing piece is channelized routing of upstream land
discharge through the lake** (issue #19 / `f_route_lake` enablement). Until that
ships, lake-rich basins cannot be honestly fit — the model either pegs substrate
timescales (no lake), or uses an unrelated fast-pathway mechanism as a proxy
(multipath was never meant for this).

## Empirical motivation — Crow Wing v1 / v2 results

Crow Wing River at Nimrod (USGS 05244000, ~2 335 km², ~30 % lake area). Single
decade 2001–2010, otherwise identical configs except for the multipath mechanism:

| Test | Land cascade | KGE | τ_matrix (land) | Lake outflow params | Reading |
|------|--------------|----:|----------------:|---------------------|---------|
| v6 etscale baseline (no lake) | 1-res + multipath | 0.597 | 956 d ⚠ | — | substrate stretched to lake-like timescales |
| v1 lake + multipath | 1-res + multipath, lake sub-catchment | 0.703 | 321 d | `a` and `H_sill` at bounds | τ relaxed — but partly via multipath, which is the wrong physics for this basin |
| **v2 lake, NO multipath** | 1-res only, lake sub-catchment | **0.626** | **974 d** ⚠ | `a` = 0.0014, `H_sill` = 275 mm (both interior) | τ sprang back to baseline — lake alone can't carry the storage-and-routing role |

The v2 result is the honest one. Removing multipath isolates the lake's actual
contribution: `Q_gw` is doing real work (the lake's own parameters are interior
and physically sensible), but it cannot relax the till substrate's recession
timescale on its own. **The land cascade still needs an outlet through which to
release its storm-generated discharge with the right buffering**, and that outlet
is the lake.

Re-enabling multipath as a calibration crutch would be the wrong call: Crow Wing
is a sandy outwash / lake-rich till basin with **no agricultural tile drainage**,
so the "multipath" mechanism — conceptualized for tile drains — is the wrong
physical model. The right model is `f_route_lake > 0`: land discharge is
partially routed into the lake, the lake buffers and re-releases it via its
outflow law, and the basin discharge is the lake-buffered fraction plus the
direct-channel fraction.

This matches exactly the third row of the diagnostic table in
`crow_wing_river/HANDOFF_lake_support.md`:

> *"KGE up but τ still pinned at bounds → Lake only added flexibility; suspect
> the missing **channelized routing** (lake on main stem)."*

## What `f_route_lake` needs to do (water balance)

The DESIGN doc already specifies the physics — restating the v1 → v1.x increment
for completeness:

- **Lake water balance (per unit lake area):**
  $$\frac{dH_\text{lake}}{dt} = (P - E_\text{lake}) + Q_\text{gw} + f_\text{route} \cdot \frac{a_\text{land}}{a_\text{lake}} Q_\text{land} - Q_\text{out}$$
  where the area ratio mirrors the `Q_gw` exchange convention (volume conserved
  across the land/lake area difference).
- **Basin discharge:**
  $$Q_\text{basin} = (1 - f_\text{route}) \cdot a_\text{land} Q_\text{land} + a_\text{lake} Q_\text{out}$$
- `Q_land` is the land sub-catchment's discharge before the routing fraction is
  diverted. The diverted fraction enters the lake reservoir at the inflow step
  and exits via `Q_out = a(H - H_\text{sill})^b` (unchanged).
- **`f_route_lake = 1`** routes all land discharge through the lake (lake fully
  on main stem, no bypass channel). **`f_route_lake = 0`** is the current v1
  behavior. Crow Wing is somewhere in between, probably high (~0.6–0.9), per the
  lake chain's position on the upper main stem.

## Where the change lands in the code

*Now implemented (commit `b63720a`); the pointers below describe where it landed.
The `NotImplementedError` was replaced by a `[0, 1]` bound check, and the basin-Q
split / routed-inflow term are wired through both the pure-Python `update()` (now
two-pass: land, then lakes) and the JIT.*

- `mnished.py:1486` — the `NotImplementedError` site. The check stays as a
  sanity bound (`f_route_lake ∈ [0, 1]`) but lets non-zero through.
- `_advance_lake()` (or the lake-step equivalent) — adds the routing-inflow term
  using the current step's `Q_land`. Per the design's *"operator splitting,
  start-of-step heads"* convention used for `Q_gw`, apply this consistently
  (probably right after `Q_gw` exchange, before the lake's own outflow step).
- Basin-Q aggregation (`mnished.py:~2574` and the area-fractional discharge
  loop) — reduce the land contribution by `(1 − f_route_lake)` so total mass is
  conserved.
- JIT mirror (`_jit_run`, ~`mnished.py:78+`) — same logic. Today's session also
  surfaced that the JIT cache (`__pycache__/*.nbi`/`.nbc`) goes stale across
  signature changes; the test matrix that exercises lakes in both pure-Python
  and JIT paths is the right place to catch regressions.
- Calibration interface (`run_and_score`, `sub_catchments` override) — likely no
  change needed if `f_route_lake` is sourced from the YAML config (per the
  design's "data-derived, never calibrated" stance).

## Validation

The cleanest end-to-end test is a re-run of today's Crow Wing v2 single-decade
(2001–2010), in this order:

1. **f_route_lake = 0** — must reproduce today's v2 numbers (KGE = 0.626,
   τ_matrix = 974 d). Sanity check that the new code path is a no-op when
   disabled.
2. **f_route_lake = 0.5** — first non-zero test. Expected: KGE improves over
   0.626, τ_matrix begins to relax from 974 d toward till-basin range
   (few-hundred d).
3. **f_route_lake = 0.9** — high-routing test, motivated by Crow Wing's
   on-stem lake-chain topology. Expected: KGE approaches or exceeds v1's 0.703
   *without* multipath, τ_matrix in the few-hundred-d range, lake outflow
   params still interior.
4. **f_route_lake = 1** — degenerate-bypass test. All land discharge through
   the lake. Expected: mass balance closes; KGE may be slightly worse than 0.9
   because some real basins do have a bypass channel.

The headline diagnostic remains **τ_matrix relaxation**, not KGE. If τ relaxes
to ~few-hundred-d at f_route_lake = 0.9, the channelized-routing piece is doing
the right work and Crow Wing can resume on the calibration side.

Mass-conservation unit test: run the basin a synthetic step (`P` impulse, no `E`,
no `Q_gw`) and verify $\int Q_\text{basin}\,dt + \Delta\text{storage} = \int P\,dt$
across both sub-catchments at `f_route_lake ∈ {0, 0.5, 1}`.

## Interface — where does `f_route_lake` come from?

Per `DESIGN_lakes.md`: **data-derived, never calibrated**, ultimately from the
GRASS-fluvial-profiler "prior factory" (DEM → channel network → lake position).
For shipping `f_route_lake > 0` *now*, before the prior factory is in place, the
minimum interface is the existing config block:

```yaml
sub_catchments:
  - name: lake
    kind: lake
    area_fraction: 0.30
    lake:
      outflow_coefficient: 0.005
      sill_storage__mm:    275.0
      outflow_exponent:    1.6667
      gw_partner:          land
      f_route_lake:        0.9          # NEW: data-derived; user-supplied per basin
    initial_conditions:
      lake_storage__mm:    400.0
```

That keeps the calibration-side API simple: the calibration session pins
`f_route_lake` from data (NHD network position + sub-basin contributing-area
fraction draining into the lake), and the model handles it. When the prior
factory ships, the same field gets populated automatically from the DEM
pipeline.

## Open design questions (carried from `DESIGN_lakes.md` and relevant here)

These don't all need to be settled to ship `f_route_lake > 0`, but they shape
how far this slice goes:

1. **On-channel vs. terminal lakes.** Crow Wing is on-channel; Devils Lake /
   prairie potholes are terminal. `f_route_lake = 0` already covers the
   terminal case (no surface inflow). Worth ensuring the new code path doesn't
   break the terminal-lake semantics.
2. **Constant area vs. `A(H)` hypsometry.** v1 assumed constant area; the
   routing increment doesn't strictly require `A(H)` either, but shallow lakes
   in the upper Midwest probably will eventually. Note for the docs, not a
   blocker.
3. **Operator-splitting order** within a timestep: `Q_gw` exchange, then routing
   inflow, then outflow? Or routing inflow first, then `Q_gw`, then outflow?
   Probably matters little at daily resolution but worth committing to a
   convention in the docs.
4. **Multipath separation.** This handoff and today's calibration session both
   make the case that "multipath" (threshold-activated fast pathway) is a
   distinct mechanism from lake routing — not a substitute for it. Worth one
   docs paragraph on when to use which, to prevent the same conflation by
   future contributors who haven't lived this calibration history.

## Why this is the right next slice

From the calibration side, `f_route_lake > 0` unblocks lake-rich basins (Crow
Wing now, Otter Tail, Mississippi headwaters, much of the central Minnesota
basins) and tests the lake architecture in its full form rather than its
deliberately-incomplete v1 form. It is small relative to the v1 lake feature
itself: the underlying lake mechanism, `Q_gw` exchange, JIT mirror, and
calibration wiring all exist on `master`. This is one routing term plus a JIT
mirror plus a docs paragraph.

The deeper "physical priors from DEM geometry" effort (recession priors, lake
network position, contributing-area `f_route` derivation) is the larger v4.0
direction and is correctly out of scope here. What this handoff asks for is just
the code path that lets the user supply a non-zero `f_route_lake` from config —
data-derived for now, prior-factory-derived later.

## References

- `DESIGN_lakes.md` (MNiShed root) — full design, including the
  data-derived-never-calibrated stance for `f_route_lake`
- `HANDOFF_lake_support.md` (MNiShed root) — original science case for lake
  support; the "ET is not the lever" argument
- `crow_wing_river/HANDOFF_lake_support.md` (in `~/dataanalysis/Wickert2026-hydroRaVENS-decadal-optimization/`) — the calibration-side bridge handoff with the v1 test plan and the prediction this handoff confirms
- `mnished.py:1486` — current `NotImplementedError` site
- Issue #19 — GitHub tracking issue
- Today's Crow Wing test runs:
  - `crow_wing_river/backbone_runs/2026-06-25_152500_v6_lake_2001` (v1: lake + multipath)
  - `crow_wing_river/backbone_runs/2026-06-25_155721_v6_lake_2001_no_mp` (v2: lake, no multipath)
- Memory: `project_crow_wing_status.md`

## Constraints (per `~/.claude/CLAUDE.md`)

Same as the other handoffs: stewardship-mode work, explicit per-step
authorization required for commits, tags, pushes, releases. This handoff is
descriptive, not prescriptive — implement what makes sense for MNiShed's
architecture and user community.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
