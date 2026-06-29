# MNiShed roadmap

A prioritized timeline of upcoming work, ordered by when it should happen. This
is a living document; the authoritative per-item detail lives in the linked
GitHub issues and the `notes/` design/handoff docs.

**Status (2026-06-29):** v3.2.0 released (Zenodo concept DOI 6787390). The
v3.3.0 content is staged in `CHANGELOG.md`'s `[Unreleased]` section — `mnished.io`
input contract, photoperiod senescence + `leafout_GDD_from_date`,
`KGE_logKGE_seasonal` objective, the numpy≥2 floor, and the #36 water-balance /
flux-partition fixes.

---

## Phase 1 — Ship v3.3.0 (immediate)

Each irreversible/outward step needs explicit authorization.

1. **Push** the local `master` commits.
2. **Pre-release `/code-review`** of the `[Unreleased]` diff; fix findings. (A
   review must precede the tag — it caught a real bug before the v3.2.0 DOI.)
3. **Close #37** (code complete).
4. **Version bump** 3.2.0 → 3.3.0 (`mnished/_version.py`, `CITATION.cff` version +
   `date-released`) and **roll the CHANGELOG** `[Unreleased] → [3.3.0]` with date +
   compare link; open a fresh `[Unreleased]`.
5. **Tag `v3.3.0` → draft release notes → approval → publish** the GitHub Release
   (mints the Zenodo version DOI).

## Phase 2 — Round out the phenology story (post-release; no external blockers)

6. **Multi-basin latitude-transfer validation** of photoperiod senescence — the
   real test the single-basin Crow Wing wash could not do (it is structurally
   blind to cross-latitude transfer). A natural job for `mnished-builder` across
   the study basins at different latitudes.
7. *(optional)* Configurable per-season weights for `KGE_logKGE_seasonal`, if a
   basin ever needs to privilege a season (default stays equal-weight).

## Phase 3 — Priors-from-data arc

**Gated** on the GRASS `r.fluvial.channelheads` tool emitting the full channel
network (not just heads). The GIS producers live in
`MNiMORPH/GRASS-fluvial-profiler` / `rivernetworkx`; MNiShed stays GIS-free.

8. **#29** hillslope / flow-routing geometry extraction → **#31** recession τ and
   b priors from geometry + soil. *(This is where the "where does the τ,b
   transfer-function live" decision lands.)*
9. **#32** diffusion-length / matrix-τ identifiability diagnostic + prior.
10. **#30** POLARIS soil-property consumption (producer `r.in.polaris` exists);
    **#22** land-cover-derived degree-day-factor priors.
11. **mnished-builder** recession-prior stage that orchestrates 8–10 and folds the
    derived priors into the generated MNiShed config.

## Phase 4 — Lakes from terrain

**Gated**; needs newer algorithms.

12. **#33** terrain-derived lake routing (`f_route_lake` / lake-network position),
    **#34** lake outflow coefficient from outlet geometry, **#19** generalize
    `SubCatchment` into a basin-element abstraction.

## Phase 5 — v4.0.0 (next major)

13. **#18** drop the K=1 special-casing → uniform per-sub-catchment state *(wanted
    sooner)*, **#20** make in-process SPOTPY the default engine, **#24** rename
    `decades:` → `windows:`.

## Standalone (anytime)

14. **#26** NDVI/phenology-driven ET factor (the data-rich complement to the
    parsimonious phenology); **#13** nested-gauge analysis utilities.
