# Handoff — Lake support for MNiShed

*Written 2026-06-24. From a calibration session in `~/dataanalysis/`.
Documents the science case for adding lake-storage and lake-outflow
physics to MNiShed, the empirical evidence motivating it from a lake-rich
till basin (Crow Wing River at Nimrod), the **key design constraint that
ET should NOT be added as a new lake-specific formulation**, and a
discussion of whether lake ice cover is a separately worthwhile feature.
A thinking-piece for the MNiShed session to decide whether and how to
act on; not an implementation spec.*

## The science case (one paragraph)

A lake-rich basin's hydrograph contains physics MNiShed v3 cannot
represent: direct precipitation onto the open water surface (bypassing
the soil cascade), surface storage on an open water body, and outflow
controlled by a river channel, sill, or dam (asymptotically
$Q \propto H^{5/3}$ for the common case of a friction-controlled river
outlet, but the exponent in nature varies with outlet geometry — see
"Outflow law" below) rather than by porous-medium recession. When a basin like
**Crow Wing River at Nimrod (USGS 05244000, ~30% lake area)** is forced
through the v3 1-reservoir multipath architecture with the regional ET
multiplier, the calibration fails dramatically — KGE = −0.55 with two
parameters at upper bounds. Even after freeing `et_scale` (which
recovers KGE to 0.60), three structural parameters peg at extreme values
(matrix τ = 956 d, multipath τ_mp at upper bound = 98 d, H_thr = 384 mm)
demonstrating that the *storage and outflow physics*, not the ET
treatment, is the missing piece. See the empirical detail below.

## ⚠ Critical design constraint: ET is NOT the lake lever

This is the single most important point and the easiest one to get
wrong.

On the timescales the calibration sees (monthly to annual integrals
of daily Q), **Thornthwaite ET and Penman lake-evaporation are
functionally collinear**:

- Both depend on temperature (Thornthwaite directly; Penman through
  $\Delta(T)$ and $e_s(T)$)
- Both depend on photoperiod / net radiation (Thornthwaite through
  $L(\phi,t)$; Penman through $R_n$, which is correlated)
- Penman's wind, VPD, and radiation extras integrate out on
  monthly-to-annual windows when those forcings aren't supplied
  separately

For a basin with lake fraction $f$:

$$E_{\text{basin}} = (1-f)\,a_{\text{land}}\,E_T \;+\; f\,K\,E_T = E_T\bigl[(1-f)a_{\text{land}} + fK\bigr]$$

That's *the same functional shape* as $E_T$ multiplied by a basin-specific
constant. Calibrating $a_{\text{land}}$ and $K$ separately is
mathematically identical to calibrating a single `et_scale`. A separate
`Lake_ET` module without new forcing inputs (radiation, wind, VPD,
ice-cover dates) is a **parameter-degeneracy machine**. Don't build one.

The right ET treatment for a Lake reservoir is **one extra scalar
multiplier** on the existing Thornthwaite-Chang ET (`et_scale_lake`).
Calibratable, but not a new functional form. No degeneracy.

The Crow Wing empirical test confirmed this. With `et_scale` fixed at
0.755 (BLR's value), KGE = −0.55. With `et_scale` free, KGE = 0.60 and
`et_scale` landed at **0.56** — *lower* than the BLR value, not higher
as a "lake ET undercount" hypothesis would predict. The negative-KGE
recovery was a climate-mismatch fix (Crow Wing at ~47°N is colder than
BLR at ~44°N), not a lake-ET signal. The lake signature surfaced in
storage and outflow parameters, not in ET.

## What IS functionally different about lakes

Three pieces of physics whose functional forms are distinct from the
soil-cascade machinery:

1. **Direct precipitation onto the water surface.** Bypasses canopy
   interception, soil infiltration, and the entire reservoir cascade.
   Goes straight to lake storage and then to outflow. This is a
   parallel-path mechanism, not a scale factor on existing terms.
2. **River-outlet (or sill-controlled) outflow.** A power-law
   stage-discharge relation $Q = a\,(H - H_{\text{sill}})^b$ with
   default $b = 5/3$ (friction-controlled river outlet; Manning) and
   typical range $1.5{-}2.5$ depending on outlet morphology (see
   "Outflow law" below for the physics behind the exponent). Or
   piecewise/threshold for dam-controlled outlets. The *integration*
   of this rating vs porous-medium drainage gives different
   hydrograph shapes — the calibration data can distinguish them.
3. **Lake storage geometry.** A water body sitting on the landscape
   with surface area = lake area, depth proportional to outlet head.
   Direct precipitation contributes to lake stage; outflow depends on
   stage above the sill. Different geometric scaling from porous-medium
   storage.

## Empirical evidence from Crow Wing

USGS 05244000, drainage area 2,335 km² as delineated by the pipeline,
substantial lake cover in upper basin (Park Rapids lakes + the "11 Crow
Wing Lakes" chain), largely unregulated above Nimrod.

| Run | et_scale | KGE | Notable |
|---|---|---|---|
| v6 1-res multipath, `et_scale` fixed at 0.755 (BLR value) | 0.755 fixed | **−0.555** | PDD pegged at 10 (upper bound), H_thr at 1000 mm (upper bound — multipath effectively disabled) |
| v6 1-res multipath, `et_scale` calibrated | **0.56** | **0.597** | matrix τ = 956 d (2× typical till); τ_mp pegged at 98 d (upper bound; vs 13–23 d in pure-till basins); H_thr = 384 mm (3× typical) |

Same architecture works on pure-till basins at KGE ≈ 0.78 (Blue Earth,
Le Sueur) and ≈ 0.62 (Cottonwood — the lower one limited by 1930s
drought / sparse early forcing data). Crow Wing's 0.60 with parameters
stretched far from physical norms is the structural signature of
missing lake physics: the optimizer is force-fitting lake-like
storage and slow-recession behavior by extending matrix and multipath
timescales beyond reasonable substrate values.

## Suggested design: a `Lake` reservoir class (or lake-zone sub-catchment)

Two implementation paths, both compatible with v3.1's sub-catchments
feature:

### Option A — `Lake` reservoir class

A new top-level class distinct from `Reservoir`, with:

- **Surface area** as a parameter (`lake_area_km2` or
  `lake_area_fraction` of basin).
- **Direct-precipitation accounting**: at each timestep, precipitation
  × lake_area is added directly to lake storage (not to the basin
  cascade).
- **ET**: `Thornthwaite-Chang × et_scale_lake` (one scalar parameter).
  No new functional form.
- **Outflow law**: $Q_{\text{out}} = a \cdot (H - H_{\text{sill}})^b$
  when $H > H_{\text{sill}}$, else 0. Three parameters: $a$ (outflow
  coefficient), $H_{\text{sill}}$ (sill/outlet elevation above some
  datum), and $b$ (exponent, default $5/3$ for a friction-controlled
  river outlet but should be calibratable — see "Outflow law"
  subsection below for the physics behind the choice). Optionally a
  piecewise or step-function form for dam-controlled outlets.
- **Optional stage-area-volume curve** if data exists.

This is cleaner conceptually but more code intrusion.

### Option B — Lake-zone sub-catchment

Use v3.1's sub-catchments machinery: declare a "lake zone" sub-catchment
with `area_fraction = lake_area / basin_area`. Within that
sub-catchment, override the reservoir-drainage law to use a weir
formula instead of matrix recession; route precipitation directly to
storage; apply a sub-catchment-specific `et_scale`.

This piggybacks on existing architecture and may be the cleaner first
move. The sub-catchment's reservoir(s) become functionally a Lake
under the hood.

### Either way

- **One ET parameter per lake zone, not a separate formula.** Just
  `et_scale_lake`. Critical.
- **The weir outflow is the functional novelty** — make sure the
  drainage law there is not a matrix recession.
- Backward-compat: when no lake zone / Lake reservoir is declared,
  behavior is unchanged.

### Outflow law — exponent physics

Most natural lakes drain through a river outlet — a meandering channel
where Manning friction sets the stage-discharge rating, not a discrete
sill where critical flow controls upstream stage. For a wide
rectangular outlet channel, Manning gives:

- $u = (1/n) R^{2/3} S^{1/2}$, and for wide channels $R \approx h$
- so $u \propto h^{2/3}$
- $A = w \cdot h \propto h$
- $Q = u\,A \propto h^{5/3}$

This is the default. The **5/3** comes from friction along the outlet
reach, not from a sill control.

**Where other exponents apply:**

| Physics | Exponent | Where it applies |
|---|---|---|
| Manning friction, wide rectangular channel ($u\propto h^{2/3}$, $A \propto h$) | **5/3** | **Default — most natural lake outlets (friction-controlled river reach)** |
| Broad-crested weir, critical flow (rectangular sill) — Bernoulli + Froude=1: $h_c = (2/3)H$, $q = \sqrt{g}\,h_c^{3/2}$ | **3/2** | Rapids/falls/sill at outlet; dam spillways |
| Manning friction, triangular channel ($u\propto h^{2/3}$, $A \propto h^2$) | **8/3** | Outlet channel with strongly triangular/wedge cross-section |
| Sharp-crested V-notch weir, Bernoulli ($u\propto h^{1/2}$, $A \propto h^2$) | **5/2** | Engineered V-notch; rare in nature |
| Dam spillway, engineered | **3/2** with structure-specific $C_d$ | Regulated outlets |

Empirically, stage-discharge ratings for real lake outlets fit
$Q = a(H - H_{\text{sill}})^b$ with $b$ landing in roughly **1.5–2.5**
depending on outlet geometry — 5/3 ≈ 1.67 is in the middle of that
range and is the right central value for natural river outlets;
3/2 = 1.5 is the floor (sill-controlled), 8/3 ≈ 2.67 the ceiling
(strongly triangular channel). Recommendation:

- **Default $b = 5/3$** — friction-controlled river outlet,
  representative of most natural lakes.
- **Suggested calibration range: 1.5–2.5**, with hard bounds 1.0–3.0
  as a sanity envelope.
- **Make $b$ configurable / calibratable**, not hardcoded. A single
  lake-rich basin can drive $b$ to a basin-specific value during
  calibration; cross-basin transfer should be evaluated empirically.
- **A note on naming:** don't call the outflow law a "weir formula"
  in code or docs, because for the default friction-controlled regime
  it isn't one — call it `outflow_law` or `stage_discharge` with $b$
  as the exponent. "Weir" is the right name only for the $b = 3/2$
  sill-controlled and $b = 5/2$ V-notch cases.

## Lake ice — separate question, probably second-order

In northern Minnesota, lakes are ice-covered ~5 months per year
(approximately mid-November through mid-April, with substantial
year-to-year variability and a documented climate-driven shortening
trend). Ice cover changes the hydrology in three ways:

1. **Open-water ET → 0 during ice-on.** But Penman ET in winter is small
   anyway (low net radiation, cold temperatures), so the annual
   integral is dominated by the open-water season. Likely a
   second-order effect on annual water balance.
2. **Precipitation accumulates as snow on ice.** This goes into the
   snowpack module — basin-mean snowpack already does this OK at the
   bulk level.
3. **Ice-melt pulse in spring.** When ice goes out, any liquid water
   that was stored under the ice (from rain-on-snow or partial melt
   events) releases. This is a phase signature in the hydrograph.

### Is ice cover identifiable from streamflow?

A clean test would be:

1. Run a calibrated lake-supported model (whichever option above) on
   Crow Wing.
2. Extract monthly residuals (model_Q − observed_Q) over the calibration
   window.
3. Compute the climatology of those residuals (monthly mean over all
   years).
4. Compare to a published ice-on / ice-off climatology for the Crow
   Wing region (Minnesota DNR has them; MN State Climate Office tracks
   ice cover on major lakes).

If the monthly residual climatology shows a strong dip around ice-out
(April–May) or peak near freeze-up (October–November), and the timing
aligns with ice climatology, **then ice timing carries identifiable
signal in the streamflow record**. Otherwise it's a second-order effect
that can be absorbed into seasonal `et_scale` variation without an
explicit ice module.

My intuition: it'll be borderline detectable. The ice-out spring pulse
is real but in agricultural / forested basins with substantial snow
storage, it gets overprinted by the basin-wide snowmelt response. The
freeze-up suppression of lake ET is real but small (annual integral
dominated by summer). A simple model with no ice module + a Lake
reservoir with the right storage and outflow physics will probably
capture 90% of the signal. Worth the test, not worth a feature unless
the data demands it.

### If ice support IS warranted

It would naturally extend MNiShed's existing `frozen_ground` (FGI)
module:

- Track a per-lake-zone heat balance (or use air-temperature freezing-
  degree-day threshold) to decide ice-on / ice-off dates.
- During ice-on: set the lake-zone open-water ET to zero; route
  lake-surface precipitation through the snowpack module.
- During ice-off: revert to normal lake physics.

This is a minor extension of code; the bigger work is deciding whether
the data justifies it. Test first, build only if needed.

## Suggested triage

1. **Build lake-storage and weir-outflow support** (Option A or B
   above). This is the meaningful new physics. Ship as a feature.
2. **Add `et_scale_lake` as a single scalar.** Don't build a separate
   ET formulation. Document the degeneracy reason.
3. **Defer ice support** until empirical evidence (monthly-residual
   climatology test on Crow Wing) shows the signal is identifiable.

In rough effort terms: probably one PR for lake-storage + weir outflow
+ `et_scale_lake` (≈ multipath-PR scope), and a separate later PR for
ice cover if/when the test motivates it.

## Connections to existing MNiShed features

- **Sub-catchments (v3.1)** is the natural home for Option B.
- **Multipath** is for engineered fast-drainage paths in soil
  reservoirs; not a substitute for weir outflow.
- **`frozen_ground` (FGI)** is the natural extension point for ice
  cover if needed later.
- **`baseflow_Q`** is a constant additive flux; could in principle
  hold a constant lake-outflow contribution but that's a crude
  workaround, not lake physics.

## Open design questions

These need product-design judgment more than coding:

1. **Where does `lake_area_fraction` come from?** User-supplied (from
   surface-hydrology / National Hydrography Dataset / state lake
   inventories), or computed from a basin polygon and a lake-cover
   raster?
2. **One bulk Lake reservoir per basin, or one per major lake?** Most
   basins will have many small lakes; lumping them is probably the
   right first move.
3. **Is the weir coefficient calibrated, or derived from a
   stage-discharge relation if one exists?** Calibrated for first
   pass; refine if the test basin has external data.
4. **How does this interact with stage-based downstream routing?**
   Probably no interaction; lake outflow just feeds the channel as a
   discharge boundary condition.
5. **Storage-water-balance accounting at the basin level.** Direct
   precipitation accounting means the basin water budget needs a
   slightly different bookkeeping; make sure that's wired right and
   tested in unit tests.

## Documentation work needed alongside this feature

The "ET is not the lake lever" claim is the load-bearing reason for
the design choice of `et_scale_lake` over a separate Penman-style
formulation. A one-paragraph assertion in this handoff is not enough;
the argument should be written into the MNiShed architecture guide (or
a dedicated "ET treatment for open water" appendix) as part of the
lake-support PR. Without it, the design decision is one good-faith
"why not Penman, it's more physical?" question away from being
reversed by a future contributor who hasn't seen this context.

The docs analysis should contain:

1. **The algebra of the collinearity.** The basin-ET decomposition
   $E_{\text{basin}} = E_T[(1-f)a_{\text{land}} + fK]$ and the
   observation that $a_{\text{land}}$ and $K$ enter as a single
   effective multiplier — they are mathematically indistinguishable
   without additional information.
2. **The information-content argument.** With only $T$, $P$, and
   photoperiod $L(\phi,t)$ as forcings, Thornthwaite and Penman
   reduce to the same functional shape times a scalar (their
   $T$-dependences differ only weakly through $\Delta(T)$ and
   $e_s(T)$). New separation of $a_{\text{land}}$ vs $K$ requires
   *new forcings*: wind, net radiation, VPD, and/or ice-cover dates.
3. **A worked example.** Compute Thornthwaite and Penman ET for a
   representative MN lake on monthly means; show that after a single
   scalar fit they nearly overlap on the *amount* axis but disagree
   on *phase* (see point 4).
4. **The honest limit of the collinearity argument.** The decomposition
   above assumes $K$ has the same temporal shape as $E_T$. This is
   true at annual timescales but breaks at monthly timescales because
   of lake **thermal mass**: real lake ET peaks in August (not July)
   and is suppressed in spring even after air warms — a phase lag of
   ~1 month relative to air-T-driven $E_T$. A single `et_scale_lake`
   captures the *magnitude* correction but cannot reproduce this
   phase lag. The natural extension when needed is a one-state
   thermal-lag layer on the lake reservoir (lake-water temperature
   tracked with a simple heat budget); flag this as a deferred
   feature.
5. **Literature anchors**: Henderson-Sellers (1986) for the lake heat
   budget framework; Lenters et al. (2005) for Great Lakes
   Penman-vs-T comparisons; Tanny et al. (Lake Kinneret) for
   thermal lag observations; Finch & Hall (2005) for a survey of
   lake-evaporation formulations.
6. **Conclusion**: at the forcing resolution MNiShed uses ($T$, $P$,
   $L(\phi,t)$), a single `et_scale_lake` scalar on Thornthwaite-Chang
   is the *correct* treatment for magnitude. A separate Penman
   formulation without the extra forcings would be a parameter-
   degeneracy machine. Thermal lag is a separate, later feature when
   phase fidelity matters.

This belongs in the docs, not buried in this handoff, because future
designers will need to see the formal argument to trust the design
choice.

## Where to find the full reasoning

Detailed memory entry: 
`/home/awickert/.claude/projects/-home-awickert-dataanalysis/memory/feedback_lake_support_needs.md`

That entry contains the empirical evidence, the diagnostic signatures
to look for in lake-rich basins, the degeneracy argument for ET, and
the separate climate-gradient finding that `et_scale` doesn't transfer
across the Upper Midwest's latitudinal climate gradient (BLR's 0.755
at 44°N is too high for Crow Wing at 47°N — methodological correction
worth carrying forward independently of lake support).

## Constraints (per `~/.claude/CLAUDE.md`)

Same as the other handoffs: stewardship-mode work, explicit per-step
authorization required for commits, tags, pushes, releases. This
handoff is descriptive, not prescriptive — implement what makes sense
for MNiShed's user community, in whatever form makes sense, when the
time is right.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
