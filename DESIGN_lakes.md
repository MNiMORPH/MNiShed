# Design ideas — lakes in MNiShed

*Thinking-piece, 2026-06-24. Not a spec. Captures ideas for representing lakes
(and open water generally) in MNiShed, in the parsimony-first, field-data-driven
spirit of the model: **minimum free parameters, maximum use of basic field data,
explain as much of nature as possible.** Researched against PRMS for inspiration;
links to the MNiMORPH open-water Penman implementation.*

---

## ✅ Committed design (2026-06-25 session)

*This section is the agreed design and supersedes parts of the original
exploration below. The exploration is **retained as history**, but two of its
premises are now changed:*

1. *The lake is **not** a parallel land sub-catchment — it is a **series element
   at the basin outlet**, plumbed into the land subsurface reservoir by a
   bidirectional groundwater flux.*
2. *The lake does **not** take a separate open-water (Penman) evaporation forcing.
   It reuses the basin-wide Thornthwaite `et_scale`. (See
   `HANDOFF_lake_support.md`: at MNiShed's forcing resolution — `T`, `P`,
   photoperiod — Penman and Thornthwaite are collinear, so a separate lake-ET
   formulation is a parameter-degeneracy machine. The Crow Wing test confirmed it:
   freeing `et_scale` sent it **lower**, 0.56, not higher.)*

### Topology — series element + bidirectional groundwater exchange

A lake is a lumped storage **in series at the outlet**, coupled to the land
**subsurface reservoir** `h_s` by a bidirectional, head-driven groundwater flux.
Three nodes (subsurface `h_s`, lake `H_lake`, outlet), three fluxes:

```
   recharge (soil cascade)              P − E (direct, on lake surface)
        │                                    │
        ▼                                    ▼
   ┌──────────┐        Q_gw            ┌──────────┐
   │ subsurf. │◄══════════════════════►│   LAKE   │
   │   h_s    │  a_sub·sgn(Δh)|Δh|^b_sub│  H_lake  │
   └────┬─────┘   Δh = h_s − H_lake    └────┬─────┘
        │ Q_bf = a_sub·h_s^b_sub            │ Q_out = a·(H_lake−H_sill)^(5/3)
        ▼                                   ▼
   ══════════════ OUTLET / gauge  Q ══════════════►
```

- **GW → outlet** (unchanged from current model):
  `Q_bf = a_sub · h_s^{b_sub}` — the existing subsurface recession.
- **Lake → outlet**: `Q_out = a · max(H_lake − H_sill, 0)^{5/3}` — surface
  stage–discharge (Manning, friction-controlled river outlet).
- **GW ↔ lake**: `Q_gw = a_sub · sign(h_s − H_lake) · |h_s − H_lake|^{b_sub}` —
  the **same recession coefficient and exponent** as GW→outlet, made bidirectional
  by the `sign(·)|·|^b` form (reduces to linear Darcy when `b_sub = 1`).

Surface inputs: recharge to `h_s` from the soil cascade; direct `P − E` onto the
lake surface, with `E` from the global `et_scale` (no separate lake-ET term).

### Parameters — five collapsed to two calibrated

| Parameter | Decision | Why |
|---|---|---|
| `b` (lake outflow exponent) | **fixed 5/3, not calibrated** | Manning wide-rectangular river outlet; an outlet property, unaffected by lake area. Assumes a **linear storage→stage conversion** (constant-area); lake spreading is unmodeled, not a distortion of the exponent. |
| open-water ET | **reuse global `et_scale`** (no new param) | Penman/Thornthwaite collinear at this forcing resolution; add a lake-specific scalar only if open-water residuals demand it. |
| `C` (lake↔GW conductance) | **eliminated** — reuse `a_sub`, `b_sub` | Same aquifer drains to stream and lake; the GW exchange rides the existing recession law. Inherits the Ksat / recession-from-geometry prior automatically. |
| `a` (lake outflow coeff) | **calibrated** | Lumped **effective** coeff: hydraulics × storage→stage × storage→area translation — *not* a clean Manning-n value. |
| `H_sill` (sill / dead-pool) | **calibrated** | Threshold in *conceptual storage units*; sets the dead pool that keeps exchanging `P`, `E`, GW below the outlet. |
| `lake_area_fraction` | **from data** (NHD / lake inventory), not fit | Scales direct `P − E`. |
| `f_route_lake` (fraction of basin routed *through* the lake) | **held at 0 for first runs**; later **data-derived** (DEM contributing area), never calibrated | Channelized surface inflow `Q_in`. 0 = lake hydrologically disconnected from river inflow. |

Net: **two new calibrated parameters (`a`, `H_sill`)** for storage + dead pool +
surface outflow + bidirectional groundwater exchange.

### v1 scope — surface inflow deferred (`f_route_lake = 0`)

For the first implementation the lake is **hydrologically disconnected from river
inflow**: its only water sources are direct `P − E` on its surface and the
bidirectional groundwater exchange `Q_gw`. This is deliberately *wrong* (a real
lake's mass is dominated by channelized inflow) but useful:

- It is a **clean controlled test of the `Q_gw` mechanism alone.** With no surface
  routing, any buffering improvement on a lake-rich basin (Crow Wing) comes purely
  from the groundwater capacitor — so the first run tells us whether the
  lake↔subsurface coupling, by itself, relaxes the stretched substrate τ. If it
  does, the GW pathway is doing real work; if it barely moves, that localizes the
  action to channelized routing (the bigger piece below).
- The lake is **not inert** even disconnected: the aquifer still drains partly
  into it and the lake re-releases slowly, so the seasonal store/release buffering
  is present via the subsurface route — just missing the (larger) main-stem flood
  routing.

The `f_route_lake` parameter is included (default 0) so the hook exists, but the
**routing machinery and the lake-network position are deferred** to the larger
planned work: **channel network → drainage density → soil Ksat priors** (the
GRASS-fluvial-profiler "prior factory", `DESIGN_recession_priors_from_geometry.md`
in that repo). That same DEM→network pipeline that delivers recession priors and
contributing areas is where lakes get inserted into the drainage network with a
data-derived `f_route_lake` (and main-stem / tributary / terminal classification).
Lake routing and the recession-from-geometry work are one body of work, not two.

### Key decisions and why

- **No partition coefficient** between GW→lake and GW→outlet. The split is set
  dynamically by the relative head differences (`h_s − H_lake` vs `h_s − 0`); a
  fixed fraction would freeze it and kill the self-regulating flow-through
  behavior. (The nearest knob is an optional dimensionless conductance *ratio*
  `γ` on `Q_gw`, default **1, off** — a ratio `0..∞`, not a partition `0..1`.)
- **`a` and `H_sill` both free, not linked.** They share the storage→elevation
  factor α (`a = k·α^b`, `H_sill = η_sill/α`), but each also carries an
  independent unknown (outlet conductance `k`; real sill height `η_sill`), so
  collapsing them would smuggle in an unjustified assumption. With `b` fixed, the
  rating is a clean, well-identified 2-parameter curve.
- **Conceptual datum.** MNiShed storages are conceptual depths, not surveyed
  elevations; `H_lake`, `H_sill`, and `h_s − H_lake` live on conceptual scales, so
  `a` and `C` are effective (absorb the storage↔elevation translation). Only `b`
  and `et_scale` carry usable physical priors.
- **Reduces to the current model when no lake is declared** — no `H_lake`, no
  `Q_gw`, no `Q_out`; `h_s` drains to the outlet exactly as today. Bit-identical
  backward compatibility; the lake is purely additive.

### What it buys — lake buffering

Adding the lake gives the groundwater reservoir a **second exit**, and water that
takes it re-emerges slowly through the lake's outlet — so the lake is a storage
element (capacitor) inserted *into* the subsurface flow network, not a bucket
catching only rain. The single bidirectional `Q_gw` flips sign on its own:

- **wet / high water table** (`h_s > H_lake`): GW flows *into* the lake — it fills
  and **attenuates the flood peak**;
- **dry / recession** (`h_s < H_lake`): the flux reverses — the lake feeds the
  aquifer and spills its sill, **sustaining baseflow**.

This store-in-spring / release-in-summer is the classic lake signature, emergent
from one term with no seasonal switching. It is the direct fix for the Crow Wing
diagnosis (`HANDOFF_lake_support.md`): with no lake, the optimizer manufactured
this slow store-and-release by stretching matrix τ to 956 d and multipath τ to its
bound. **Validation target: KGE recovers *while* those substrate timescales relax
back toward till-basin values** — proof the lake does the right physical work
rather than merely adding flexibility.

---

## Why lakes are different from land

A lake is not a soil column. Four things make it distinct, and each maps onto
something MNiShed already has — or nearly has:

1. **Open-water evaporation, not ET.** A lake evaporates at (near) the potential
   rate set by energy and the atmosphere — no stomatal/soil-moisture limit, no
   Thornthwaite land-ET surrogate. It needs an *open-water* potential
   evaporation forcing (Penman), distinct from the land `ET for model` column.
2. **Direct precipitation onto the surface.** Rain on the lake goes straight to
   lake storage; it does not infiltrate.
3. **A surface area / fraction of the basin.** Both evaporation and direct
   precipitation scale with lake area, so a lake carries an area weight — exactly
   what sub-catchment `area_fraction` already is.
4. **Outflow set by the outlet geometry, not by drainage recession.** Lake stage
   drives outflow through a weir/orifice/level-pool relationship — a *measurable*
   structure, not a fitted timescale.

## Lake water balance (per unit lake area)

```
dH_lake/dt = P  +  Q_in/A_lake  −  E_open  −  Q_out(H_lake)/A_lake
```

where `P` is precip on the surface, `Q_in` is catchment inflow, `E_open` is
open-water evaporation, and `Q_out` is stage-driven outlet discharge. This is the
classic level-pool (lumped reservoir) lake.

## What PRMS does (for inspiration)

PRMS represents a lake as a distinct **lake-type HRU** wired into the stream
network, with a dedicated **stream-and-lake flow-routing** module
(`muskingum_lake`; documented in USGS TM 6-B8). The lake HRU receives
precipitation on its surface and loses **open-water evaporation**, takes inflow
from upstream segments, and computes **outflow by one of several user-selected
methods**, including:

- **Modified-Puls / level-pool routing** — storage–outflow relation solved on
  the continuity equation (the standard reservoir-routing approach).
- **Linear reservoir** — outflow linear in storage (a single coefficient).
- **Broad-crested weir** — `Q ∝ (H − H_crest)^{3/2}`.
- **Gate / orifice opening** — controlled-flow equation.
- **Measured / specified** stage–discharge (or stage–area–volume) relationship.

(Exact method names/numbering should be verified against the USGS TM 6-B8
documentation before quoting them anywhere user-facing.) HBV and Raven take a
similar line: a lake is a storage with potential-rate evaporation and a rating
curve for outflow; Raven in particular exposes explicit reservoir/lake elements
with management rules.

The take-away for MNiShed: **none of this needs distributed physics.** A lake is
a lumped storage with (a) open-water evaporation, (b) surface precipitation, and
(c) a stage-driven outlet rule. MNiShed already has most of the pieces.

## How a lake maps onto MNiShed's existing machinery

This is the appealing part — a lake is close to a special case of structures we
already built:

| Lake feature | MNiShed mechanism that nearly covers it |
|---|---|
| Area fraction of the basin | sub-catchment `area_fraction` |
| Own forcing (open-water evap) | **per-sub-catchment forcing** (the deferred `forcing:` hook — currently `NotImplementedError`) |
| Stage-driven weir outflow `Q ∝ H^{3/2}` | **power-law recession** with `b ≈ 1.5` (broad-crested weir is exactly a power law in head!) |
| Dead storage below the outlet | **threshold junction** (`H_threshold`) — water below the sill does not discharge |
| Level-pool storage | a reservoir with `H_threshold` + power-law outflow |
| Direct precip to storage, no infiltration | a reservoir that takes `P` directly and does not subtract land ET |

So the **minimal** lake is: *a sub-catchment whose forcing is open-water
evaporation instead of land ET, containing a single threshold-power-law reservoir
(weir).* Strikingly, that is mostly assembly of existing parts — the one genuinely
missing capability is **per-sub-catchment forcing**, which is already on the
radar as a v4.0/post-v3.1 candidate. Lakes are a strong *motivating use case* for
that feature.

What's still missing or worth adding, in increasing order of effort:

1. **Open-water evaporation forcing** — an optional input column (or per-zone
   forcing file) carrying Penman PE; the lake sub-catchment uses it in place of
   the land-ET term. Direct precip already flows to the top reservoir.
2. **A weir/threshold outlet preset** — convenience so a user writes
   `outlet: broad_crested_weir` and gets `b = 1.5` + a threshold at the sill
   stage, with the weir coefficient from outlet width/geometry rather than a
   fitted τ.
3. **A `Lake` (or `OpenWater`) sub-catchment kind** — explicit, so the
   bookkeeping (surface precip, no infiltration, open-water E, area weighting)
   and the diagnostics are clear, rather than a land sub-catchment in disguise.
4. **Stage–area–volume (hypsometry)** — if lake area varies with stage (most do),
   evaporation and surface precip should scale with `A(H)`. This is the first
   place a lake genuinely needs more than the existing reservoir math. Often
   derivable from a bathymetric survey or a DEM + shoreline.

## Open-water evaporation and the Penman repo

`MNiMORPH/TerraClimate-potential-open-water-evaporation` already computes
potential open-water (lake/ocean) evaporation via the **Penman** equation, in
modular Python (`penman.py`, `net_radiation.py`, `vapor_pressure.py`,
`wind_shear_velocity.py`, `atmospheric_parameters.py`, `sunpos.py`), driven by
TerraClimate gridded meteorology; GPL-3.0. This is exactly the forcing a lake
module needs.

Options for how it relates to MNiShed (these are not exclusive):

- **(A) Make it a freestanding, installable package.** A focused, reusable
  open-water-evaporation library is valuable on its own (and citable). It does
  one physical thing well — squarely in the MNiShed philosophy of small,
  field-data-driven components. *Recommended regardless.*
- **(B) MNiShed accepts open-water PE as a forcing input.** Lowest coupling:
  MNiShed reads an `Open-water evaporation [mm/day]` column (analogous to the
  `datafile` ET method) for a lake sub-catchment, computed by whatever tool the
  user likes — including (A). No hard dependency.
- **(C) MNiShed depends on it (optional extra).** `pip install mnished[lakes]`
  pulls the Penman package; MNiShed can compute open-water PE internally from the
  same met inputs it already reads. Convenient, more coupling.

**Lean:** do (A) — cut the Penman code as its own package — and wire MNiShed via
(B), with (C) as an optional extra later. That keeps MNiShed parsimonious and the
evaporation physics reusable and independently testable. The Penman method (vs.
Thornthwaite for land) also nicely matches the "use real physics where the field
data supports it" goal.

## Why this fits the v4.0 / "min params, max field data" thesis

A lake is the cleanest possible example of the philosophy: its parameters are
**measurable field data, not free knobs** —

- surface area / fraction → maps, shoreline + DEM, remote sensing;
- outlet sill stage and width → survey of the outlet structure;
- weir exponent → physics (`3/2` for a broad-crested weir);
- evaporation → Penman from met data.

So a lake can add real physical structure to a basin model while adding **zero or
near-zero calibrated parameters**. That is the same move as the K-threshold /
identifiability work (substrate properties from soil-hydrology data set what's
constrainable) and the same move per-sub-catchment forcing enables. Lakes,
substrate-K identifiability, and field-data priors are three faces of one v4.0
direction: *physics and field data set the structure and the bounds; calibration
only touches what the data can actually constrain.*

## Open design questions

1. **On-channel vs. off-channel/terminal lakes.** An on-channel lake routes
   inflow→outflow (level-pool); a terminal/closed lake (no outlet) only balances
   P − E and integrates to a stage (e.g. Devils Lake, prairie potholes). Both are
   interesting in the upper-Midwest basins; the terminal case is the simpler
   water balance but needs the evaporation right.
2. **Where does surface area come from, and does it vary with stage?** Constant
   area is the minimal model; `A(H)` hypsometry is the honest one for shallow
   lakes/wetlands. Decide whether v1 is constant-area.
3. **Is a lake a sub-catchment, or its own object?** Reusing sub-catchments is
   parsimonious but conflates "parallel land zone" with "open water." A distinct
   `Lake`/`OpenWater` kind is clearer but more code. (Leaning: start as a
   sub-catchment with open-water forcing; promote to its own kind if the
   bookkeeping demands it.)
4. **Coupling depth to the Penman package** — input column (B) vs. optional
   dependency (C).
5. **Wetlands / depression storage as the same machinery?** Prairie-pothole
   depression storage (cf. PRMS surface-depression storage) is arguably the same
   lumped open-water-with-threshold-outlet object at smaller scale — worth
   keeping in view so we don't build two things.

## Suggested triage

- **First, freestanding Penman package (A)** — useful immediately, independent of
  the rest, on-philosophy.
- **Then, lake-as-sub-catchment with open-water forcing (B)** — which rides on the
  *per-sub-catchment forcing* feature; lakes are a good reason to prioritize that.
- **Defer** stage-area-volume hypsometry and a dedicated `Lake` class until a real
  basin needs them.

## Sources

- USGS PRMS — stream and lake flow routing (TM 6-B8):
  https://www.usgs.gov/data/dynamic-parameter-water-use-stream-and-lake-flow-routing-and-two-summary-output-modules-and
- PRMS-IV documentation (TM 6-B7): https://pubs.usgs.gov/tm/6b7/pdf/tm6-b7.pdf
- PRMS software: https://www.usgs.gov/software/precipitation-runoff-modeling-system-prms
- Open-water Penman implementation:
  https://github.com/MNiMORPH/TerraClimate-potential-open-water-evaporation
