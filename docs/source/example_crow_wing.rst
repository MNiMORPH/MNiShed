Worked Example: Crow Wing River
===============================

The Crow Wing River (north-central Minnesota; ~2335 km², sandy outwash, mixed
forest, many lakes and wetlands) is MNiShed's **process-structure** example.
Where the Cannon example (:doc:`calibration`) is about calibration and
identifiability, this one shows how three optional process modules — **lakes**,
**parallel sub-catchments**, and **growing-degree-day phenology** — combine, with
a physical reason for each, to reproduce a northern-forest hydrograph across
eight decade-long windows (1931–2020).

The runnable files are in ``examples/crow_wing/``::

    cd examples/crow_wing
    conda activate mnished-jit                  # numba JIT + spotpy
    python calibrate.py sceua [reps]            # best-fit (SCE-UA)
    python calibrate.py dream [reps] [iid|ar1]  # posterior / UQ (DREAM)

``calibrate.py`` is the *same* generic, config-driven runner as the Cannon
example — no per-basin Python. The basin is defined entirely by the ``target:``
mappings in ``params.yml`` (see :class:`~mnished.Calibrator`).

The three modules, and why each is here
---------------------------------------

**Lakes** (:ref:`lakes`).
  A single ``kind: lake`` sub-catchment is a series outlet element: it stores
  water above a sill and releases it as
  :math:`Q = a\,(H - H_\text{sill})^{5/3}` (Manning wide-channel spillover),
  exchanging groundwater bidirectionally with its partner land zone. It buffers
  the part of the basin that drains through it.

**Sub-catchments** (:ref:`parallel-sub-catchments`).
  Two parallel land zones — ``direct_land`` (drains straight to the gauge) and
  ``lake_basin_land`` (drains *through* the lake, fraction ``f_route_lake``) —
  each a two-reservoir cascade: a fast "soil" store over a slow "groundwater"
  store. Splitting the surface is what lets one part of the basin be buffered by
  the lake while the other responds directly.

**Phenology** (the vegetation coefficient in :doc:`configuration`).
  A growing-degree-day vegetation coefficient (:math:`K_c`) on ET. The single
  calibrated knob is ``leafout_GDD``; the rest of the curve is fixed from
  regional phenology.

**Frozen ground** supplies the snowmelt freshet pathway: when the frozen-ground
index exceeds ``fdd_threshold`` the soil's top reservoir sheds all snowmelt as
direct runoff, with infiltration resuming on thaw. For this sandy basin — which
only sheds water when frozen — that is the dominant overland-flow control. (The
PDM saturation-excess store is disabled: it fired on summer rain and was rejected
by calibration.)

The freshet problem, and why phenology fixes it
-----------------------------------------------

Thornthwaite ET is temperature-only: it ramps up with spring warmth regardless of
whether the canopy has leafed out. In a northern mixed forest the canopy does not
transpire until **mid-to-late May**, so a temperature-only ET *evaporates the
April snowmelt freshet* that should reach the gauge — spring under-produces and,
through the low ``et_scale`` needed to close the annual balance, fall
over-produces.

The growing-degree-day :math:`K_c` holds ET near ``dormant_Kc`` until thermal-time
leaf-out (``leafout_GDD``), then ramps to ``full_Kc`` and back down through autumn
senescence. Calibrating ``leafout_GDD`` lets the data place green-up: for
north-central Minnesota it settles near **~200 GDD ≈ a late-May onset**, matching
regional phenology — *later* than the model's generic ~100-GDD (southern-MN)
default.

.. note::

   **The GDD leaf-out prior is latitude-dependent.** Crow Wing (~46–47°N)
   calibrates to a later green-up than the generic default tuned to the Cannon
   basin further south. A regional spring-index (Schwartz et al. 2013; see
   :doc:`references`) is the natural source for the ``leafout_GDD`` prior — a
   candidate for the phenology defaults discussion (MNiMORPH/MNiShed#26).

The seasonal mass-balance diagnostic
------------------------------------

The phenology fix was *found*, not assumed. A spring-under / fall-over residual
on this basin first looked like a lake-routing problem; what reframed it as an
ET-phasing problem was a **seasonal mass-balance decomposition** — the method
shipped as :class:`~mnished.SeasonalMassBalance` (see
:doc:`calibration`). For a calibrated run it tabulates, per season, the basin
water balance (P, ET, storage change, discharge) with **discharge split by
source** (fast/event, slow/baseflow, lake outlet) against observed Q::

    result = run_and_score('crow_wing_config.yml', ..., store_fluxes=True)
    print(SeasonalMassBalance(result.buckets).report())

On Crow Wing it showed the April freshet being *evaporated* (spring ET ≈ P at
melt) rather than mis-routed — distinguishing an ET-phasing error from a routing,
storage, or recession error. The method generalizes to any basin; it is the
reusable tool behind this example's narrative.

Running it and what to expect
-----------------------------

SCE-UA reaches a composite ``KGE_logKGE`` ≈ **0.77** over the eight decades, with
``leafout_GDD`` at a physically-correct late-May green-up and the melt factor at a
forest-physical value:

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - parameter
     - value
     - note
   * - ``leafout_GDD``
     - ~205 GDD
     - green-up ~late May (physical for ~46–47°N)
   * - ``PDD_melt_factor``
     - ~2.2
     - mm SWE °C⁻¹ day⁻¹ — forest-physical
   * - ``et_scale``
     - ~0.83
     - free (no water-balance rescaling)

Calibrating ``leafout_GDD`` (vs. fixing it at the generic default) tightens the
spring freshet — top-20 observed-peak mod/obs goes from ~0.67 to ~1.04 — and pulls
the fall recession into line (SON ~1.2 → ~1.0).

.. note::

   **This is a research example, not a polished operational calibration.** Two
   residuals remain visible: a winter↔fall trade-off (``KGE_logKGE`` is fairly
   flat across the seasonal shape, so runs can swap a slightly high winter for a
   slightly high fall at nearly the same score — a seasonally-weighted objective
   would resolve it), and a gentle groundwater recession whose shape contributes
   part of the residual fall flow. They are documented rather than hidden.

See ``examples/crow_wing/README.md`` for the full parameter set, the eight decade
windows (the 1981–1990 decade is omitted — its observed discharge is ~8% present),
and the per-season mod/obs table.
