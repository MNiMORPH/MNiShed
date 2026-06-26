Calibration
===========

MNiShed calibration is built around :func:`~mnished.run_and_score`,
which runs the model with a given parameter set, applies optional Nash-cascade
routing and a scoring window, and returns a :class:`~mnished.CalibResult`
containing all goodness-of-fit metrics.  The recommended workflow drives this
**in-process** via :class:`~mnished.Calibrator`, a config-driven, build-once
model setup (see *In-process calibration with* ``Calibrator`` below).  An
external optimizer (e.g. `Dakota <https://dakota.sandia.gov>`_) can be used
instead via a thin driver that reads parameter values from the optimizer's
output file and calls ``run_and_score``.

.. contents:: On this page
   :local:
   :depth: 2

run_and_score()
---------------

.. autofunction:: mnished.run_and_score
   :no-index:

CalibResult
-----------

.. autoclass:: mnished.CalibResult
   :no-inherited-members:
   :no-index:

The fields of ``CalibResult`` are:

.. list-table::
   :widths: 20 15 55
   :header-rows: 1

   * - Field
     - Type
     - Description
   * - ``score``
     - float
     - Composite goodness-of-fit score from the metric chosen in
       ``metric=``.  Higher is better (all metrics return higher values
       for better fits).  ``np.nan`` if no observed discharge falls in
       the scoring window.
   * - ``aic``
     - float
     - Akaike Information Criterion computed on log-transformed discharge.
       Lower is better.  See :ref:`aic-model-comparison`.
   * - ``bfi_obs``
     - float
     - Observed baseflow index (Eckhardt recursive-filter method,
       α = 0.98, BFI_max = 0.80).
   * - ``bfi_mod``
     - float
     - Modeled baseflow index, computed identically to ``bfi_obs``.
   * - ``kge_logfdc``
     - float
     - KGE on the log-transformed flow-duration curve.  Sensitive to
       reservoir partitioning and long-term storage behaviour.
   * - ``fdc_obs``
     - pd.Series
     - Observed FDC at 99 exceedance percentiles (log-space).
   * - ``fdc_mod``
     - pd.Series
     - Modeled FDC at 99 exceedance percentiles (log-space).
   * - ``final_states``
     - dict
     - End-of-run storage: reservoir depths (mm), snowpack SWE, frozen-ground
       index, and carried deficit.  Flat for a single-zone basin; nested under
       ``'sub_catchments'`` for a partitioned basin.  Pass as
       ``initial_states=`` in the next decade run to chain simulations.
   * - ``buckets``
     - :class:`~mnished.Buckets`
     - The fully initialised and run model object.  Access
       ``result.buckets.hydrodata`` for the full simulated time series.

Goodness-of-Fit Metrics
------------------------

The ``metric`` argument to :func:`~mnished.run_and_score` selects the
composite score maximised during calibration and returned as
``CalibResult.score``.  The AIC is always computed and returned separately
regardless of the chosen metric.

Available metrics
^^^^^^^^^^^^^^^^^

.. list-table::
   :widths: 30 15 45
   :header-rows: 1

   * - ``metric`` key
     - Score range
     - Best suited for
   * - ``KGE``
     - −∞ to 1
     - Overall hydrograph shape; balanced sensitivity to timing, magnitude,
       and variability.
   * - ``NSE``
     - −∞ to 1
     - Peak-flow emphasis; dominated by high-flow events.
   * - ``logKGE``
     - −∞ to 1
     - Low-flow and baseflow behaviour; down-weights large peak events.
   * - ``KGE_logKGE``
     - −∞ to 1
     - Balanced calibration target: average of KGE and logKGE.
       Recommended for most temperate catchments.
   * - ``KGE_logKGE_logFDC``
     - −∞ to 1
     - As above but also includes the flow-duration curve.  Useful
       when the full range of flows is of interest.
   * - ``KGE_logKGE_logFDC_BFI``
     - −∞ to 1
     - As above plus a baseflow-index match term
       (:math:`1 - |\text{BFI}_\text{mod}/\text{BFI}_\text{obs} - 1|`),
       averaged over all four.  Use when the groundwater partitioning
       (baseflow fraction) must be matched alongside the hydrograph shape.
   * - ``logKGE_logFDC_BFI``
     - −∞ to 1
     - Low-flow focus: averages logKGE, the log-FDC KGE, and the
       baseflow-index match, with no high-flow KGE term.  Use for
       baseflow- and drought-dominated studies.

**Recommendation:** ``KGE_logKGE`` is the default and recommended metric for
calibrating multi-reservoir models in temperate climates.  It weights peak and
low flow equally and does not over-emphasise the largest flood events the way
NSE or raw KGE do.  When the flow-duration curve (groundwater contribution
to dry-season flow) is a specific focus, use ``KGE_logKGE_logFDC``.

Baseflow index (BFI) is reported in all ``CalibResult`` objects as a
diagnostic regardless of the chosen calibration metric.  Use it to verify
that the modeled groundwater contribution is physically plausible.

.. _aic-model-comparison:

AIC and model comparison
^^^^^^^^^^^^^^^^^^^^^^^^

The Akaike Information Criterion allows comparison of models with different
numbers of free parameters:

.. math::

    \text{AIC} = n \ln\!\left(\frac{\text{RSS}_{\ln Q}}{n}\right) + 2k

where :math:`n` is the number of scored days, :math:`\text{RSS}_{\ln Q}` is the
residual sum of squares on log-transformed discharge, and :math:`k` is the count
of free calibration parameters.  Lower AIC indicates a better trade-off between
fit quality and model complexity.

The AIC is evaluated on **log-transformed discharge** because the residuals of
raw discharge are dominated by large flood events, which gives a misleading
picture of overall model skill.  Log-transformation gives roughly equal weight
to all flow magnitudes.

**Counting free parameters (k):**
Only parameters that are actively calibrated (not fixed) contribute to :math:`k`.
In ``run_and_score``, :math:`k` is accumulated automatically for each parameter
that is passed as a non-trivial argument (e.g. each element of ``recession_coeff`` and
``f_to_discharge`` that is free).  Fixed structural choices (e.g. ``b_deep = 1``
by convention, not optimization) do not add to :math:`k`.

Some arguments are lists that may hold a mix of fixed and calibrated entries
— the junction parameters (``leakance_R``, ``H_threshold``), the multipath
drain (``multipath_threshold``, ``multipath_timescale``), and the
``recession_exponents``.  For these, the free count cannot be inferred from
the values alone, so it is given explicitly through a companion
``*_calibrated`` argument.  Set ``leakance_R_calibrated``,
``H_threshold_calibrated``, ``multipath_calibrated``, and
``recession_exponents_calibrated`` to the number of entries in the
corresponding list that the optimizer is free to vary; only that many are
added to :math:`k`.

The ET scaling factor added by ``enforce_water_balance='water-year'`` is a
hidden degree of freedom: one multiplier per water year is computed from the
data, effectively adding :math:`N_{\text{years}}` free parameters.
``enforce_water_balance='global'`` computes a single multiplier for the full
record, adding only 1 degree of freedom (not counted in :math:`k` because it is
derived analytically from the data, not optimized).  For multi-model AIC
comparisons, use ``enforce_water_balance='global'`` to avoid this hidden
per-year penalty.

**Interpreting ΔAIC:**

.. list-table::
   :widths: 20 60
   :header-rows: 1

   * - :math:`\Delta\text{AIC}`
     - Interpretation
   * - 0–2
     - Substantial support for both models; the simpler one is preferred
       on parsimony grounds.
   * - 2–4
     - Moderate evidence favouring the lower-AIC model.
   * - 4–10
     - Considerably less support for the higher-AIC model.
   * - > 10
     - The higher-AIC model is effectively ruled out.

In-process calibration with ``Calibrator``
------------------------------------------

The recommended way to calibrate MNiShed is **in-process**: drive the model
directly from Python with the Numba JIT warm, rather than launching a fresh
process per evaluation. This is roughly two orders of magnitude faster per
evaluation than a fork-based external optimizer, and one configuration then
serves both best-fit optimization and Bayesian uncertainty quantification.

The standard model setup is :class:`~mnished.Calibrator`, in which the *run
method is declared in config, not code*. A ``params.yml`` names the free
parameters and their bounds, and each parameter's ``target`` declares where its
value maps in the model — so one generic runner calibrates any basin with no
per-basin Python.

Declaring the parameter mapping
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each calibrated parameter carries a ``target`` naming its destination in
:func:`~mnished.run_and_score`:

.. code-block:: yaml

   driver:
     config_template: model.yml             # the model configuration
     metric:          KGE_logKGE_logFDC
     spin_up_cycles:  3
     routing_N:       2

   parameters:
     log__t_recession_shallow: {lower: 0.3, upper: 2.0, initial: 1.2, target: recession_coeff[0]}
     f_exfiltration_shallow:   {lower: 0.01, upper: 0.99, initial: 0.8, target: f_to_discharge[0]}
     PDD_melt_factor:          {lower: 0.1, upper: 10.0, initial: 1.0, target: melt_factor}

The ``target`` grammar:

* ``name`` — a scalar keyword (e.g. ``melt_factor``).
* ``name[i]`` — element *i* of a list keyword (e.g. ``recession_coeff[0]``);
  parameters sharing a keyword are grouped into one list.
* ``sub_catchments[I].key`` / ``sub_catchments[I].key[j]`` — an override on
  sub-catchment(s) *I*, where *I* is an index or comma-list (``0,1`` for a
  parameter shared across zones, e.g. two land zones with common physics, or a
  lake's ``f_route_lake``).
* A ``log__`` name prefix applies ``10**`` to the sampled value.

Untargeted list positions are left at their config value, so a parameter need
only name the elements it calibrates.

The Calibrator
^^^^^^^^^^^^^^

.. autoclass:: mnished.Calibrator
   :members: from_yaml, score, score_windows, run_kwargs
   :no-index:

``Calibrator`` builds the model **once** (via :class:`~mnished.ScoringModel`)
and reuses it every evaluation; :meth:`~mnished.Calibrator.score` is
bit-identical to the equivalent ``run_and_score`` call. It is sampler-agnostic —
point any optimizer or sampler at ``score`` and ``parameter_set``:

.. code-block:: python

   from mnished import Calibrator

   cal = Calibrator.from_yaml('params.yml')
   result = cal.score({'log__t_recession_shallow': 1.4,
                       'f_exfiltration_shallow': 0.6})   # or a vector ordered as cal.names
   result.score                                # the metric (result is a CalibResult)

Calibration windows
^^^^^^^^^^^^^^^^^^^^

The scoring span is set in one place — the driver's ``decades:`` key, a list of
``{start, end}`` windows (``None`` = full record):

.. code-block:: yaml

   driver:
     decades:
       - {start: '1991-01-01', end: '2000-12-31'}
       - {start: '2001-01-01', end: '2010-12-31'}

:meth:`~mnished.Calibrator.score_windows` scores a parameter set on each window
and returns one ``CalibResult`` per window; how to aggregate them — a mean score
for an optimizer, concatenated residuals for a likelihood — is the caller's
choice, so the Calibrator stays sampler-agnostic. The single-window
``decade_start`` / ``decade_end`` driver keys are a shorthand for a one-element
``decades:`` list, and :meth:`~mnished.Calibrator.score` scores that first (or
only) window; both spellings flow through the same mechanism. (``decades:`` is
named for its original decadal-backbone use; it generalizes to any windows and
will be renamed ``windows:`` in v4.0 — MNiMORPH/MNiShed#24.)

Build once, score many
^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: mnished.ScoringModel
   :no-index:

A plain ``run_and_score`` call re-reads the forcing and reconstructs the model
each time; for thousands of evaluations that rebuild dominates the wall-clock
(and a multi-window objective re-reads the forcing once per window per eval).
``ScoringModel`` — which ``Calibrator`` uses internally — does the build once.

Running a calibration
^^^^^^^^^^^^^^^^^^^^^^

The ``examples/cannon_inverse/calibrate.py`` runner wires ``Calibrator`` to
`SPOTPY <https://spotpy.readthedocs.io>`_, exposing both an optimizer (SCE-UA,
for best-fit) and a Bayesian sampler (DREAM, for the posterior):

.. code-block:: bash

   python calibrate.py sceua 5000          # best-fit
   python calibrate.py dream 20000 ar1     # posterior (AR(1) log-flow likelihood)

Run with the Numba JIT active (``pip install 'mnished[jit]'``; the JIT requires
``numpy < 2.3``), and serially: each evaluation returns a long simulation
vector, so multiprocessing loses more to inter-process transfer than it saves
on compute.

After a fit, characterise *how well the data constrained each parameter* with
the ``mnished.identifiability`` tools (per-parameter profiles and a curvature
eigenspectrum that names stiff vs. degenerate parameter combinations), or run
DREAM for the full posterior. ``log_flow_residual_terms`` exposes the scored
log-flow residuals used by the Bayesian likelihood.

For an **external** optimizer instead (e.g. Dakota, for expensive models or
cluster runs), write a thin driver that calls ``run_and_score`` directly; see
``examples/cannon_inverse/`` for both paths.

Practical Calibration Workflow
-------------------------------

Scoring window
^^^^^^^^^^^^^^

By default, ``run_and_score`` scores over all time steps that have observed
discharge.  Use ``start`` and ``end`` (YYYY-MM-DD strings or ``datetime``
objects) to restrict the scoring window:

.. code-block:: python

    result = run_and_score(
        'config.yml',
        recession_coeff=[671, 22, 50881],
        f_to_discharge=[0.553, 0.446, 1.0],
        recession_exponents=[4.62, 2.203, 1.0],
        start='1993-01-01',   # exclude first two years of spin-up uncertainty
        end='2011-09-30',
    )

Leaving a 1–2 year buffer at the start of the record is good practice when the
deep reservoir is slow (:math:`\tau` ≫ record length) and spin-up may not fully
equilibrate it.

Spin-up
^^^^^^^

``spin_up_cycles=None`` (the default) triggers automatic calculation:
``ceil(τ_max / record_length)``.  This is adequate when initial conditions are
set to the analytical steady-state (the default behaviour in ``run_and_score``).
Set ``spin_up_cycles=0`` when chaining decade runs — the ``final_states`` from
the previous run supply well-initialised starting conditions, so additional
spin-up is not needed.

Nash-cascade routing
^^^^^^^^^^^^^^^^^^^^

The model outputs point-source discharge from each reservoir.  To simulate
channel travel time (appropriate for large basins), ``routing_N`` and
``routing_K`` apply Nash-cascade attenuation after the reservoir cascade:

.. code-block:: python

    result = run_and_score(
        'config.yml',
        ...,
        routing_N=2,     # number of routing reservoirs (integer ≥ 1)
        routing_K=0.5,   # routing timescale (days)
    )

For small catchments (<500 km²), routing adds negligible improvement and
should be omitted.  For large basins, fit one routing step (N=1–3) with a
timescale of order (basin_length / wave_speed).  Routing parameters are
calibration parameters and add to :math:`k` if non-trivial.

Calibrating parallel sub-catchments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When the config partitions the basin into parallel sub-catchments (see
:ref:`sub-catchments-config`), pass a ``sub_catchments`` argument — one dict per
sub-catchment, in config order — to override that zone's ``area_fraction`` and
per-reservoir parameters by position.  The config owns the structure (how many
sub-catchments and reservoirs); the argument overrides the values each
evaluation:

.. code-block:: python

    result = run_and_score(
        'wild_rice.yml',
        sub_catchments=[
            {'area_fraction': 0.55, 'recession_coeff': [50, 500]},   # till uplands
            {'area_fraction': 0.45, 'recession_coeff': [1500]},      # clay lowlands
        ],
        melt_factor=4.0,        # snow and ET parameters stay basin-level
    )

This is mutually exclusive with the flat per-reservoir arguments
(``recession_coeff``, ``multipath_threshold``, …), which apply only to a
single-cascade config.  Each overridden value counts as one free parameter for
the AIC, plus ``n_sub_catchments − 1`` when the area fractions are calibrated.
Because ``final_states`` carries each zone's storage separately, decade chaining
(below) works unchanged for a partitioned basin.

Calibrating a lake
^^^^^^^^^^^^^^^^^^

A lake sub-catchment (see :ref:`lake-config`) calibrates through the same
``sub_catchments`` override.  Its outlet is the lake reservoir's
``recession_coeff`` and ``H_threshold``, so the two free lake parameters map as

* outlet coefficient :math:`a` → ``recession_coeff`` :math:`= 1/a`
  (small :math:`a` = slow outlet = large ``recession_coeff``; set bounds in
  ``recession_coeff`` space accordingly), and
* sill :math:`H_\text{sill}` → ``H_threshold``.

.. code-block:: python

    result = run_and_score(
        'crow_wing.yml',
        et_scale=0.6,                  # basin-level; keep free for lake-rich basins
        sub_catchments=[
            {'recession_coeff': [100],          # land cascade
             'multipath_threshold': [50.0], 'multipath_timescale': [5.0]},
            {'recession_coeff': [1 / 0.05],     # lake: a = 0.05
             'H_threshold': [200.0]},           # lake: H_sill = 200 mm
        ],
    )

The outlet exponent :math:`b` stays fixed at its config value (``5/3`` by
default) — it is not exposed to the override.  The bidirectional groundwater
exchange :math:`Q_\text{gw}` adds **no** calibrated parameter: it reuses the
partner land reservoir's own ``recession_coeff`` / ``recession_exponent``, so it
tracks the substrate you are already fitting.  ``f_route_lake`` is **not**
calibrated — set it from data in the YAML (lake position in the drainage
network).  Open-water evaporation likewise has no separate knob; it reuses the
basin ``et_scale``.  So a lake adds at most two free parameters (:math:`a`,
:math:`H_\text{sill}`) to the calibration.

Chaining decade runs
^^^^^^^^^^^^^^^^^^^^

To simulate non-stationary behaviour (changing land use, long-term climate
trends), chain runs across decades, passing end-state storage from one decade
as the initial state of the next:

.. code-block:: python

    result_d1 = run_and_score('config_1991_2000.yml', ..., spin_up_cycles=1)
    result_d2 = run_and_score(
        'config_2001_2010.yml',
        ...,
        initial_states=result_d1.final_states,
        spin_up_cycles=0,   # no spin-up — initial states already equilibrated
    )

Treat ``final_states`` as an opaque token to pass straight back as
``initial_states``: it carries the reservoir depths, snowpack, frozen-ground
index, and carried deficit at the end of the run.  For a partitioned basin it
is nested per sub-catchment, so each zone's snowpack and storage chain
independently.

Calibrating initial storage after spin-up
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``initial_states`` sets the reservoir storage *before* spin-up, which then
re-equilibrates it.  When a decade's initial condition is itself poorly
constrained — for example, when little forcing is available ahead of the
scored window to drive a meaningful spin-up — it can instead be treated as a
calibration target.  ``post_spinup_states`` injects reservoir depths *after*
spin-up completes and *before* the scored run begins, overriding the spin-up
end state:

.. code-block:: python

    result = run_and_score(
        'config.yml',
        ...,
        start='2001-01-01', end='2010-12-31',                  # decade mode
        post_spinup_states={'reservoirs': [None, None, H0_deep]},  # H0_deep from optimizer
        post_spinup_k=1,                                       # one free initial-storage parameter
    )

Only the ``'reservoirs'`` list is used, and any entry left as ``None`` keeps
that reservoir at its spin-up end state — here, only the deep reservoir's
initial storage is calibrated.  ``post_spinup_states`` is applied only in
decade mode (when ``start`` is set) and is ignored for a full-record run.
Set ``post_spinup_k`` to the number of injected depths the optimizer varies,
so they are counted as degrees of freedom in the AIC.

Suggested Parameter Sets
-------------------------

The following configurations are illustrative starting points derived from
calibration studies on temperate agricultural catchments in the upper Midwest
(USA).  They should be treated as starting points for a calibration, not as
universal defaults — the correct parameter values depend on the specific basin.

Three-reservoir temperate agricultural basin
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This structure (soil → intermediate → deep) is appropriate for glaciated
lowland catchments with tile drainage, shallow Quaternary aquifers, and
regional bedrock groundwater.  It is the recommended starting structure for
AIC-based exploration.

.. code-block:: yaml

    reservoirs:
        recession_coefficients:
            - 671      # soil: calibrate; log10 bounds [1.0, 4.5]
            - 22       # intermediate: calibrate; log10 bounds [1.0, 4.5]
            - 50000    # deep: calibrate; log10 bounds [3.5, 5.0]
        exfiltration_fractions:
            - 0.55     # soil → stream vs. percolating to intermediate
            - 0.45     # intermediate → stream vs. recharging deep
            - 1.0      # deep → all to stream
        maximum_effective_depths__mm:
            - .inf
            - .inf
            - .inf
        recession_exponents:
            - 4.62     # soil: calibrate [1.5, 6.0]; encodes tile-drain nonlinearity
            - 2.203    # intermediate: fix at B–N lower-envelope value (basin-specific)
            - 1.0      # deep: fix at 1.0 (linear; confined aquifer)

    snowmelt:
        PDD_melt_factor: 5.1   # mm SWE °C⁻¹ day⁻¹; calibrate [0.1, 10.0]
        snow_insulation_k: 0.0
        fdd_threshold: .inf    # frozen ground disabled; b_soil absorbs the signal

    modules:
        snowpack:          true
        frozen_ground:     false   # b_soil > 1 absorbs spring-pulse signal
        rain_on_snow:      true
        direct_runoff:     false
        dtr_fgi_decay:     false
        et_water_stress:   false
        et_reservoir_draw: true

    general:
        enforce_water_balance: global   # single WB multiplier; no hidden per-year d.o.f.
        et_alpha: 1.0

**Key design choices and rationale:**

- ``frozen_ground: false`` — the calibrated soil nonlinearity (:math:`b \approx 3`–5
  for tile-drained basins) reproduces the same spring fast-flow signal as frozen
  infiltration blockage.  Activating both simultaneously leads to equifinality.
  Use frozen ground only when independent soil temperature or frost-depth data
  warrant it; see :doc:`model_description` for the equifinality argument.

- ``recession_b_intermediate`` fixed at the basin-specific Brutsaert–Nieber
  lower-envelope value — the B–N lower envelope represents the slowest observable
  recession, which corresponds to the intermediate Quaternary zone (MRT of days to
  weeks).  See :doc:`recession_analysis` for the derivation and
  :class:`~mnished.BrutsaertNieber` for how to compute the estimate.

- ``recession_b_deep = 1.0`` fixed — the deep reservoir (MRT of decades to
  centuries) produces discharge too small and too slow to constrain :math:`b`
  from streamflow.  The linear approximation (Darcy flow in a confined or
  semi-confined matrix) is physically appropriate.

- ``enforce_water_balance: global`` — a single ET multiplier applied to the
  full record.  Recommended for AIC-based model comparison because
  ``'water-year'`` mode adds one hidden degree of freedom per year.

- ``et_reservoir_draw: true`` — ET drawn from post-cascade soil storage, giving
  a temporal buffer equal to the soil MRT without any additional free parameter.

**Free parameters for this structure (k = 7):**

.. list-table::
   :widths: 35 20 35
   :header-rows: 1

   * - Parameter
     - Suggested bounds
     - Notes
   * - ``PDD_melt_factor``
     - [0.1, 10.0]
     - Sensitive to climate; 1–6 mm SWE °C⁻¹ d⁻¹ typical
   * - ``log10(τ_soil)``
     - [1.0, 4.5]
     - Unconstrained; MRT more informative than raw τ
   * - ``recession_b_soil``
     - [1.5, 6.0]
     - High b encodes tile drainage or surface nonlinearity
   * - ``f_exfiltration_soil``
     - [0.01, 0.99]
     - —
   * - ``log10(τ_intermediate)``
     - [1.0, 4.5]
     - Unconstrained; allows model to self-discover timescale ordering
   * - ``f_exfiltration_intermediate``
     - [0.01, 0.99]
     - —
   * - ``log10(τ_deep)``
     - [3.5, 5.0]
     - ~1000–100,000 days; well-constrained range for bedrock GW

.. _tile-drain-degeneracy:

Tile drain parameter degeneracy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When the explicit tile-drain module (``f_tile``, ``tau_tile``) is used instead of
letting :math:`b_{\text{soil}} > 1` absorb tile-drain signals implicitly, a
steady-state parameter degeneracy limits interpretation.

**At steady state**, the tile-drain flux is independent of ``tau_tile``:

.. math::

   Q_{\text{tile}} = f_{\text{tile}} \cdot Q_{\text{to\_next}}

where :math:`Q_{\text{to\_next}}` is the downward percolation from the soil
reservoir.  The tile sub-reservoir reaches steady state with ``tau_tile``
controlling only the *lag* (and thus the shape of the recession limb), not
the long-term mean flux.

**Consequence for calibration**: a water-balance calibration that targets
decadal-mean fluxes cannot distinguish ``f_tile`` from the product
:math:`f_{\text{tile}} \cdot \tau_{\text{tile}}`.  Specifically:

.. math::

   f_{\text{tile}} \cdot \tau_{\text{tile}} \;\propto\; A_t

where :math:`A_t` is the areal fraction of the watershed covered by tile drains
(Hooghoudt 1940).  The calibrated ``f_tile`` alone encodes both areal extent
*and* the inverse of drain spacing squared — it is not a clean estimate of
tiled area unless ``tau_tile`` is independently constrained from site data
(e.g. drain spacing surveys, water-table monitoring, or hydrograph recession
analysis targeting the tile-drain recession component).

**Time-varying interpretation**: agricultural tile drainage in the Upper Midwest
expanded and intensified progressively through the 20th century.  Two physical
processes drive this:

- **Areal expansion** (new tiling): increases ``f_tile`` — more of the watershed
  routes water through drains.
- **Intensification** (closer drain spacing): decreases ``tau_tile`` ∝ drain
  spacing² (Hooghoudt 1940), accelerating the tile-drain response without
  necessarily changing ``f_tile``.

A per-decade transient calibration that holds ``tau_tile`` fixed and varies only
``f_tile`` therefore conflates these two processes: a trend in calibrated
``f_tile`` may reflect new tiling, closer spacing, or both.  Interpreting the
trend physically requires additional constraint on drain spacing or tile
installation records.

**Practical guidance**: unless detailed tile-drain records are available,
treat calibrated ``f_tile`` as a bulk "effective tile-drain efficiency" rather
than an areal coverage fraction.  Document ``tau_tile`` as a fixed structural
assumption when reporting results.

**Breaking the degeneracy — a path to physical interpretability**: if
``tau_tile`` can be independently constrained (from tile-drain mapping,
water-table monitoring, or recession decomposition), the degeneracy collapses
and ``f_tile`` recovers its intended physical meaning — the areal fraction of
the catchment that drains through tiles.  A calibration run that fixes
``tau_tile`` at a physically motivated value and calibrates only ``f_tile``
per decade would then produce a direct, comparable time series of effective
tile-drain coverage, with trend slopes interpretable as rates of agricultural
drainage expansion.  This is a scientifically meaningful improvement over the
implicit nonlinear-:math:`b` proxy approach, provided the ``tau_tile``
constraint is defensible.

Note that the parameter names ``f_tile`` and ``tau_tile`` in the model code
retain the same names regardless of whether they are calibrated or externally
constrained.  When ``tau_tile`` is fixed, ``f_tile`` gains the interpretation
of areal coverage; when both are free, they are individually unidentifiable
from water-balance data alone.  The documentation (or the params.yml
``description`` field) should state which regime applies for a given run.

Two-reservoir baseline
^^^^^^^^^^^^^^^^^^^^^^

A simpler structure for data-sparse basins or as a structural null hypothesis
before adding a third reservoir.  AIC comparison against the three-reservoir
structure indicates whether the additional reservoir is warranted.

.. code-block:: yaml

    reservoirs:
        recession_coefficients:
            - 30       # soil/shallow: calibrate
            - 5000     # groundwater: calibrate
        exfiltration_fractions:
            - 0.6      # calibrate
            - 1.0
        maximum_effective_depths__mm:
            - .inf
            - .inf
        recession_exponents:
            - 2.0      # soil: starting point from B–N; calibrate
            - 1.0      # groundwater: fix linear

    general:
        enforce_water_balance: global

This structure has k = 5 free parameters (τ × 2, f × 1, b_soil × 1,
PDD_melt_factor × 1).  If the three-reservoir structure gives ΔAIC < 2 compared
to this baseline, the additional reservoir is not justified by the data.
