Calibration
===========

HydroRaVENS calibration is built around :func:`~hydroravens.run_and_score`,
which runs the model with a given parameter set, applies optional Nash-cascade
routing and a scoring window, and returns a :class:`~hydroravens.CalibResult`
containing all goodness-of-fit metrics.  Integration with an external optimizer
(e.g. `Dakota <https://dakota.sandia.gov>`_) is achieved by writing a thin
driver that reads parameter values from the optimizer's output file and calls
``run_and_score`` with those values.

.. contents:: On this page
   :local:
   :depth: 2

run_and_score()
---------------

.. autofunction:: hydroravens.run_and_score

CalibResult
-----------

.. autoclass:: hydroravens.CalibResult
   :no-inherited-members:

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
     - Reservoir water depths (mm) at the end of the run.  Pass as
       ``initial_states=`` in the next decade run to chain simulations.
   * - ``buckets``
     - :class:`~hydroravens.Buckets`
     - The fully initialised and run model object.  Access
       ``result.buckets.hydrodata`` for the full simulated time series.

Goodness-of-Fit Metrics
------------------------

The ``metric`` argument to :func:`~hydroravens.run_and_score` selects the
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
that is passed as a non-trivial argument (e.g. each element of ``t_efold`` and
``f_to_discharge`` that is free).  Fixed structural choices (e.g. ``b_deep = 1``
by convention, not optimization) do not add to :math:`k`.

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
        t_efold=[671, 22, 50881],
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

The ``final_states`` dict has keys matching ``reservoir_order`` in the driver
config (e.g. ``{'soil': H_soil, 'intermediate': H_int, 'deep': H_deep}``).

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
        e_folding_residence_times__days:
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
  :class:`~hydroravens.BrutsaertNieber` for how to compute the estimate.

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

Two-reservoir baseline
^^^^^^^^^^^^^^^^^^^^^^

A simpler structure for data-sparse basins or as a structural null hypothesis
before adding a third reservoir.  AIC comparison against the three-reservoir
structure indicates whether the additional reservoir is warranted.

.. code-block:: yaml

    reservoirs:
        e_folding_residence_times__days:
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
