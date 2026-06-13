Tutorial: From Data to Calibrated Model
========================================

This tutorial walks through a complete HydroRaVENS workflow, from raw
daily discharge and precipitation data to a calibrated model with
physically informed parameter starting points.

.. contents:: On this page
   :local:
   :depth: 2

Overview
--------

The recommended workflow has four stages:

1. **Estimate priors** — analyse the discharge record to get data-driven
   starting points for timescales, recession exponents, and initial
   storage depths (:func:`~hydroravens.suggest_priors`).
2. **Write a configuration file** — populate a YAML file using the
   suggested priors.
3. **Calibrate** — run the optimizer to refine parameters.
4. **Evaluate** — inspect goodness-of-fit metrics and the best-fit
   hydrograph.

Input Data Format
-----------------

HydroRaVENS reads a CSV file with at minimum these columns:

.. code-block:: text

    Date,Precipitation [mm/day],Discharge [m^3/s]

Optional columns activate additional model processes:

.. code-block:: text

    Evapotranspiration [mm/day]  # measured ET; enables evapotranspiration_method: datafile
    Mean Temperature [C]         # activates snowpack and frozen-ground modules
    Minimum Temperature [C]      # activates DTR-based FGI decay and ThorntwaiteChang2019
    Maximum Temperature [C]      # activates DTR-based FGI decay and ThorntwaiteChang2019
    Photoperiod [hr]             # required for ThorntwaiteChang2019 ET method

All columns must be at a daily timestep.  Missing values in the
discharge column are allowed and are excluded from scoring.

Stage 1 — Estimate Priors
--------------------------

Before running or calibrating the model, use :func:`~hydroravens.suggest_priors`
to extract data-driven parameter starting points directly from the
observed discharge record.

.. code-block:: python

    import pandas as pd
    from hydroravens import suggest_priors

    df = pd.read_csv('data/input.csv', parse_dates=['Date'])
    Q  = df['Specific Discharge [mm/day]'].values   # or compute from m³/s
    P  = df['Precipitation [mm/day]'].values

    pr = suggest_priors(Q, P=P, n_reservoirs=3)
    pr.summary()

The output looks like::

    ============================================================
    HydroRaVENS data-driven priors
    ============================================================

    E-folding residence times (fastest → slowest):
      soil    : 28.4 days
      karst   : 847.3 days
      deep    : could not be estimated — use calibration default

    Power-law recession exponents (fastest → slowest):
      soil    : 2.841  (B–N data-driven estimate — calibrate)
      karst   : 2.203  (theoretical B–N 2.203 — consider fixing)
      deep    : 1.000  (linear (b=1) — deep reservoir default)

    Initial storage depths [mm] (fastest → slowest):
      soil    : 45.2 mm
      karst   : 312.8 mm
      deep    : 18.4 mm

    Calibration bounds (log10 days) for params.yml:
      log__t_efold_soil:   initial=1.453,  lower=0.953,  upper=1.953
      log__t_efold_karst:  initial=2.928,  lower=2.428,  upper=3.428

    Brutsaert–Nieber recession cloud:
      slope n   = 1.4817
      coeff a   = 0.0052
      → b_HR    = 2.841  (used as b_soil prior)
    ============================================================

You can also visualise the recession cloud and get a YAML snippet:

.. code-block:: python

    pr.bn.plot()                   # log-log recession cloud with fit
    print(pr.to_yaml_snippet())    # paste into your config file

The YAML snippet gives a complete ``reservoirs:`` and
``initial_conditions:`` block ready to paste into your configuration
file.

**What the priors mean:**

- **Timescales** come from the spectral decomposition of the hydrograph
  (:class:`~hydroravens.HydrographSeparation`). Use them as starting
  points for calibration, not fixed values — the optimizer will refine
  them.
- **b_soil** (the B–N estimate) is a data-driven starting point.
  Calibrated values for agricultural catchments are typically *b* ≈ 3–4.
  Do not fix it.
- **b_karst = 2.203** is the theoretical Brutsaert & Nieber (1977)
  long-time baseflow value. Consider fixing it rather than calibrating
  — the data rarely constrain it independently. See
  :doc:`recession_analysis` for the theory.
- **b_deep = 1.0** — the deep reservoir is effectively a slow linear
  store; the linear approximation is adequate and avoids adding a free
  parameter.

Stage 2 — Write a Configuration File
--------------------------------------

Create a YAML configuration file. Start from the snippet produced by
``pr.to_yaml_snippet()`` and fill in the remaining sections:

.. code-block:: yaml

    timeseries:
        datafile: data/input.csv

    initial_conditions:
        water_reservoir_effective_depths__mm:
            - 45.2   # soil
            - 312.8  # karst
            - 18.4   # deep
        snowpack__mm_SWE: 0

    catchment:
        drainage_basin_area__km2: 3800
        evapotranspiration_method: ThorntwaiteChang2019
        water_year_start_month: 10

    general:
        spin_up_cycles: null        # auto: ceil(τ_max / record_length)
        enforce_water_balance: global

    reservoirs:
        e_folding_residence_times__days:
            - 28.4   # soil
            - 847.3  # karst
            - 10000  # deep — placeholder; calibrate or fix
        exfiltration_fractions:
            - 0.6    # soil — placeholder; calibrate
            - 0.5    # karst — placeholder; calibrate
            - 1.0    # deep — all to stream
        maximum_effective_depths__mm:
            - .inf
            - .inf
            - .inf
        recession_exponents:
            - 2.841  # soil — from B–N prior; calibrate
            - 2.203  # karst — theoretical B–N; consider fixing
            - 1.0    # deep — linear

    snowmelt:
        PDD_melt_factor: 2.0        # mm SWE / °C / day; starting point
        fdd_threshold: .inf         # °C·day; .inf disables frozen ground
        snow_insulation_k: 0.0      # mm⁻¹ SWE; 0 disables insulation

    modules:
        snowpack:          true
        frozen_ground:     true
        rain_on_snow:      true
        direct_runoff:     false
        dtr_fgi_decay:     false
        et_water_stress:   false
        et_reservoir_draw: true

See :doc:`configuration` for a full reference of all available options.

Stage 3 — Quick Check Before Calibrating
-----------------------------------------

Run the model once with the prior starting points to verify the
configuration is valid and the fit is in the right ballpark:

.. code-block:: python

    from hydroravens import Buckets

    model = Buckets()
    model.initialize('config.yml')
    model.run()
    model.compute_NSE(verbose=True)
    model.plot()

If the hydrograph shape is qualitatively reasonable (seasonal timing is
correct, magnitudes are in the right range), proceed to calibration.

**Diagnosing a poor Stage 3 fit:**

.. list-table::
   :widths: 40 55
   :header-rows: 1

   * - Symptom
     - Likely cause and fix
   * - Modeled peak timing shifted by weeks
     - Melt factor too high or low; soil τ too short or long.
       Adjust ``PDD_melt_factor`` and ``log__t_efold_soil``.
   * - Baseflow too high in summer
     - Intermediate or deep τ too short, or f_exfiltration too high for
       those layers.  Extend timescales; reduce exfiltration fractions.
   * - Rising limb too slow, flashy peaks absent
     - Soil τ too long; try shorter timescale or higher recession exponent.
   * - Annual water balance far off (model wet or dry)
     - Check ``drainage_basin_area__km2`` and discharge units.
       Verify ``enforce_water_balance`` is set; check ET data for unit errors.
   * - NSE near 0 or negative
     - Check input data alignment — is precipitation spatially appropriate
       for this basin?  Are dates in the correct format?

Stage 4 — Calibrate
--------------------

HydroRaVENS uses `Dakota <https://dakota.sandia.gov>`_ for calibration.
The calibration workflow is driven by a ``params.yml`` file that
specifies which parameters to calibrate and their bounds.  Use the
bounds from ``pr.log_t_efold_bounds`` as the starting point for
``log__t_efold_*`` parameters:

.. code-block:: yaml

    parameters:

      log__t_efold_soil:
        lower:   0.953    # from pr.log_t_efold_bounds
        upper:   1.953
        initial: 1.453

      log__t_efold_karst:
        lower:   2.428
        upper:   3.428
        initial: 2.928

      recession_b_soil:
        lower:   1.5
        upper:   6.0
        initial: 2.841    # from B–N prior

After calibration, compare results using AIC to evaluate model structure.
The AIC is lower for better-fitting models and penalises additional free
parameters.  A ΔAIC > 2 relative to a simpler model is required before
the more complex structure is considered justified.  See :doc:`calibration`
for the full AIC guide and metric selection.

Stage 5 — Evaluate
-------------------

After calibration, inspect the best-fit parameters and metrics:

.. code-block:: python

    from hydroravens import run_and_score

    result = run_and_score(
        'config.yml',
        t_efold           = [28.4, 847.3, 10000],
        f_to_discharge    = [0.65, 0.48],
        recession_exponents = [3.4, 2.203, 1.0],
        enforce_water_balance = 'global',
    )
    print(f"Score = {result.score:.3f}")   # composite metric (KGE_logKGE here)
    print(f"AIC   = {result.aic:.1f}")
    print(f"BFI   = {result.bfi_mod:.3f}  (obs: {result.bfi_obs:.3f})")

Key metrics to examine:

- **Score (KGE_logKGE)** — overall and low-flow performance combined.
- **AIC** — penalises extra parameters; use to compare model structures.
- **BFI** — check that modeled baseflow index matches observations.
- **KGE_logFDC** — flow duration curve match; sensitive to reservoir
  partitioning.

**Reporting mean residence time (MRT):**

Raw e-folding times (τ) are not physically comparable across reservoirs
or experiments when the recession exponent *b* differs.  After obtaining
best-fit parameters, compute MRT for each reservoir using
:meth:`~hydroravens.Reservoir.mean_residence_time` with the long-term mean
flux attributed to that layer as ``Q_ref``:

.. code-block:: python

    buckets = result.buckets

    # Mean annual discharge attributed to each reservoir.
    # Approximate as: total mean Q × exfiltration-weighted fraction.
    # For a quick cross-experiment comparison, use mean total Q as Q_ref.
    import numpy as np
    Q_mean = buckets.hydrodata['Specific Discharge (modeled) [mm/day]'].mean()

    for i, res in enumerate(buckets.reservoirs):
        mrt = res.mean_residence_time(Q_ref=Q_mean)
        print(f"Reservoir {i}: τ={res.t_recession:.1f}d, b={res.recession_exponent:.3f}, "
              f"MRT={mrt:.1f}d  (Q_ref={Q_mean:.3f} mm/d)")

MRT should be reported alongside τ and *b* whenever comparing calibrated
models across experiments or basins.  τ alone is misleading when *b* > 1,
because the optimizer can trade τ against *b* while holding MRT nearly
constant.  See :ref:`Mean Residence Time <mean-residence-time>` in
:doc:`model_description` for the derivation and interpretation.

See Also
--------

- :doc:`configuration` — full YAML reference
- :doc:`recession_analysis` — Brutsaert & Nieber theory and API
- :doc:`calibration` — scoring functions, AIC, and metric details
- :func:`~hydroravens.suggest_priors` — full API reference
- :class:`~hydroravens.HydrographSeparation` — timescale decomposition
