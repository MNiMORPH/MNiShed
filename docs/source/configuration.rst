Configuration Reference
========================

MNiShed is configured entirely through a YAML file. This page documents all 
available options.

Configuration File Structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

    timeseries:
        # Data input section
    
    initial_conditions:
        # Starting state of reservoirs and snowpack
    
    catchment:
        # Basin properties
    
    general:
        # Simulation settings
    
    reservoirs:
        # Reservoir cascade configuration
    
    snowmelt:
        # Required section; snowpack processes are only active if
        # Mean Temperature [C] is present in the input CSV

The ``timeseries`` section
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 15 50
   :header-rows: 1

   * - Option
     - Type
     - Description
   * - ``datafile``
     - string
     - Path to CSV input file with daily time series

Example:

.. code-block:: yaml

    timeseries:
        datafile: data/streamflow_2010_2020.csv

The ``initial_conditions`` section
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 35 15 40
   :header-rows: 1

   * - Option
     - Type
     - Description
   * - ``water_reservoir_effective_depths__mm``
     - list of floats
     - Initial water depth (mm) in each reservoir, listed top to bottom
   * - ``snowpack__mm_SWE``
     - float
     - Initial snowpack depth in mm snow-water equivalent (SWE)

Example:

.. code-block:: yaml

    initial_conditions:
        water_reservoir_effective_depths__mm:
            - 150    # Top (soil) reservoir starts at 150 mm
            - 400    # Bottom (groundwater) starts at 400 mm
        snowpack__mm_SWE: 50  # 50 mm SWE initially

The ``catchment`` section
~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 35 15 40
   :header-rows: 1

   * - Option
     - Type
     - Description
   * - ``drainage_basin_area__km2``
     - float
     - Basin area (km²); used for discharge-to-depth conversion
   * - ``evapotranspiration_method``
     - string
     - ``datafile`` or ``ThornthwaiteChang2019``. See note below.
   * - ``water_year_start_month``
     - int
     - Month (1–12) when water year begins (10 = October for USGS)
   * - ``baseflow_Q``
     - float
     - Constant regional groundwater import (mm/day). Added to modeled
       discharge after routing; not mass-balanced against P or ET.
       Default ``0.0`` (disabled). Use only when independent
       hydrogeological evidence supports an external groundwater source.

.. note::

    ``ThornthwaiteChang2019`` requires long-term monthly temperature
    normals to compute the Thornthwaite heat index :math:`I` and
    exponent :math:`a`. If normals are not supplied via
    ``T_monthly_normals`` in the :class:`~mnished.Buckets`
    constructor, they are computed automatically from the mean monthly
    temperatures in the input record. A ``UserWarning`` is raised when
    the record is shorter than 20 years, as a short period may not
    represent long-term climatology. For best results, supply normals
    from a 30-year reference period (e.g. WMO 1991–2020 normals).

    Water years with no discharge observations receive an ET multiplier
    of 1.0 (raw ET, no water-balance correction) with a ``UserWarning``
    naming the affected years.

Example:

.. code-block:: yaml

    catchment:
        drainage_basin_area__km2: 3800
        evapotranspiration_method: datafile
        water_year_start_month: 10
        baseflow_Q: 0.0   # disabled; set > 0 for regional groundwater import

The ``general`` section
~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 25 15 50
   :header-rows: 1

   * - Option
     - Type
     - Description
   * - ``spin_up_cycles``
     - int or ``null``
     - Number of complete passes through data before the main run. Use 0 to
       skip spin-up (e.g. when supplying ``initial_states`` for chained
       decade runs). ``null`` (or omitting the key when calling
       :func:`~mnished.calibration.run_and_score`) triggers automatic
       calculation: ``ceil(τ_max / record_length)``, where ``τ_max`` is the
       longest reservoir e-folding time. Because initial conditions are set
       to analytical steady-state depths, one e-folding time is sufficient
       to resolve seasonal and inter-annual climate memory.
   * - ``et_alpha``
     - float
     - Fraction (0–1) of potential ET drawn from the top (soil) reservoir
       when ``et_reservoir_draw: true``. The remainder ``1 − et_alpha`` is
       drawn from the second reservoir. Default ``1.0`` (all ET from top
       reservoir). Has no effect when ``et_reservoir_draw: false``. See
       :ref:`et-modules`.
   * - ``direct_runoff_fraction``
     - float
     - Fast-bypass fraction (0–1) of positive daily recharge that exits
       directly as runoff, bypassing the reservoir cascade. Active only
       when ``direct_runoff: true`` in the ``modules`` block. Default
       ``0.0`` (disabled). Also settable as a calibration parameter via
       :func:`~mnished.calibration.run_and_score`.
   * - ``enforce_water_balance``
     - string
     - Controls how ET is scaled to close the water balance.
       Default ``'water-year'`` if the key is absent. Accepted values:

       ``'water-year'`` — scale ET by a per-water-year multiplier so that
       P − Q − ET = 0 over each water year.

       ``'global'`` — scale ET by a single multiplier computed from the full
       record (sum(P − Q) / sum(ET_raw)). No per-year overfitting; does not
       add hidden degrees of freedom to AIC comparisons. Recommended when
       calibrating with a long record and comparing models by AIC.

       ``'none'`` — use raw ET without any correction. Appropriate only when
       supplying trusted measured ET (e.g. eddy covariance). Using ``'none'``
       with ``ThornthwaiteChang2019`` raises a warning because Thornthwaite ET
       carries large systematic biases. Also appropriate when ``et_scale`` is
       a free calibration parameter (see :ref:`et-modules`).

       Legacy boolean values are accepted: ``true`` maps to ``'water-year'``
       and ``false`` maps to ``'none'``.

Example:

.. code-block:: yaml

    general:
        spin_up_cycles: 2                  # Run through data twice to initialize
        enforce_water_balance: water-year  # default; omit to accept the default
        et_alpha: 1.0                      # fraction of ET from top reservoir (et_reservoir_draw only)
        direct_runoff_fraction: 0.0        # fast-bypass fraction (direct_runoff module only)

The ``reservoirs`` section
~~~~~~~~~~~~~~~~~~~~~~~~~~

This section defines the cascade of reservoirs (1 or more); each may be
linear or nonlinear (power-law, via ``recession_exponents``).
All lists must have the same length.

.. list-table::
   :widths: 35 15 40
   :header-rows: 1

   * - Option
     - Type
     - Description
   * - ``recession_coefficients``
     - list of floats
     - Recession coefficient per reservoir. For a **linear** reservoir
       (``recession_exponents`` entry = 1) this is the e-folding drainage
       timescale in days; for a **nonlinear** reservoir (exponent > 1) it is a
       drainage coefficient (units day·mm^(b−1)), **not** a timescale — use
       :meth:`~mnished.Reservoir.mean_residence_time` for a comparable
       timescale (see :ref:`mean residence time <mean-residence-time>`).
   * - ``exfiltration_fractions``
     - list of floats
     - Fraction (0–1) of drainage exiting as discharge. Used by the
       ``fraction`` and ``threshold`` junction types; ignored for
       ``leakance``.
   * - ``junction_types``
     - list of strings
     - Routing rule at each reservoir's outlet to the next-deeper reservoir:
       ``fraction`` (default), ``leakance``, or ``threshold``. See
       :ref:`reservoir-junctions`. Default: all ``fraction`` (the
       ``exfiltration_fractions`` split).
   * - ``leakance_R__days``
     - list of floats or ``null``
     - Leakance resistance :math:`R` (days) for reservoirs whose
       ``junction_types`` entry is ``leakance``. Downward flow is
       :math:`\max(H_i - H_{i+1}, 0)/R`, capped at the total drainage;
       larger :math:`R` impedes it. Required for a ``leakance`` junction;
       ignored otherwise. Default: all ``null``.
   * - ``H_threshold__mm``
     - list of floats
     - Dead-storage depth (mm) for reservoirs whose ``junction_types`` entry
       is ``threshold``: only :math:`\max(H - H_\text{thr}, 0)` drains, so
       storage below the threshold is retained indefinitely. Default: all
       ``0.0`` (no threshold; full storage drains).
   * - ``maximum_effective_depths__mm``
     - list of floats
     - Storage capacity (mm) per reservoir; use ``.inf`` for unlimited
   * - ``pdm_H0__mm``
     - list of floats or ``null``
     - PDM characteristic depth (mm) per reservoir. When set for a
       reservoir, a fraction :math:`f_\text{sat} = 1 - e^{-H / H_0}` of
       incoming recharge becomes saturation-excess runoff. ``null`` (or
       omitting the entry) disables PDM for that reservoir. Default: all
       ``null`` (PDM off). Mutually exclusive with a finite
       ``maximum_effective_depths__mm`` entry for the same reservoir.
   * - ``tile_fractions``
     - list of floats
     - Fraction (0–1) of infiltrating water diverted to a fast tile-drain
       sub-reservoir, per main reservoir. Default: all ``0.0`` (no tile
       drainage). When non-zero, the corresponding
       ``tile_residence_times__days`` entry must also be set.
   * - ``tile_residence_times__days``
     - list of floats or ``null``
     - E-folding residence time (days) of the tile-drain sub-reservoir per
       main reservoir. Required when the corresponding ``tile_fractions``
       entry is greater than zero; ignored otherwise. Default: all ``null``.
   * - ``multipath_thresholds__mm``
     - list of floats or ``null``
     - Storage depth (mm) above which a *parallel* fast drainage path
       activates for each reservoir. ``null`` disables multipath for that
       reservoir. When set, adds
       :math:`Q_\text{mp} = \max(0, H - H_\text{thr})/\tau_\text{mp}` to
       discharge alongside the primary recession. Requires the matching
       ``multipath_timescales__days`` entry. Distinct from the
       ``tile_fractions``/``tile_residence_times__days`` mechanism (see
       below). Default: all ``null`` (multipath off).
   * - ``multipath_timescales__days``
     - list of floats or ``null``
     - E-folding timescale (days) of the parallel multipath drain.
       Required when the corresponding ``multipath_thresholds__mm`` entry
       is set; ignored otherwise. Default: all ``null``.
   * - ``recession_exponents``
     - list of floats
     - Power-law recession exponent :math:`b` per reservoir. :math:`b = 1`
       recovers the standard linear reservoir; :math:`b > 1` produces
       concave recession limbs consistent with subsurface flow theory
       (Brutsaert & Nieber 1977; Kirchner 2009). The theoretical value
       for catchment-integrated baseflow recession is :math:`b \approx 2.2`
       (Brutsaert & Nieber 1977); calibrated soil-zone values are typically
       larger (:math:`b \approx 3`–5). Default: all ``1.0`` (linear).
       Also overridable per calibration run via
       :func:`~mnished.calibration.run_and_score`. See
       :doc:`model_description` for theory and :doc:`recession_analysis`
       for estimating *b* from observed streamflow.

Example:

.. code-block:: yaml

    reservoirs:
        recession_coefficients:
            - 16      # Interflow: fast lateral drainage, days to weeks
            - 200     # Soil zone: seasonal storage, months
            - 3650    # Groundwater: multi-year storage
        exfiltration_fractions:
            - 0.8     # 80% to stream, 20% percolates to soil zone
            - 0.1     # 10% to stream, 90% recharges groundwater
            - 1.0     # All exits as baseflow
        maximum_effective_depths__mm:
            - .inf
            - .inf
            - .inf

No reservoir is fixed to a particular process; physical meaning is set
by the parameters. Successive reservoirs naturally span progressively
longer storage timescales — from interflow (days) to soil moisture
(months) to groundwater (years) — but that mapping is the user's choice.

The ``snowmelt`` section
~~~~~~~~~~~~~~~~~~~~~~~~

This section is **always required** in the configuration file. Snowpack
*processes* are only activated when ``Mean Temperature [C]`` is present in
the input CSV, but the ``PDD_melt_factor`` key must be present regardless.

.. list-table::
   :widths: 25 15 50
   :header-rows: 1

   * - Option
     - Type
     - Description
   * - ``PDD_melt_factor``
     - float
     - Melt rate (mm SWE per °C per day)
   * - ``fgi_decay_coeff``
     - float
     - Maximum daily FGI decay rate :math:`A_0` (default ``0.97``,
       Molnau & Bissell 1983). When ``Minimum Temperature [C]`` and
       ``Maximum Temperature [C]`` are present in the input CSV, the
       effective per-day coefficient :math:`A_t` is computed from the
       diurnal temperature range: :math:`A_t = 1 - (1-A_0)\,f_{\text{above}}`,
       where :math:`f_{\text{above}}` is the fraction of the day above
       0°C (linear diurnal cycle). On all-below-zero days :math:`A_t = 1`
       (no decay); on days that straddle 0°C, :math:`A_t` decreases toward
       :math:`A_0`. This makes FGI accumulation persistent in continental
       climates and decay-bounded in maritime ones. Without T_min/T_max,
       falls back to constant :math:`A_0` (original M&B behaviour).
   * - ``fdd_threshold``
     - float
     - Frozen ground index (FGI) threshold (°C·day) above which the top
       reservoir's lateral drainage is blocked and all water exits as
       direct runoff (Molnau & Bissell 1983). The FGI accumulates freezing
       degree-days and decays during thaw; when it exceeds this threshold
       the soil is considered frozen. Default ``inf`` (frozen ground never
       triggers). Requires ``frozen_ground: true`` in the ``modules``
       block and temperature data. Typically calibrated; literature values
       for agricultural soils range from ~50–200 °C·day. Also overridable
       via :func:`~mnished.calibration.run_and_score`. **Do not
       calibrate simultaneously with** ``snow_insulation_k`` — see
       :doc:`model_description`.
   * - ``snow_insulation_k``
     - float
     - Snow insulation decay constant (mm⁻¹ SWE). Scales effective air
       temperature reaching the soil as
       :math:`T_\text{eff} = T \cdot e^{-k \cdot \text{SWE}}`, reducing
       FGI accumulation under deep snowpack. Default ``0.0`` (no
       insulation). The exponential form originates in Molnau & Bissell
       (1983); LISFLOOD uses ~0.057 mm⁻¹ (van der Knijff et al. 2010)
       as a literature starting point. **Do not calibrate simultaneously
       with** ``fdd_threshold`` — the two parameters are correlated and
       will produce equifinal solutions. Fix one from independent data
       before calibrating the other; see :doc:`model_description`.

Example:

.. code-block:: yaml

    snowmelt:
        PDD_melt_factor: 1.0     # 1 mm SWE melts per °C per day
        fdd_threshold: 83.0      # °C·day; calibrate or fix from soil data
        snow_insulation_k: 0.057 # LISFLOOD default; set 0.0 to disable

The ``modules`` section
~~~~~~~~~~~~~~~~~~~~~~~

Optional process modules can be enabled or disabled through the ``modules``
block. All keys default to the values shown below if the block is absent.

.. list-table::
   :widths: 25 10 10 45
   :header-rows: 1

   * - Key
     - Type
     - Default
     - Description
   * - ``snowpack``
     - bool
     - ``true``
     - Accumulation, melt (PDD), and rain-on-snow. Requires
       ``Mean Temperature [C]`` in the input CSV; silently inactive
       if the column is absent regardless of this flag.
   * - ``frozen_ground``
     - bool
     - ``true``
     - Frozen ground index (FGI; Molnau & Bissell 1983). Accumulates
       freezing degree-days; blocks shallow-to-soil infiltration when
       the index exceeds ``fdd_threshold``. Requires ``snowpack: true``
       and temperature data.
   * - ``rain_on_snow``
     - bool
     - ``true``
     - Sensible-heat contribution of warm rain falling on snowpack.
       Only active when ``snowpack: true``.
   * - ``direct_runoff``
     - bool
     - ``false``
     - Hortonian-inspired fast-bypass fraction. A fixed fraction
       (``f_direct_runoff``) of positive daily recharge bypasses the
       reservoir cascade entirely. Conceptually motivated by
       infiltration-excess overland flow, but not a rigorous physical
       representation at the daily timestep.
   * - ``dtr_fgi_decay``
     - bool
     - ``true``
     - DTR-based FGI decay coefficient. When enabled and
       ``Minimum Temperature [C]`` / ``Maximum Temperature [C]``
       columns are present, the per-day FGI decay coefficient
       :math:`A_t` varies with the diurnal temperature range rather
       than being held constant at ``fgi_decay_coeff``. Disable to
       revert to constant :math:`A_0` (original Molnau & Bissell
       1983 behaviour). Has no effect when T_min/T_max columns are
       absent or when ``frozen_ground: false``.
   * - ``et_water_stress``
     - bool
     - ``false``
     - Soil-moisture-dependent ET. Actual ET is multiplied by
       :math:`1 - e^{-H_0 / H_{\text{pdm}}}`, where :math:`H_0` is
       the current storage in the top reservoir and
       :math:`H_{\text{pdm}}` is its PDM characteristic depth
       (``pdm_H0__mm``). ET approaches potential ET as the reservoir
       fills and drops to zero as it empties. Mutually exclusive with
       ``et_reservoir_draw``; if both are set, ``et_reservoir_draw``
       takes precedence.
   * - ``et_reservoir_draw``
     - bool
     - ``false``
     - Post-cascade reservoir ET extraction. Instead of subtracting
       ET from precipitation before it enters the reservoir cascade,
       ET is drawn directly from reservoir storage after drainage.
       Partitioned between the top two reservoirs by ``et_alpha``
       (see :ref:`et-modules`). Supports a wilting-point threshold
       (``wp_soil``, ``wp_soil_sigma``) below which soil-reservoir ET
       extraction is reduced or blocked. Mutually exclusive with
       ``et_water_stress``.

Example:

.. code-block:: yaml

    modules:
        snowpack:          true
        frozen_ground:     true
        rain_on_snow:      true
        direct_runoff:     false   # off by default
        dtr_fgi_decay:     true    # on by default; disable for continental climates with A=1
        et_water_stress:   false   # off by default; enable for soil-moisture-limited ET
        et_reservoir_draw: false   # off by default; enable for post-cascade reservoir ET

When a module is disabled, its associated parameter (e.g.
``PDD_melt_factor`` when ``snowpack: false``) has no effect and need not
be calibrated.

.. _et-modules:

ET Module Parameters
~~~~~~~~~~~~~~~~~~~~

The parameters below are calibration parameters associated with the ET
modules. They are passed to :func:`~mnished.calibration.run_and_score`
rather than read directly from the YAML file (except ``et_alpha``, which
may also be set in the ``general`` block).

.. list-table::
   :widths: 25 15 50
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``et_scale``
     - float
     - Global ET multiplier. Scales potential ET up or down uniformly
       across the entire record. Active (calibrated) only when
       ``et_water_stress: true`` or ``et_reservoir_draw: true`` *and*
       ``enforce_water_balance: 'none'`` in the driver. Incompatible with
       ``'water-year'`` or ``'global'`` water-balance enforcement (those
       modes compute their own multiplier). Default ``1.0``.
   * - ``et_alpha``
     - float
     - Fraction (0–1) of potential ET extracted from the top (soil)
       reservoir; ``1 − et_alpha`` is extracted from the second reservoir.
       Active only when ``et_reservoir_draw: true``. Default ``1.0``.
       May also be set as ``et_alpha`` in the ``general`` YAML block.
   * - ``wp_soil``
     - float
     - Wilting-point threshold (mm water depth) for the soil reservoir
       when ``et_reservoir_draw: true``. ET extraction from the soil
       reservoir is reduced or blocked when storage falls below this
       value. Default ``0.0`` (no wilting point). When ``wp_soil_sigma``
       is also set, a smooth Gaussian-CDF transition replaces the hard
       threshold.
   * - ``wp_soil_sigma``
     - float
     - Standard deviation (mm) of the spatially distributed wilting-point
       threshold. When ``> 0``, the fraction of the catchment below
       wilting point is modeled as a Gaussian CDF centered on ``wp_soil``
       with standard deviation ``wp_soil_sigma``, and ET extraction is
       reduced proportionally. Default ``0.0`` (hard threshold). Requires
       ``wp_soil > 0``.

.. _sub-catchments-config:

Parallel Sub-catchments (Optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default a configuration describes a single vertical reservoir cascade for
the whole basin. A basin can instead be partitioned into **parallel
sub-catchments** — spatially distinct hydraulic zones (for example till
uplands versus lake-clay lowlands) that drain to the same channel in parallel
rather than in a vertical cascade. Each sub-catchment has its own reservoir
cascade, snowpack, and frozen-ground state; basin discharge is the
area-weighted mean of the sub-catchments. See
:ref:`parallel-sub-catchments` in the model description for the concept.

Replace the top-level ``reservoirs`` (and the reservoir part of
``initial_conditions``) with a ``sub_catchments`` list. Each entry mirrors the
single-cascade form: a unique ``name``, an ``area_fraction``, its own
``reservoirs`` block (same keys as above), and an optional
``initial_conditions`` block for that zone's starting reservoir depths and SWE.

.. code-block:: yaml

    sub_catchments:
      - name: till_uplands
        area_fraction: 0.55
        reservoirs:
          recession_coefficients:        [50, 500]   # a two-reservoir cascade
          exfiltration_fractions:        [0.6, 1.0]
          maximum_effective_depths__mm:  [.inf, .inf]
          multipath_thresholds__mm:      [100.0, null]
          multipath_timescales__days:    [10.0, null]
        initial_conditions:                            # optional; defaults to 0
          water_reservoir_effective_depths__mm: [5, 300]
          snowpack__mm_SWE: 0
      - name: clay_lowlands
        area_fraction: 0.45
        reservoirs:
          recession_coefficients:        [1500]       # a single reservoir
          exfiltration_fractions:        [1.0]
          maximum_effective_depths__mm:  [.inf]

Rules:

* ``area_fraction`` values must sum to 1 (within ``1e-6``).
* ``name`` is required and must be unique.
* Each sub-catchment must have at least one reservoir.
* The number of reservoirs may differ between sub-catchments.
* A ``forcing`` block per sub-catchment is reserved for a future release (for
  zone-specific precipitation/ET/temperature) and currently raises
  ``NotImplementedError``; forcing is shared across sub-catchments for now.

The legacy single ``reservoirs`` block is equivalent to one sub-catchment with
``area_fraction: 1.0`` and reproduces the original behaviour exactly, so
existing configurations need no changes.

To calibrate a partitioned basin, pass a ``sub_catchments`` argument to
:func:`~mnished.calibration.run_and_score` (one dict per sub-catchment, in
config order) to override ``area_fraction`` and the per-reservoir parameters by
position; snow and ET parameters remain basin-level.

Complete Example
~~~~~~~~~~~~~~~~

.. code-block:: yaml

    timeseries:
        datafile: data/input.csv
    
    initial_conditions:
        water_reservoir_effective_depths__mm: [2, 400]
        snowpack__mm_SWE: 0

    catchment:
        drainage_basin_area__km2: 3800
        evapotranspiration_method: ThornthwaiteChang2019  # normals auto-computed if not supplied
        water_year_start_month: 10
        baseflow_Q: 0.0   # mm/day regional groundwater import; 0 = disabled

    general:
        spin_up_cycles: null  # auto: ceil(tau_max / record_length)

    reservoirs:
        recession_coefficients: [16, 2000]
        exfiltration_fractions: [0.8, 1.0]
        maximum_effective_depths__mm: [.inf, .inf]

    snowmelt:
        PDD_melt_factor: 1.0
        fgi_decay_coeff: 0.97       # default; Molnau & Bissell (1983)
        fdd_threshold: .inf         # °C·day; .inf = frozen ground never triggers
        snow_insulation_k: 0.0      # mm⁻¹ SWE; 0 = disabled

    modules:
        snowpack:          true
        frozen_ground:     true
        rain_on_snow:      true
        direct_runoff:     false
        dtr_fgi_decay:     true
        et_water_stress:   false   # soil-moisture-dependent ET
        et_reservoir_draw: false   # post-cascade reservoir ET extraction

Illustrative Starting Points
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These are rough starting points only, not calibrated or validated parameter
sets. Actual values depend strongly on the specific catchment.

**Alpine watershed with snowpack:**

.. code-block:: yaml

    reservoirs:
        recession_coefficients: [8, 1500]
        exfiltration_fractions: [0.7, 1.0]
    
    snowmelt:
        PDD_melt_factor: 1.5
    
    catchment:
        evapotranspiration_method: ThornthwaiteChang2019

**Humid lowland catchment:**

.. code-block:: yaml

    reservoirs:
        recession_coefficients: [25, 3000]
        exfiltration_fractions: [0.9, 1.0]
    
    catchment:
        evapotranspiration_method: datafile

**Arid/semi-arid watershed:**

.. code-block:: yaml

    reservoirs:
        recession_coefficients: [3, 500]
        exfiltration_fractions: [0.5, 1.0]
    
    snowmelt:
        PDD_melt_factor: 0.5
