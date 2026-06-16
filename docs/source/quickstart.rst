Quick Start
===========

This guide will get you running your first MNiShed model.

Prepare Your Data
~~~~~~~~~~~~~~~~~

MNiShed requires daily time series data in CSV format. Required and optional columns:

.. list-table:: Input CSV Columns
   :widths: 35 15 40
   :header-rows: 1

   * - Column Name
     - Units
     - When required
   * - ``Date``
     - YYYY-MM-DD
     - Always; must be continuous daily data with no gaps
   * - ``Precipitation [mm/day]``
     - mm/day
     - Always
   * - ``Discharge [m^3/s]``
     - m³/s
     - Always; used for scoring
   * - ``Mean Temperature [C]``
     - °C
     - Snowpack and frozen-ground modules
   * - ``Minimum Temperature [C]``
     - °C
     - DTR-based frozen-ground decay; ``ThorntwaiteChang2019`` ET
   * - ``Maximum Temperature [C]``
     - °C
     - DTR-based frozen-ground decay; ``ThorntwaiteChang2019`` ET
   * - ``Photoperiod [hr]``
     - hours
     - ``ThorntwaiteChang2019`` ET method
   * - ``Evapotranspiration [mm/day]``
     - mm/day
     - ``datafile`` ET method

Example input (first few rows):

.. code-block:: text

    Date,Precipitation [mm/day],Discharge [m^3/s],Mean Temperature [C],Minimum Temperature [C],Maximum Temperature [C],Photoperiod [hr]
    2010-01-01,0.0,15.2,-2.5,-6.1,1.2,9.1
    2010-01-02,2.1,16.8,-1.3,-4.8,2.3,9.2
    2010-01-03,0.5,15.9,0.2,-3.1,3.6,9.3

Create a Configuration File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MNiShed is configured through a YAML file (``config.yml``).  The example
below uses a three-reservoir structure (soil → intermediate → deep), which
is a recommended starting point for temperate catchments.  See
:doc:`configuration` for all available options and :doc:`calibration` for
guidance on choosing a model structure.

.. code-block:: yaml

    timeseries:
        datafile: input_data.csv

    initial_conditions:
        water_reservoir_effective_depths__mm:
            - 50       # Soil reservoir
            - 500      # Intermediate (shallow GW / outwash)
            - 10000    # Deep groundwater
        snowpack__mm_SWE: 0

    catchment:
        drainage_basin_area__km2: 3800
        evapotranspiration_method: ThorntwaiteChang2019
        water_year_start_month: 10

    general:
        spin_up_cycles: null          # auto: ceil(τ_max / record_length)
        enforce_water_balance: global

    reservoirs:
        recession_timescales__days:
            - 200      # Soil (placeholder; calibrate)
            - 30       # Intermediate (placeholder; calibrate)
            - 50000    # Deep GW (placeholder; calibrate)
        exfiltration_fractions:
            - 0.5      # Soil: fraction to stream vs. percolating deeper
            - 0.5      # Intermediate: fraction to stream vs. recharging deep
            - 1.0      # Deep: all exits as baseflow
        maximum_effective_depths__mm:
            - .inf
            - .inf
            - .inf
        recession_exponents:
            - 3.0      # Soil: nonlinear (calibrate; typical 2–6)
            - 2.2      # Intermediate: fix near B–N lower-envelope value
            - 1.0      # Deep: linear (fix)

    snowmelt:
        PDD_melt_factor: 2.0          # mm SWE °C⁻¹ day⁻¹ (calibrate)
        snow_insulation_k: 0.0
        fdd_threshold: .inf           # disabled; set a finite value to enable

    modules:
        snowpack:          true
        frozen_ground:     false
        rain_on_snow:      true
        direct_runoff:     false
        dtr_fgi_decay:     false
        et_water_stress:   false
        et_reservoir_draw: true

Run the Model
~~~~~~~~~~~~~

**Using the Python API:**

.. code-block:: python

    import mnished

    model = mnished.Buckets()
    model.initialize('config.yml')
    model.run()
    model.compute_NSE(verbose=True)
    model.plot()

**Using the command-line interface:**

.. code-block:: bash

    mnished -y config.yml

Adjust Parameters
~~~~~~~~~~~~~~~~~

Model performance depends on the reservoir parameters:

**Residence times** (``recession_timescales__days``)
  Larger values = slower response. Order of magnitude ranges:

  * Soil zone: days to weeks (fast lateral drainage)
  * Intermediate (shallow GW / outwash): weeks to months
  * Deep groundwater: years to centuries

  Use :func:`~mnished.suggest_priors` to estimate data-driven starting
  points from the observed discharge record before manual adjustment.

**Exfiltration fractions** (``exfiltration_fractions``)
  Fraction of each reservoir's drainage going directly to the river.

  * Higher = more direct runoff from that reservoir
  * Lower = more infiltration to the next-deeper reservoir
  * Bottom reservoir must be 1.0 for mass conservation

**Recession exponents** (``recession_exponents``)
  Controls the nonlinearity of storage–discharge in each reservoir.
  :math:`b = 1` is a standard linear reservoir; :math:`b > 1` gives a concave
  (faster-draining-when-full) response.  Use :class:`~mnished.BrutsaertNieber`
  to estimate the appropriate exponent from observed streamflow recession.

  * Soil zone: calibrate (typical 2–6; higher values for tile-drained basins)
  * Intermediate: fix at Brutsaert–Nieber lower-envelope estimate (~2.0–2.5)
  * Deep: fix at 1.0 (linear)

**Water balance closure** (``enforce_water_balance``)
  ``'global'`` uses a single ET scaling factor over the full record and is
  recommended when comparing model structures by AIC.  ``'water-year'`` fits
  a separate multiplier each year, which can overfit interannual ET variability.

Next Steps
~~~~~~~~~~

* Use :func:`~mnished.suggest_priors` for data-driven parameter starting points
* Read :doc:`model_description` for the theory behind each component
* Read :doc:`calibration` for guidance on metric selection, AIC, and parameter sets
* Explore :doc:`configuration` for all configuration options
* See :doc:`recession_analysis` for estimating recession exponents from data
