Input Contract
==============

MNiShed consumes two inputs: a **forcing CSV** (a daily time series) and a
**config YAML**. The :mod:`mnished.io` module is the authoritative, code-grounded
specification of both — it *hosts the data spec* — and provides a fast pre-flight
validator that checks a config/forcing pair against it and reports every problem
at once.

This is a **contract**-level check: are the right columns and config sections
present, given the options you chose? It complements, and does not replace, the
deeper model-consistency validation that :meth:`mnished.Buckets.initialize`
performs (sub-catchment area fractions, lake parameters, junction types, mass
conservation). Run the validator first for a quick, complete list of input
problems; ``initialize`` then checks the model is internally consistent.

Validating inputs
-----------------

.. code-block:: python

    import mnished

    report = mnished.validate_inputs("crow_wing_config.yml")
    print(report)            # human-readable list of errors + warnings
    report.raise_if_errors()  # raise ValueError if anything is wrong

    if report.ok:
        b = mnished.Buckets()
        b.initialize("crow_wing_config.yml")

:func:`~mnished.validate_inputs` loads the config, resolves
``timeseries.datafile`` relative to it, loads the forcing, and runs both checks.
Use :func:`~mnished.validate_config` or :func:`~mnished.validate_forcing`
directly to check one in isolation (``validate_forcing`` takes an optional
``config`` so it knows which conditional columns to require).

Errors vs. warnings
~~~~~~~~~~~~~~~~~~~~~

* An **error** is something that will make ``initialize`` fail or give wrong
  results — a missing required column, an unknown ET method, a malformed cascade.
* A **warning** is an input MNiShed accepts but with degraded behaviour — e.g. a
  missing ``Mean Temperature [C]`` silently disables snowpack; missing
  min/max temperature falls back to a constant frozen-ground-index decay.

The forcing-CSV spec
--------------------

The programmatic source of truth is :data:`mnished.io.FORCING_COLUMNS`.

.. list-table::
   :header-rows: 1
   :widths: 30 12 58

   * - Column
     - Units
     - Required
   * - ``Date``
     - ISO date
     - always (continuous daily series, exactly 1-day intervals)
   * - ``Precipitation [mm/day]``
     - mm/day
     - always
   * - ``Discharge [m^3/s]``
     - m³/s
     - always (converted to mm/day via ``drainage_basin_area__km2``)
   * - ``Mean Temperature [C]``
     - °C
     - the daily mean drives snowmelt and frozen ground — **error** if an active
       ``snowmelt.fdd_threshold``, **warning** if ``modules.snowpack`` — but
       supply *either* this column *or* both Min and Max (the model derives the
       mean from their midpoint). Thornthwaite ET uses Min/Max directly, not this.
   * - ``Minimum Temperature [C]``
     - °C
     - **error** for ``ThornthwaiteChang2019`` (the Chang effective temperature
       uses the diurnal range); with Max also supplies the snowmelt mean and the
       ``dtr_fgi_decay`` range
   * - ``Maximum Temperature [C]``
     - °C
     - **error** for ``ThornthwaiteChang2019`` (see Minimum Temperature); with
       Min also supplies the snowmelt mean and the ``dtr_fgi_decay`` range
   * - ``Photoperiod [hr]``
     - hours
     - **error** if ``evapotranspiration_method`` is ``ThornthwaiteChang2019``
   * - ``Evapotranspiration [mm/day]``
     - mm/day
     - **error** if ``evapotranspiration_method`` is ``datafile``

The config-YAML spec
--------------------

The programmatic source of truth is :data:`mnished.io.CONFIG_SECTIONS`.

.. list-table::
   :header-rows: 1
   :widths: 22 18 60

   * - Section
     - Required
     - Required keys
   * - ``timeseries``
     - yes
     - ``datafile``
   * - ``catchment``
     - yes
     - ``drainage_basin_area__km2``, ``evapotranspiration_method``
       (``datafile`` | ``ThornthwaiteChang2019``), ``water_year_start_month``
   * - ``general``
     - yes
     - ``spin_up_cycles``
   * - ``reservoirs``
     - one of ``reservoirs`` / ``sub_catchments``
     - ``recession_coefficients``, ``exfiltration_fractions``,
       ``maximum_effective_depths__mm`` (equal-length lists)
   * - ``initial_conditions``
     - with a top-level ``reservoirs``
     - ``water_reservoir_effective_depths__mm`` (length = number of reservoirs)
   * - ``sub_catchments``
     - one of ``reservoirs`` / ``sub_catchments``
     - a non-empty list of land/lake zones
   * - ``snowmelt``, ``modules``, ``phenology``
     - no
     - — (see :doc:`configuration` for their options and defaults)

API
---

.. autofunction:: mnished.validate_inputs

.. autofunction:: mnished.validate_config

.. autofunction:: mnished.validate_forcing

.. autoclass:: mnished.ValidationReport
   :members: ok, raise_if_errors, extend
   :member-order: bysource
