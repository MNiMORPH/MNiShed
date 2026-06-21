CSDMS Basic Model Interface
============================

MNiShed includes a `CSDMS Basic Model Interface (BMI)
<https://bmi.csdms.io/>`_ wrapper that enables it to be driven by any
BMI-compliant coupling framework and to exchange variables with other
BMI models.

.. contents:: On this page
   :local:
   :depth: 2

Overview
--------

The BMI wrapper exposes MNiShed as a scalar (lumped) model with a
single grid of rank 0 and size 1.  All variables are scalars representing
catchment-integrated quantities.

The wrapper supports two usage modes:

**File-driven** (standard workflow)
    The YAML configuration file points to a CSV containing all forcing
    data.  The framework calls :meth:`~mnished.BmiMNiShed.update`
    repeatedly to step through the record.  No :meth:`~mnished.BmiMNiShed.set_value`
    calls are needed.

**Online coupled**
    An upstream model provides forcing each timestep via
    :meth:`~mnished.BmiMNiShed.set_value` before calling
    :meth:`~mnished.BmiMNiShed.update`.  The CSV file still
    provides the initial time series (used for spin-up and as a default
    if a variable is not overridden).

Installation
------------

The BMI wrapper requires `bmipy <https://github.com/csdms/bmi-python>`_.
Install it with the ``bmi`` optional-dependency group:

.. code-block:: bash

    pip install 'MNiShed[bmi]'

Usage
-----

File-driven
~~~~~~~~~~~

.. code-block:: python

    from mnished import BmiMNiShed

    bmi = BmiMNiShed()
    bmi.initialize("config.yml")

    while bmi.get_current_time() < bmi.get_end_time():
        bmi.update()

    bmi.finalize()

Use :meth:`~mnished.BmiMNiShed.update_until` to advance to a
specific time without writing the loop yourself::

    bmi.update_until(365.0)   # advance one year

Online coupled
~~~~~~~~~~~~~~

.. code-block:: python

    import numpy as np
    from mnished import BmiMNiShed

    bmi = BmiMNiShed()
    bmi.initialize("config.yml")        # CSV values loaded for spin-up

    while bmi.get_current_time() < bmi.get_end_time():
        # Override forcing from an upstream model
        bmi.set_value(
            "atmosphere_water__liquid_equivalent_precipitation_rate",
            np.array([p_from_upstream])
        )
        bmi.set_value("atmosphere__temperature", np.array([t_from_upstream]))
        bmi.update()

        # Pass volumetric discharge [m³ s⁻¹] to a downstream channel or
        # sediment model.  Specific discharge [mm d⁻¹] is also available
        # via "land_surface_water__runoff_volume_flux".
        q_m3s = np.empty(1, dtype=np.float64)
        bmi.get_value("channel_exit_water__volume_flow_rate", q_m3s)
        downstream_model.set_value("channel_exit_water__volume_flow_rate",
                                   q_m3s)

    bmi.finalize()

Converting specific discharge to volumetric flow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``land_surface_water__runoff_volume_flux`` is area-normalised specific
discharge in mm d⁻¹.  To convert to volumetric discharge Q [m³ s⁻¹]:

.. code-block:: python

    area_km2 = bmi._model.drainage_basin_area__km2
    Q_m3s = q_mm_d * 1e-3 * area_km2 * 1e6 / 86400  # mm→m, km²→m², day→s

Exposed Variables
-----------------

All variables are scalar (grid rank 0, size 1, location ``node``).
Types are ``float64``; time unit is ``d`` (days).

Input variables
~~~~~~~~~~~~~~~

These variables are read from the CSV by default.  In online-coupled
mode, call :meth:`~mnished.BmiMNiShed.set_value` before each
:meth:`~mnished.BmiMNiShed.update` to override them.

Temperature and ET inputs are declared even when those columns are absent
from the CSV.  Calling :meth:`~mnished.BmiMNiShed.set_value` for
an absent column raises :exc:`KeyError`.

.. list-table::
   :widths: 55 15 30
   :header-rows: 1

   * - CSDMS Standard Name
     - Units
     - MNiShed column
   * - ``atmosphere_water__liquid_equivalent_precipitation_rate``
     - mm d⁻¹
     - ``Precipitation [mm/day]``
   * - ``atmosphere__temperature``
     - °C
     - ``Mean Temperature [C]``
   * - ``atmosphere__minimum_temperature``
     - °C
     - ``Minimum Temperature [C]``
   * - ``atmosphere__maximum_temperature``
     - °C
     - ``Maximum Temperature [C]``
   * - ``land_surface_water__potential_evapotranspiration_volume_flux``
     - mm d⁻¹
     - ``Evapotranspiration [mm/day]``

Output variables
~~~~~~~~~~~~~~~~

These variables are updated by each call to
:meth:`~mnished.BmiMNiShed.update` and retrieved via
:meth:`~mnished.BmiMNiShed.get_value`.

.. list-table::
   :widths: 55 15 30
   :header-rows: 1

   * - CSDMS Standard Name
     - Units
     - Source
   * - ``land_surface_water__runoff_volume_flux``
     - mm d⁻¹
     - Modelled specific discharge (area-normalised)
   * - ``channel_exit_water__volume_flow_rate``
     - m³ s⁻¹
     - Volumetric discharge (specific discharge × catchment area)
   * - ``snowpack__liquid_equivalent_depth``
     - mm
     - Snowpack SWE; 0.0 if no snowpack
   * - ``subsurface_water__depth``
     - mm
     - Total subsurface storage (all reservoirs)
   * - ``land_surface_water__evapotranspiration_volume_flux``
     - mm d⁻¹
     - Model evapotranspiration flux (after water-balance scaling)
   * - ``land_surface_water__direct_runoff_volume_flux``
     - mm d⁻¹
     - Hortonian-style fast bypass (``direct_runoff`` module)
   * - ``land_surface_water__baseflow_volume_flux``
     - mm d⁻¹
     - Constant regional baseflow (``baseflow_Q``); see note below
   * - ``land_surface_water__tile_drain_volume_flux``
     - mm d⁻¹
     - Tile-drain sub-reservoir discharge (``tile_fractions``)
   * - ``land_surface_water__multipath_drain_volume_flux``
     - mm d⁻¹
     - Threshold-activated parallel drain (``multipath_thresholds__mm``)
   * - ``land_surface__frozen_ground_index``
     - degC d
     - Frozen-ground index (FGI) state
   * - ``subsurface_water_reservoir_0__depth``
     - mm
     - Reservoir 0 storage (shallowest)
   * - ``subsurface_water_reservoir_1__depth`` … ``subsurface_water_reservoir_9__depth``
     - mm
     - Reservoirs 1–9 storage (deeper); ``nan`` for indices ≥ number of
       configured reservoirs

.. note::

   **Flux partition.** The four ``*_volume_flux`` components above
   (direct runoff, baseflow, tile drain, multipath drain) decompose the
   fast-flow contributions to discharge and are recorded by each
   :meth:`~mnished.BmiMNiShed.update`.  They are diagnostic: the
   primary cascade discharge is reported by
   ``land_surface_water__runoff_volume_flux``.  ``baseflow`` is the
   constant regional-import term (``baseflow_Q``); like ``run_and_score``,
   it is applied as an output-layer addition and is **not** part of the
   routed cascade, so it is exposed separately and is *not* folded into
   ``land_surface_water__runoff_volume_flux``.  A coupled model that wants
   total discharge including regional import should add the two.

   The CSDMS Standard Names for evapotranspiration, direct runoff, and
   baseflow extend the registered ``land_surface_water__…_volume_flux``
   family; ``tile_drain``, ``multipath_drain``, and ``frozen_ground_index``
   are MNiShed-specific quantities named to follow the same convention.

.. note::

   MNiShed itself places no limit on the number of reservoirs — you
   can add as many as you like to the ``reservoirs:`` block in the YAML
   configuration.  The BMI wrapper caps *exposed* reservoir outputs at 10
   (indices 0–9) because the BMI specification requires variable names to
   be fixed at import time.  If you configure more than 10 reservoirs,
   :meth:`~mnished.BmiMNiShed.initialize` will raise a
   :exc:`ValueError` with instructions pointing to the four constants in
   ``mnished/bmi.py`` that need updating:
   ``_OUTPUT_VAR_NAMES``, ``_VAR_UNITS``, ``_RESERVOIR_DEPTH_NAMES``, and
   ``_BMI_MAX_RESERVOIRS``.  The total subsurface storage across all
   reservoirs is always available via ``subsurface_water__depth``
   regardless of reservoir count.

Grid and Time
-------------

.. list-table::
   :widths: 40 60
   :header-rows: 0

   * - Grid ID
     - 0 (all variables share one scalar grid)
   * - Grid rank
     - 0 (scalar — no spatial dimensions)
   * - Grid size
     - 1
   * - Grid type
     - ``"scalar"``
   * - Time unit
     - ``"d"`` (days)
   * - Timestep
     - 1.0 day (MNiShed is a daily model)
   * - Start time
     - 0.0
   * - End time
     - Number of days in the input record

Shape, spacing, origin, coordinate, and connectivity methods raise
:exc:`NotImplementedError` — these are not defined for rank-0 scalar
grids.

Caveats
-------

**get_value_ptr is a snapshot, not a live pointer.** MNiShed stores
scalar state as Python floats rather than numpy arrays.
:meth:`~mnished.BmiMNiShed.get_value_ptr` therefore returns a
fresh length-1 array each call; the array does not update automatically
when the model advances.  Call
:meth:`~mnished.BmiMNiShed.get_value` after each
:meth:`~mnished.BmiMNiShed.update` to retrieve current values.

**Spin-up is internal.** Spin-up cycles specified in the YAML config run
inside :meth:`~mnished.BmiMNiShed.initialize`; the BMI time
counter starts at 0.0 after spin-up completes.

**finalize() does not plot.** Calling
:meth:`~mnished.BmiMNiShed.finalize` discards the internal model
object but does not call :meth:`~mnished.Buckets.finalize` on it,
which would trigger an NSE print and a plot pop-up unsuitable for
headless coupling runs.

**Calling update() past the end of the record raises an error.** The
file-driven loop ``while bmi.get_current_time() < bmi.get_end_time()``
terminates correctly.  Calling :meth:`~mnished.BmiMNiShed.update`
after all rows have been consumed will raise a ``KeyError`` from the
internal pandas DataFrame.  Guard against this in custom loops by checking
``get_current_time() < get_end_time()`` before each call.

API Reference
-------------

.. autoclass:: mnished.BmiMNiShed
   :members: initialize, update, update_until, finalize,
             get_component_name,
             get_input_item_count, get_output_item_count,
             get_input_var_names, get_output_var_names,
             get_var_grid, get_var_type, get_var_units,
             get_var_itemsize, get_var_nbytes, get_var_location,
             get_start_time, get_end_time, get_current_time,
             get_time_step, get_time_units,
             get_grid_rank, get_grid_size, get_grid_type,
             get_value, get_value_ptr, get_value_at_indices,
             set_value, set_value_at_indices
   :member-order: bysource
