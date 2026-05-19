CSDMS Basic Model Interface
============================

HydroRaVENS includes a `CSDMS Basic Model Interface (BMI)
<https://bmi.csdms.io/>`_ wrapper that enables it to be driven by any
BMI-compliant coupling framework and to exchange variables with other
BMI models.

.. contents:: On this page
   :local:
   :depth: 2

Overview
--------

The BMI wrapper exposes HydroRaVENS as a scalar (lumped) model with a
single grid of rank 0 and size 1.  All variables are scalars representing
catchment-integrated quantities.

The wrapper supports two usage modes:

**File-driven** (standard workflow)
    The YAML configuration file points to a CSV containing all forcing
    data.  The framework calls :meth:`~hydroravens.BmiHydroRaVENS.update`
    repeatedly to step through the record.  No :meth:`~hydroravens.BmiHydroRaVENS.set_value`
    calls are needed.

**Online coupled**
    An upstream model provides forcing each timestep via
    :meth:`~hydroravens.BmiHydroRaVENS.set_value` before calling
    :meth:`~hydroravens.BmiHydroRaVENS.update`.  The CSV file still
    provides the initial time series (used for spin-up and as a default
    if a variable is not overridden).

Installation
------------

The BMI wrapper requires `bmipy <https://github.com/csdms/bmi-python>`_.
Install it with the ``bmi`` optional-dependency group:

.. code-block:: bash

    pip install 'hydroRaVENS[bmi]'

Usage
-----

File-driven
~~~~~~~~~~~

.. code-block:: python

    from hydroravens import BmiHydroRaVENS

    bmi = BmiHydroRaVENS()
    bmi.initialize("config.yml")

    while bmi.get_current_time() < bmi.get_end_time():
        bmi.update()

    bmi.finalize()

Online coupled
~~~~~~~~~~~~~~

.. code-block:: python

    import numpy as np
    from hydroravens import BmiHydroRaVENS

    bmi = BmiHydroRaVENS()
    bmi.initialize("config.yml")        # CSV values loaded for spin-up

    while bmi.get_current_time() < bmi.get_end_time():
        # Override forcing from an upstream model
        bmi.set_value(
            "atmosphere_water__liquid_equivalent_precipitation_rate",
            np.array([p_from_upstream])
        )
        bmi.set_value("atmosphere__temperature", np.array([t_from_upstream]))
        bmi.update()

        # Pass discharge downstream
        q = np.empty(1, dtype=np.float64)
        bmi.get_value("land_surface_water__runoff_volume_flux", q)
        downstream_model.set_value("channel_entrance__discharge", q)

    bmi.finalize()

Converting specific discharge to volumetric flow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``land_surface_water__runoff_volume_flux`` is area-normalised specific
discharge in mm d⁻¹.  To convert to volumetric discharge Q [m³ s⁻¹]:

.. code-block:: python

    area_km2 = bmi._model.drainage_basin_area__km2
    Q_m3s = q_mm_d * area_km2 * 1e3 / 86400

Exposed Variables
-----------------

All variables are scalar (grid rank 0, size 1, location ``node``).
Types are ``float64``; time unit is ``d`` (days).

Input variables
~~~~~~~~~~~~~~~

These variables are read from the CSV by default.  In online-coupled
mode, call :meth:`~hydroravens.BmiHydroRaVENS.set_value` before each
:meth:`~hydroravens.BmiHydroRaVENS.update` to override them.

Temperature and ET inputs are declared even when those columns are absent
from the CSV.  Calling :meth:`~hydroravens.BmiHydroRaVENS.set_value` for
an absent column raises :exc:`KeyError`.

.. list-table::
   :widths: 55 15 30
   :header-rows: 1

   * - CSDMS Standard Name
     - Units
     - HydroRaVENS column
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
:meth:`~hydroravens.BmiHydroRaVENS.update` and retrieved via
:meth:`~hydroravens.BmiHydroRaVENS.get_value`.

.. list-table::
   :widths: 55 15 30
   :header-rows: 1

   * - CSDMS Standard Name
     - Units
     - Source
   * - ``land_surface_water__runoff_volume_flux``
     - mm d⁻¹
     - Modelled specific discharge
   * - ``snowpack__liquid_equivalent_depth``
     - mm
     - Snowpack SWE; 0.0 if no snowpack
   * - ``subsurface_water__depth``
     - mm
     - Total subsurface storage (all reservoirs)
   * - ``subsurface_water_reservoir_0__depth``
     - mm
     - Reservoir 0 storage (shallowest)
   * - ``subsurface_water_reservoir_1__depth``
     - mm
     - Reservoir 1 storage
   * - ``subsurface_water_reservoir_2__depth``
     - mm
     - Reservoir 2 storage (deepest); ``nan`` if fewer than 3 reservoirs

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
     - 1.0 day (HydroRaVENS is a daily model)
   * - Start time
     - 0.0
   * - End time
     - Number of days in the input record

Shape, spacing, origin, coordinate, and connectivity methods raise
:exc:`NotImplementedError` — these are not defined for rank-0 scalar
grids.

Caveats
-------

**get_value_ptr is a snapshot, not a live pointer.** HydroRaVENS stores
scalar state as Python floats rather than numpy arrays.
:meth:`~hydroravens.BmiHydroRaVENS.get_value_ptr` therefore returns a
fresh length-1 array each call; the array does not update automatically
when the model advances.  Call
:meth:`~hydroravens.BmiHydroRaVENS.get_value` after each
:meth:`~hydroravens.BmiHydroRaVENS.update` to retrieve current values.

**Spin-up is internal.** Spin-up cycles specified in the YAML config run
inside :meth:`~hydroravens.BmiHydroRaVENS.initialize`; the BMI time
counter starts at 0.0 after spin-up completes.

**finalize() does not plot.** Calling
:meth:`~hydroravens.BmiHydroRaVENS.finalize` discards the internal model
object but does not call :meth:`~hydroravens.Buckets.finalize` on it,
which would trigger an NSE print and a plot pop-up unsuitable for
headless coupling runs.

API Reference
-------------

.. autoclass:: hydroravens.BmiHydroRaVENS
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
