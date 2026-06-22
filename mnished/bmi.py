"""
mnished.bmi
~~~~~~~~~~~~~~~
CSDMS Basic Model Interface (BMI) wrapper for MNiShed.

Wraps the :class:`~mnished.Buckets` class so that MNiShed can be
driven by any CSDMS-compliant framework or coupled online with other BMI
models (e.g. a channel-routing or sediment-transport model receiving daily
streamflow from a gauged watershed).

References
----------
Peckham, S.D., Hutton, E.W.H., and Norris, B. (2013). A component-based
approach to integrated modeling in the geosciences: The design of CSDMS.
Computers & Geosciences, 53, 3-12. https://doi.org/10.1016/j.cageo.2012.04.002
"""

import numpy as np
import pandas as pd

try:
    from bmipy import Bmi
except ImportError as exc:
    raise ImportError(
        "bmipy is required for the BMI wrapper. "
        "Install it with:  pip install 'MNiShed[bmi]'"
    ) from exc

from .mnished import Buckets

# ---------------------------------------------------------------------------
# CSDMS Standard Names and variable metadata
# ---------------------------------------------------------------------------

_INPUT_VAR_NAMES = (
    "atmosphere_water__liquid_equivalent_precipitation_rate",
    "atmosphere__temperature",
    "atmosphere__minimum_temperature",
    "atmosphere__maximum_temperature",
    "land_surface_water__potential_evapotranspiration_volume_flux",
)

_OUTPUT_VAR_NAMES = (
    "land_surface_water__runoff_volume_flux",
    "channel_exit_water__volume_flow_rate",
    "snowpack__liquid_equivalent_depth",
    "subsurface_water__depth",
    "land_surface_water__evapotranspiration_volume_flux",
    "land_surface_water__direct_runoff_volume_flux",
    "land_surface_water__baseflow_volume_flux",
    "land_surface_water__tile_drain_volume_flux",
    "land_surface_water__multipath_drain_volume_flux",
    "land_surface__frozen_ground_index",
    "subsurface_water_reservoir_0__depth",
    "subsurface_water_reservoir_1__depth",
    "subsurface_water_reservoir_2__depth",
    "subsurface_water_reservoir_3__depth",
    "subsurface_water_reservoir_4__depth",
    "subsurface_water_reservoir_5__depth",
    "subsurface_water_reservoir_6__depth",
    "subsurface_water_reservoir_7__depth",
    "subsurface_water_reservoir_8__depth",
    "subsurface_water_reservoir_9__depth",
)

# Hard cap: reservoirs 0–9 are declared as BMI output variables.
# To support more than 10 reservoirs, add names to _OUTPUT_VAR_NAMES,
# entries to _VAR_UNITS and _RESERVOIR_DEPTH_NAMES, and raise this constant.
_BMI_MAX_RESERVOIRS = 10

_ALL_VAR_NAMES = frozenset(_INPUT_VAR_NAMES + _OUTPUT_VAR_NAMES)

_VAR_UNITS = {
    "atmosphere_water__liquid_equivalent_precipitation_rate": "mm d-1",
    "atmosphere__temperature":                               "degC",
    "atmosphere__minimum_temperature":                       "degC",
    "atmosphere__maximum_temperature":                       "degC",
    "land_surface_water__potential_evapotranspiration_volume_flux": "mm d-1",
    "land_surface_water__runoff_volume_flux":                "mm d-1",
    "channel_exit_water__volume_flow_rate":                  "m3 s-1",
    "snowpack__liquid_equivalent_depth":                     "mm",
    "subsurface_water__depth":                               "mm",
    "land_surface_water__evapotranspiration_volume_flux":    "mm d-1",
    "land_surface_water__direct_runoff_volume_flux":         "mm d-1",
    "land_surface_water__baseflow_volume_flux":              "mm d-1",
    "land_surface_water__tile_drain_volume_flux":            "mm d-1",
    "land_surface_water__multipath_drain_volume_flux":       "mm d-1",
    "land_surface__frozen_ground_index":                     "degC d",
    "subsurface_water_reservoir_0__depth":                   "mm",
    "subsurface_water_reservoir_1__depth":                   "mm",
    "subsurface_water_reservoir_2__depth":                   "mm",
    "subsurface_water_reservoir_3__depth":                   "mm",
    "subsurface_water_reservoir_4__depth":                   "mm",
    "subsurface_water_reservoir_5__depth":                   "mm",
    "subsurface_water_reservoir_6__depth":                   "mm",
    "subsurface_water_reservoir_7__depth":                   "mm",
    "subsurface_water_reservoir_8__depth":                   "mm",
    "subsurface_water_reservoir_9__depth":                   "mm",
}

# MNiShed DataFrame column name for each input variable
_INPUT_COLUMNS = {
    "atmosphere_water__liquid_equivalent_precipitation_rate": "Precipitation [mm/day]",
    "atmosphere__temperature":        "Mean Temperature [C]",
    "atmosphere__minimum_temperature": "Minimum Temperature [C]",
    "atmosphere__maximum_temperature": "Maximum Temperature [C]",
    "land_surface_water__potential_evapotranspiration_volume_flux":
        "Evapotranspiration [mm/day]",
}

# Ordered reservoir depth variable names (index = reservoir index in model)
_RESERVOIR_DEPTH_NAMES = (
    "subsurface_water_reservoir_0__depth",
    "subsurface_water_reservoir_1__depth",
    "subsurface_water_reservoir_2__depth",
    "subsurface_water_reservoir_3__depth",
    "subsurface_water_reservoir_4__depth",
    "subsurface_water_reservoir_5__depth",
    "subsurface_water_reservoir_6__depth",
    "subsurface_water_reservoir_7__depth",
    "subsurface_water_reservoir_8__depth",
    "subsurface_water_reservoir_9__depth",
)


def _check_var(name: str) -> None:
    if name not in _ALL_VAR_NAMES:
        raise KeyError(
            f"Unknown variable {name!r}. "
            f"Available: {sorted(_ALL_VAR_NAMES)}"
        )


# ---------------------------------------------------------------------------
# BMI class
# ---------------------------------------------------------------------------

class BmiMNiShed(Bmi):
    """
    CSDMS Basic Model Interface wrapper for MNiShed.

    Wraps a :class:`~mnished.Buckets` instance so that MNiShed
    can participate in a CSDMS-compliant coupling framework.

    Two usage modes are supported:

    **File-driven** (standard MNiShed workflow) — the YAML config
    points to a CSV containing all forcing data; the framework steps
    through the record by calling :meth:`update` repeatedly::

        from mnished import BmiMNiShed
        import numpy as np

        bmi = BmiMNiShed()
        bmi.initialize("config.yml")
        while bmi.get_current_time() < bmi.get_end_time():
            bmi.update()
        bmi.finalize()

    **Online coupled** — an upstream model provides forcing each step via
    :meth:`set_value` before calling :meth:`update`::

        bmi.initialize("config.yml")
        bmi.set_value(
            "atmosphere_water__liquid_equivalent_precipitation_rate",
            np.array([5.2])
        )
        bmi.set_value("atmosphere__temperature", np.array([3.1]))
        bmi.update()
        q = np.empty(1, dtype=np.float64)
        bmi.get_value("land_surface_water__runoff_volume_flux", q)
        print(f"Discharge: {q[0]:.3f} mm d-1")

    Notes
    -----
    **Specific discharge**: ``land_surface_water__runoff_volume_flux`` is
    area-normalised specific discharge in mm d⁻¹, not volumetric flux.
    To convert to volumetric discharge Q [m³ s⁻¹]:

    .. code-block:: python

        Q_m3s = q_mm_d * area_km2 * 1e3 / 86400

    where ``area_km2 = bmi._model.drainage_basin_area__km2``.

    **Difference from run_and_score discharge**: the streaming BMI reports
    the raw per-step reservoir-cascade discharge.  It does *not* apply the
    two output-layer post-processing steps that
    :func:`mnished.calibration.run_and_score` performs on the full series:
    Nash-cascade flow routing (``routing_K``), which is an inherently
    batch convolution unavailable to a per-step interface, and the constant
    regional baseflow term (``baseflow_Q``).  Baseflow is instead exposed as
    its own output (``land_surface_water__baseflow_volume_flux``) so a
    coupler can add it explicitly.  A configuration calibrated with routing
    and/or baseflow will therefore not reproduce its scored hydrograph
    through the BMI unless the coupler reapplies those terms.

    **get_value_ptr**: MNiShed stores scalar state as Python floats,
    not numpy arrays.  :meth:`get_value_ptr` therefore returns a fresh
    length-1 array rather than a live pointer into model memory.  Values
    in the returned array do not update when the model advances; call
    :meth:`get_value` after each :meth:`update` call to retrieve current
    values.

    **Per-reservoir depths**: ``subsurface_water_reservoir_0__depth``
    through ``subsurface_water_reservoir_9__depth`` correspond to
    reservoirs 0–9 (shallowest to deepest) in the configuration.
    Depths for reservoir indices that do not exist in the current
    configuration are returned as ``np.nan``.

    **Optional input variables**: the five input Standard Names are always
    declared.  Temperature and ET inputs will raise :exc:`KeyError` from
    :meth:`set_value` if the corresponding column is absent from the
    loaded time series.
    """

    def __init__(self) -> None:
        self._model: "Buckets | None" = None
        self._current_time: float = 0.0
        self._end_time: float = 0.0

    # -----------------------------------------------------------------------
    # Lifecycle
    # -----------------------------------------------------------------------

    def initialize(self, config_file: str) -> None:
        """
        Load a MNiShed YAML configuration file and prepare the model.

        Reads the configuration, loads the input CSV time series, builds
        the reservoir stack, and runs spin-up cycles.  After this call,
        :meth:`update` steps through the record one day at a time.

        Parameters
        ----------
        config_file : str
            Path to a MNiShed YAML configuration file.
        """
        self._model = Buckets()
        self._model.initialize(config_file)
        n = len(self._model.reservoirs)
        if n > _BMI_MAX_RESERVOIRS:
            raise ValueError(
                f"This model configuration has {n} reservoirs, but the BMI "
                f"wrapper declares outputs for at most {_BMI_MAX_RESERVOIRS} "
                f"(subsurface_water_reservoir_0__depth through "
                f"subsurface_water_reservoir_{_BMI_MAX_RESERVOIRS - 1}__depth). "
                f"To support more reservoirs, add names to _OUTPUT_VAR_NAMES, "
                f"_VAR_UNITS, and _RESERVOIR_DEPTH_NAMES in mnished/bmi.py "
                f"and raise _BMI_MAX_RESERVOIRS."
            )
        # Buckets.initialize() ends with _timestep_i at end-of-record if
        # spin_up_cycles > 0 (spin-up is run via run(), which exhausts the
        # index).  Reset to the start so BMI update() steps from row 0.
        self._model._timestep_i = self._model.hydrodata.index[0]
        self._current_time = 0.0
        self._end_time = float(len(self._model.hydrodata))

    def update(self) -> None:
        """
        Advance the model by one day.

        If :meth:`set_value` was called for any input variable before this
        step, those values override the CSV values for the current row.
        """
        self._model.update()
        self._current_time += self._model.dt

    def update_until(self, time: float) -> None:
        """
        Advance the model until ``get_current_time() >= time``.

        Parameters
        ----------
        time : float
            Target time [days since start of record].
        """
        while self._current_time < time:
            self.update()

    def finalize(self) -> None:
        """
        Release internal resources.

        Discards the :class:`~mnished.Buckets` instance.  Does not
        call :meth:`~mnished.Buckets.finalize` on the inner model
        (which would trigger an NSE print and a plot pop-up).
        """
        self._model = None

    # -----------------------------------------------------------------------
    # Component information
    # -----------------------------------------------------------------------

    def get_component_name(self) -> str:
        return "MNiShed"

    def get_input_item_count(self) -> int:
        return len(_INPUT_VAR_NAMES)

    def get_output_item_count(self) -> int:
        return len(_OUTPUT_VAR_NAMES)

    def get_input_var_names(self):
        return _INPUT_VAR_NAMES

    def get_output_var_names(self):
        return _OUTPUT_VAR_NAMES

    # -----------------------------------------------------------------------
    # Variable information
    # -----------------------------------------------------------------------

    def get_var_grid(self, name: str) -> int:
        _check_var(name)
        return 0

    def get_var_type(self, name: str) -> str:
        _check_var(name)
        return "float64"

    def get_var_units(self, name: str) -> str:
        _check_var(name)
        return _VAR_UNITS[name]

    def get_var_itemsize(self, name: str) -> int:
        _check_var(name)
        return 8  # float64 = 8 bytes

    def get_var_nbytes(self, name: str) -> int:
        _check_var(name)
        return 8  # scalar: 1 element × 8 bytes

    def get_var_location(self, name: str) -> str:
        _check_var(name)
        return "node"

    # -----------------------------------------------------------------------
    # Time
    # -----------------------------------------------------------------------

    def get_start_time(self) -> float:
        return 0.0

    def get_end_time(self) -> float:
        return self._end_time

    def get_current_time(self) -> float:
        return self._current_time

    def get_time_step(self) -> float:
        return 1.0

    def get_time_units(self) -> str:
        return "d"

    # -----------------------------------------------------------------------
    # Grid — scalar (rank 0, size 1, grid_id 0)
    # -----------------------------------------------------------------------

    def get_grid_rank(self, grid_id: int) -> int:
        if grid_id != 0:
            raise ValueError(f"Unknown grid_id {grid_id!r}; only grid 0 exists.")
        return 0

    def get_grid_size(self, grid_id: int) -> int:
        if grid_id != 0:
            raise ValueError(f"Unknown grid_id {grid_id!r}; only grid 0 exists.")
        return 1

    def get_grid_type(self, grid_id: int) -> str:
        if grid_id != 0:
            raise ValueError(f"Unknown grid_id {grid_id!r}; only grid 0 exists.")
        return "scalar"

    # Rank-0 grids carry no shape, spacing, origin, coordinate arrays,
    # or connectivity — raise NotImplementedError for all such methods.

    def get_grid_shape(self, grid_id: int, shape: np.ndarray) -> np.ndarray:
        raise NotImplementedError("Scalar grid (rank 0) has no shape array.")

    def get_grid_spacing(self, grid_id: int, spacing: np.ndarray) -> np.ndarray:
        raise NotImplementedError("Scalar grid (rank 0) has no spacing array.")

    def get_grid_origin(self, grid_id: int, origin: np.ndarray) -> np.ndarray:
        raise NotImplementedError("Scalar grid (rank 0) has no origin array.")

    def get_grid_x(self, grid_id: int, x: np.ndarray) -> np.ndarray:
        raise NotImplementedError("Scalar grid (rank 0) has no coordinate arrays.")

    def get_grid_y(self, grid_id: int, y: np.ndarray) -> np.ndarray:
        raise NotImplementedError("Scalar grid (rank 0) has no coordinate arrays.")

    def get_grid_z(self, grid_id: int, z: np.ndarray) -> np.ndarray:
        raise NotImplementedError("Scalar grid (rank 0) has no coordinate arrays.")

    def get_grid_node_count(self, grid_id: int) -> int:
        raise NotImplementedError("Scalar grid (rank 0) has no node count.")

    def get_grid_edge_count(self, grid_id: int) -> int:
        raise NotImplementedError("Scalar grid (rank 0) has no edge count.")

    def get_grid_face_count(self, grid_id: int) -> int:
        raise NotImplementedError("Scalar grid (rank 0) has no face count.")

    def get_grid_edge_nodes(
        self, grid_id: int, edge_nodes: np.ndarray
    ) -> np.ndarray:
        raise NotImplementedError("Scalar grid (rank 0) has no connectivity.")

    def get_grid_face_nodes(
        self, grid_id: int, face_nodes: np.ndarray
    ) -> np.ndarray:
        raise NotImplementedError("Scalar grid (rank 0) has no connectivity.")

    def get_grid_face_edges(
        self, grid_id: int, face_edges: np.ndarray
    ) -> np.ndarray:
        raise NotImplementedError("Scalar grid (rank 0) has no connectivity.")

    def get_grid_nodes_per_face(
        self, grid_id: int, nodes_per_face: np.ndarray
    ) -> np.ndarray:
        raise NotImplementedError("Scalar grid (rank 0) has no connectivity.")

    # -----------------------------------------------------------------------
    # Internal helpers
    # -----------------------------------------------------------------------

    def _scalar_output(self, name: str) -> float:
        """Return the scalar value of an output variable after the last step."""
        m = self._model
        idx = m._timestep_i - 1
        first = m.hydrodata.index[0]

        if name == "land_surface_water__runoff_volume_flux":
            if idx < first:
                return np.nan
            val = m.hydrodata.at[idx, "Specific Discharge (modeled) [mm/day]"]
            return float(val) if not pd.isna(val) else np.nan

        if name == "channel_exit_water__volume_flow_rate":
            if idx < first:
                return np.nan
            val = m.hydrodata.at[idx, "Specific Discharge (modeled) [mm/day]"]
            if pd.isna(val):
                return np.nan
            # mm→m (×1e-3) × km²→m² (×1e6) / day→s (×86400) = ×1e3/86400
            return float(val) * m.drainage_basin_area__km2 * 1e3 / 86400.0

        if name == "snowpack__liquid_equivalent_depth":
            return float(m.snowpack.Hwater) if m.has_snowpack else 0.0

        if name == "subsurface_water__depth":
            if idx < first:
                return np.nan
            val = m.hydrodata.at[idx, "Subsurface storage (modeled total) [mm]"]
            return float(val) if not pd.isna(val) else np.nan

        if name == "land_surface_water__evapotranspiration_volume_flux":
            # Model evapotranspiration flux after water-balance scaling
            # (storage-stress reduction, where enabled, applies separately).
            if idx < first:
                return np.nan
            val = m.hydrodata.at[idx, "ET for model [mm/day]"]
            return float(val) if not pd.isna(val) else np.nan

        # Flux-partition components of the total discharge, recorded by the
        # most recent update().  baseflow is the constant regional-import term
        # (mm/day); it is not part of the reservoir cascade, so the BMI keeps
        # it separate and does NOT add it to the land_surface_water__runoff_
        # volume_flux total (a coupler adds it explicitly).  Note this differs
        # from run_and_score, which DOES fold baseflow_Q (and Nash routing)
        # into its scored discharge as an output-layer post-process; see the
        # class docstring's "Difference from run_and_score discharge".
        if name == "land_surface_water__direct_runoff_volume_flux":
            return np.nan if idx < first else float(m._flux_direct_runoff)

        if name == "land_surface_water__baseflow_volume_flux":
            return np.nan if idx < first else float(m.baseflow_Q)

        if name == "land_surface_water__tile_drain_volume_flux":
            return np.nan if idx < first else float(m._flux_tile)

        if name == "land_surface_water__multipath_drain_volume_flux":
            return np.nan if idx < first else float(m._flux_multipath)

        if name == "land_surface__frozen_ground_index":
            return np.nan if idx < first else float(m._fgi)

        for i, rname in enumerate(_RESERVOIR_DEPTH_NAMES):
            if name == rname:
                if i < len(m.reservoirs):
                    return float(m.reservoirs[i].Hwater)
                return np.nan

        raise KeyError(f"Unknown output variable: {name!r}")

    # -----------------------------------------------------------------------
    # Getters
    # -----------------------------------------------------------------------

    def get_value(self, name: str, dest: np.ndarray) -> np.ndarray:
        """
        Copy the current scalar value of *name* into *dest* and return it.

        For input variables, returns the value at the pending row (the
        value that will be consumed by the next :meth:`update` call).
        For output variables, returns the value written by the most recent
        :meth:`update` call.

        Parameters
        ----------
        name : str
            CSDMS Standard Name.
        dest : numpy.ndarray
            Pre-allocated length-1 array of dtype float64.

        Returns
        -------
        numpy.ndarray
            *dest*, filled in place.
        """
        _check_var(name)
        if name in _INPUT_VAR_NAMES:
            col = _INPUT_COLUMNS[name]
            idx = self._model._timestep_i
            if (col in self._model.hydrodata.columns
                    and idx in self._model.hydrodata.index):
                val = self._model.hydrodata.at[idx, col]
                dest[0] = float(val) if not pd.isna(val) else np.nan
            else:
                dest[0] = np.nan
        else:
            dest[0] = self._scalar_output(name)
        return dest

    def get_value_ptr(self, name: str) -> np.ndarray:
        """
        Return a length-1 float64 array containing the current value.

        Returns a snapshot, not a live pointer — see class docstring.
        """
        _check_var(name)
        dest = np.empty(1, dtype=np.float64)
        return self.get_value(name, dest)

    def get_value_at_indices(
        self, name: str, dest: np.ndarray, inds: np.ndarray
    ) -> np.ndarray:
        """Get value at specific indices (scalar: only index 0 is valid)."""
        if np.any(np.asarray(inds) != 0):
            raise IndexError("Scalar variable has only index 0.")
        return self.get_value(name, dest)

    # -----------------------------------------------------------------------
    # Setters
    # -----------------------------------------------------------------------

    def set_value(self, name: str, src: np.ndarray) -> None:
        """
        Override an input variable for the current (next) timestep.

        Writes ``src[0]`` into the hydrodata DataFrame at the pending row
        (``_timestep_i``), replacing the value read from the CSV.  The
        written value is consumed by the next :meth:`update` call.  This
        is the online-coupling path; in file-driven mode this method need
        not be called.

        Parameters
        ----------
        name : str
            CSDMS Standard Name of an input variable.
        src : numpy.ndarray
            Length-1 float64 array containing the new value.

        Raises
        ------
        KeyError
            If *name* is not a recognised input variable, or if the
            corresponding DataFrame column does not exist in the loaded
            time series.
        """
        if name not in _INPUT_VAR_NAMES:
            raise KeyError(
                f"{name!r} is not a settable input variable. "
                f"Settable inputs: {_INPUT_VAR_NAMES}"
            )
        col = _INPUT_COLUMNS[name]
        if col not in self._model.hydrodata.columns:
            raise KeyError(
                f"Column {col!r} is not present in the loaded time series. "
                f"Add it to the CSV or do not call set_value for {name!r}."
            )
        self._model.hydrodata.at[self._model._timestep_i, col] = float(src[0])

    def set_value_at_indices(
        self, name: str, inds: np.ndarray, src: np.ndarray
    ) -> None:
        """Set value at specific indices (scalar: only index 0 is valid)."""
        if np.any(np.asarray(inds) != 0):
            raise IndexError("Scalar variable has only index 0.")
        self.set_value(name, src)
