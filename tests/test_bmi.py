"""
Tests for the BmiMNiShed wrapper.

Run from the repository root with::

    pytest tests/test_bmi.py -v

The tests use the Cannon River forward example in examples/cannon_forward/.
bmipy must be installed (pip install 'MNiShed[bmi]').
"""

import os

import numpy as np
import pytest

pytest.importorskip("bmipy", reason="bmipy not installed; skipping BMI tests (pip install 'MNiShed[bmi]')")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

EXAMPLE_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "examples", "cannon_forward")
)
EXAMPLE_CONFIG = "cannon_cfg.yml"

P_NAME = "atmosphere_water__liquid_equivalent_precipitation_rate"
T_NAME = "atmosphere__temperature"
Q_NAME = "land_surface_water__runoff_volume_flux"
Q_VOL_NAME = "channel_exit_water__volume_flow_rate"
SWE_NAME = "snowpack__liquid_equivalent_depth"
STOR_NAME = "subsurface_water__depth"
RES0_NAME = "subsurface_water_reservoir_0__depth"
RES1_NAME = "subsurface_water_reservoir_1__depth"
RES2_NAME = "subsurface_water_reservoir_2__depth"


@pytest.fixture
def bmi(monkeypatch):
    """Initialized BmiMNiShed using the Cannon River forward example."""
    from mnished import BmiMNiShed
    monkeypatch.chdir(EXAMPLE_DIR)
    b = BmiMNiShed()
    b.initialize(EXAMPLE_CONFIG)
    yield b
    b.finalize()


# ---------------------------------------------------------------------------
# Lifecycle
# ---------------------------------------------------------------------------

def test_smoke(bmi):
    """initialize → 5× update → finalize completes without error."""
    for _ in range(5):
        bmi.update()
    assert bmi.get_current_time() == pytest.approx(5.0)


def test_update_until(bmi):
    """update_until advances to the requested time."""
    bmi.update_until(10.0)
    assert bmi.get_current_time() == pytest.approx(10.0)


def test_finalize_clears_model(bmi):
    """finalize() sets _model to None."""
    bmi.finalize()
    assert bmi._model is None


# ---------------------------------------------------------------------------
# Component information
# ---------------------------------------------------------------------------

def test_component_name(bmi):
    assert bmi.get_component_name() == "MNiShed"


def test_item_counts(bmi):
    assert bmi.get_input_item_count() == len(bmi.get_input_var_names())
    assert bmi.get_output_item_count() == len(bmi.get_output_var_names())


def test_var_names_are_tuples(bmi):
    assert isinstance(bmi.get_input_var_names(), tuple)
    assert isinstance(bmi.get_output_var_names(), tuple)
    assert P_NAME in bmi.get_input_var_names()
    assert Q_NAME in bmi.get_output_var_names()


# ---------------------------------------------------------------------------
# Variable metadata
# ---------------------------------------------------------------------------

def test_var_units(bmi):
    assert bmi.get_var_units(P_NAME) == "mm d-1"
    assert bmi.get_var_units(T_NAME) == "degC"
    assert bmi.get_var_units(Q_NAME) == "mm d-1"
    assert bmi.get_var_units(SWE_NAME) == "mm"


def test_var_type(bmi):
    for name in bmi.get_input_var_names() + bmi.get_output_var_names():
        assert bmi.get_var_type(name) == "float64"


def test_var_grid_is_zero(bmi):
    for name in bmi.get_input_var_names() + bmi.get_output_var_names():
        assert bmi.get_var_grid(name) == 0


def test_var_itemsize_and_nbytes(bmi):
    for name in bmi.get_input_var_names() + bmi.get_output_var_names():
        assert bmi.get_var_itemsize(name) == 8
        assert bmi.get_var_nbytes(name) == 8


def test_var_location(bmi):
    for name in bmi.get_input_var_names() + bmi.get_output_var_names():
        assert bmi.get_var_location(name) == "node"


def test_unknown_var_raises(bmi):
    with pytest.raises(KeyError):
        bmi.get_var_units("not_a_real_variable")


# ---------------------------------------------------------------------------
# Time
# ---------------------------------------------------------------------------

def test_time_metadata(bmi):
    assert bmi.get_start_time() == pytest.approx(0.0)
    assert bmi.get_end_time() > 0.0
    assert bmi.get_current_time() == pytest.approx(0.0)
    assert bmi.get_time_step() == pytest.approx(1.0)
    assert bmi.get_time_units() == "d"


def test_time_advances(bmi):
    bmi.update()
    assert bmi.get_current_time() == pytest.approx(1.0)
    bmi.update()
    assert bmi.get_current_time() == pytest.approx(2.0)


# ---------------------------------------------------------------------------
# Grid (scalar)
# ---------------------------------------------------------------------------

def test_grid_metadata(bmi):
    assert bmi.get_grid_rank(0) == 0
    assert bmi.get_grid_size(0) == 1
    assert bmi.get_grid_type(0) == "scalar"


def test_unknown_grid_raises(bmi):
    with pytest.raises(ValueError):
        bmi.get_grid_rank(99)


def test_grid_shape_raises(bmi):
    with pytest.raises(NotImplementedError):
        bmi.get_grid_shape(0, np.empty(0))


# ---------------------------------------------------------------------------
# get_value / get_value_ptr
# ---------------------------------------------------------------------------

def test_get_value_discharge_nan_before_update(bmi):
    """Before any update, discharge output is nan."""
    dest = np.empty(1, dtype=np.float64)
    bmi.get_value(Q_NAME, dest)
    assert np.isnan(dest[0])


def test_get_value_discharge_finite_after_update(bmi):
    """After one update, discharge is finite and non-negative."""
    bmi.update()
    dest = np.empty(1, dtype=np.float64)
    bmi.get_value(Q_NAME, dest)
    assert np.isfinite(dest[0])
    assert dest[0] >= 0.0


def test_get_value_volumetric_discharge_after_update(bmi):
    """Volumetric discharge is finite, positive, and consistent with specific discharge."""
    bmi.update()
    dest_q = np.empty(1, dtype=np.float64)
    dest_vol = np.empty(1, dtype=np.float64)
    bmi.get_value(Q_NAME, dest_q)
    bmi.get_value(Q_VOL_NAME, dest_vol)
    assert np.isfinite(dest_vol[0])
    assert dest_vol[0] >= 0.0
    area_km2 = bmi._model.drainage_basin_area__km2
    expected = dest_q[0] * 1e-3 * area_km2 * 1e6 / 86400
    assert dest_vol[0] == pytest.approx(expected)


def test_get_value_swe_after_update(bmi):
    """After one update, SWE is finite and non-negative (snowpack enabled)."""
    bmi.update()
    dest = np.empty(1, dtype=np.float64)
    bmi.get_value(SWE_NAME, dest)
    assert np.isfinite(dest[0])
    assert dest[0] >= 0.0


def test_get_value_storage_after_update(bmi):
    """After one update, total subsurface storage is finite and positive."""
    bmi.update()
    dest = np.empty(1, dtype=np.float64)
    bmi.get_value(STOR_NAME, dest)
    assert np.isfinite(dest[0])
    assert dest[0] > 0.0


def test_get_value_reservoir_depths(bmi):
    """Per-reservoir depths are finite and non-negative after one update."""
    bmi.update()
    dest = np.empty(1, dtype=np.float64)
    for name in (RES0_NAME, RES1_NAME, RES2_NAME):
        bmi.get_value(name, dest)
        # Cannon example has 3 reservoirs; all should be finite
        assert np.isfinite(dest[0])
        assert dest[0] >= 0.0


def test_get_value_returns_dest(bmi):
    """get_value returns the dest array (BMI convention)."""
    bmi.update()
    dest = np.empty(1, dtype=np.float64)
    result = bmi.get_value(Q_NAME, dest)
    assert result is dest


def test_get_value_ptr_returns_array(bmi):
    bmi.update()
    arr = bmi.get_value_ptr(Q_NAME)
    assert isinstance(arr, np.ndarray)
    assert arr.shape == (1,)
    assert arr.dtype == np.float64


def test_get_value_at_indices(bmi):
    bmi.update()
    dest = np.empty(1, dtype=np.float64)
    bmi.get_value_at_indices(Q_NAME, dest, np.array([0]))
    assert np.isfinite(dest[0])


def test_get_value_at_indices_bad_index(bmi):
    bmi.update()
    dest = np.empty(1, dtype=np.float64)
    with pytest.raises(IndexError):
        bmi.get_value_at_indices(Q_NAME, dest, np.array([1]))


# ---------------------------------------------------------------------------
# set_value
# ---------------------------------------------------------------------------

def test_set_value_writes_to_dataframe(bmi):
    """set_value overwrites the pending row in the forcing DataFrame."""
    idx = bmi._model._timestep_i
    bmi.set_value(P_NAME, np.array([999.0]))
    col = "Precipitation [mm/day]"
    assert bmi._model.hydrodata.at[idx, col] == pytest.approx(999.0)


def test_set_value_then_update_uses_overridden_value(bmi):
    """Precipitation written via set_value is used by the next update."""
    idx = bmi._model._timestep_i
    bmi.set_value(P_NAME, np.array([999.0]))
    bmi.update()
    # The row that was just processed should have 999 mm/day precipitation
    assert bmi._model.hydrodata.at[idx, "Precipitation [mm/day]"] == pytest.approx(999.0)


def test_set_value_high_precip_increases_discharge(monkeypatch):
    """End-to-end test: very high precipitation raises discharge.

    Uses step 182 (1992-07-01, T = 21 °C) so that by construction:
    - No snowpack exists (all winter snow melted months earlier).
    - New precipitation falls entirely as liquid into the soil reservoir.
    - Snowmelt is zero (no SWE to melt) and identical in both runs.

    The only difference between the two runs is the precipitation on step 182.
    The resulting discharge difference verifies the full set_value → update
    → get_value pipeline.
    """
    from mnished import BmiMNiShed
    monkeypatch.chdir(EXAMPLE_DIR)

    SUMMER_STEP = 182  # 1992-07-01, T = 21 °C, no snowpack

    # Normal run through step SUMMER_STEP
    bmi_base = BmiMNiShed()
    bmi_base.initialize(EXAMPLE_CONFIG)
    bmi_base.update_until(float(SUMMER_STEP + 1))
    dest_normal = np.empty(1, dtype=np.float64)
    bmi_base.get_value(Q_NAME, dest_normal)
    bmi_base.finalize()

    # Same run but with very high precipitation on step SUMMER_STEP
    bmi_high = BmiMNiShed()
    bmi_high.initialize(EXAMPLE_CONFIG)
    bmi_high.update_until(float(SUMMER_STEP))
    bmi_high.set_value(P_NAME, np.array([500.0]))
    bmi_high.update()
    dest_high = np.empty(1, dtype=np.float64)
    bmi_high.get_value(Q_NAME, dest_high)
    bmi_high.finalize()

    assert dest_high[0] > dest_normal[0]


def test_set_value_bad_name_raises(bmi):
    """set_value raises KeyError for an output or unknown variable name."""
    with pytest.raises(KeyError):
        bmi.set_value(Q_NAME, np.array([1.0]))  # Q is output, not input


def test_set_value_at_indices(bmi):
    bmi.set_value_at_indices(P_NAME, np.array([0]), np.array([7.5]))
    idx = bmi._model._timestep_i
    assert bmi._model.hydrodata.at[idx, "Precipitation [mm/day]"] == pytest.approx(7.5)
