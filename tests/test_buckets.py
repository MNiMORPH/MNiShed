"""
Regression tests for the Buckets end-to-end model.

Uses the Cannon River forward example (examples/cannon_forward/).

Forcing data are derived from Livneh et al. (2015) gridded meteorology,
extracted for the Cannon River catchment by Pierce et al. (2021) via
https://github.com/MNiMORPH/LivnehPierce-hydro-extractor.

References
----------
Livneh, B., et al. (2015). A spatially comprehensive, hydrometeorological
data set for Mexico, the U.S., and Southern Canada 1950–2013. Scientific
Data, 2, 150042. https://doi.org/10.1038/sdata.2015.42

Pierce, A., et al. (2021). LivnehPierce hydro extractor.
https://github.com/MNiMORPH/LivnehPierce-hydro-extractor
"""

import os

import numpy as np
import pytest

EXAMPLE_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "examples", "cannon_forward")
)
EXAMPLE_CONFIG = "cannon_cfg.yml"


@pytest.fixture(scope="module")
def cannon():
    """Run the Cannon River forward example once; share across tests."""
    import mnished

    orig_dir = os.getcwd()
    os.chdir(EXAMPLE_DIR)
    try:
        b = mnished.Buckets()
        b.initialize(EXAMPLE_CONFIG)
        b.run()
    finally:
        os.chdir(orig_dir)
    return b


# ---------------------------------------------------------------------------
# Smoke / output columns
# ---------------------------------------------------------------------------

def test_run_completes(cannon):
    """Model runs from initialize through run() without error."""
    assert cannon.hydrodata is not None
    assert len(cannon.hydrodata) > 0


def test_output_columns_present(cannon):
    """Expected output columns are present in hydrodata."""
    df = cannon.hydrodata
    assert "Specific Discharge (modeled) [mm/day]" in df.columns
    assert "Snowpack (modeled) [mm SWE]" in df.columns
    assert "Subsurface storage (modeled total) [mm]" in df.columns
    assert "ET for model [mm/day]" in df.columns


def test_discharge_non_negative(cannon):
    """Modeled discharge is non-negative for every timestep."""
    q = cannon.hydrodata["Specific Discharge (modeled) [mm/day]"]
    assert (q >= 0).all()


def test_swe_non_negative(cannon):
    """SWE is non-negative for every timestep."""
    swe = cannon.hydrodata["Snowpack (modeled) [mm SWE]"]
    assert (swe >= 0).all()


def test_storage_non_negative(cannon):
    """Total subsurface storage is non-negative for every timestep."""
    stor = cannon.hydrodata["Subsurface storage (modeled total) [mm]"]
    assert (stor >= 0).all()


# ---------------------------------------------------------------------------
# NSE regression
# ---------------------------------------------------------------------------

def test_nse_regression(cannon):
    """NSE is within 0.01 of the known uncalibrated value (-3.6592)."""
    nse = cannon.compute_NSE(return_nse=True, verbose=False)
    assert nse == pytest.approx(-3.659183, abs=0.01)


# ---------------------------------------------------------------------------
# Water balance
# ---------------------------------------------------------------------------

def test_water_balance_closure(cannon):
    """P - Q - ET - ΔS - ΔSWe < 0.1% of total P."""
    df = cannon.hydrodata
    total_P  = df["Precipitation [mm/day]"].astype(float).sum()
    total_Q  = df["Specific Discharge (modeled) [mm/day]"].astype(float).sum()
    total_ET = df["ET for model [mm/day]"].astype(float).sum()
    stor     = df["Subsurface storage (modeled total) [mm]"].astype(float)
    swe      = df["Snowpack (modeled) [mm SWE]"].astype(float)
    delta_S   = stor.iloc[-1] - stor.iloc[0]
    delta_swe = swe.iloc[-1]  - swe.iloc[0]
    wb_error_fraction = abs(total_P - total_Q - total_ET - delta_S - delta_swe) / total_P
    assert wb_error_fraction < 0.001  # < 0.1% of P


# ---------------------------------------------------------------------------
# Basic physical checks
# ---------------------------------------------------------------------------

def test_discharge_finite(cannon):
    """All modeled discharge values are finite."""
    q = cannon.hydrodata["Specific Discharge (modeled) [mm/day]"].astype(float)
    assert np.isfinite(q).all()


def test_et_positive(cannon):
    """ET is positive over the run period."""
    et = cannon.hydrodata["ET for model [mm/day]"]
    assert et.sum() > 0


def test_snowpack_peaks_in_winter(cannon):
    """SWE peaks in winter (Dec–Feb), not in summer (Jun–Aug)."""
    import pandas as pd
    dates = pd.to_datetime(cannon.hydrodata["Date"])
    swe = cannon.hydrodata["Snowpack (modeled) [mm SWE]"].values
    winter_mask = dates.dt.month.isin([12, 1, 2])
    summer_mask = dates.dt.month.isin([6, 7, 8])
    assert swe[winter_mask].mean() > swe[summer_mask].mean()


def test_three_reservoirs_initialized(cannon):
    """Cannon River config uses 3 reservoirs."""
    assert len(cannon.reservoirs) == 3


def test_record_length(cannon):
    """Cannon River example spans 3 years (1992–1994, ~1096 days including leap year)."""
    n = len(cannon.hydrodata)
    assert 1050 <= n <= 1150


# ---------------------------------------------------------------------------
# AIC free-parameter counting
# ---------------------------------------------------------------------------

def test_hmax_inf_not_counted_in_aic(monkeypatch):
    """An ``.inf`` Hmax entry means 'no saturation cap' and must not count
    toward AIC's free-parameter k.

    Two runs that set identical effective reservoir caps but differ only by
    trailing ``.inf`` Hmax entries must yield the same AIC: same fit, and the
    same number of genuinely calibrated parameters.
    """
    from mnished import run_and_score
    monkeypatch.chdir(EXAMPLE_DIR)
    r_one = run_and_score(EXAMPLE_CONFIG, Hmax=[18.0], spin_up_cycles=1)
    r_inf = run_and_score(EXAMPLE_CONFIG,
                          Hmax=[18.0, float('inf'), float('inf')],
                          spin_up_cycles=1)
    assert np.isfinite(r_one.aic)
    assert r_one.aic == pytest.approx(r_inf.aic)


# ---------------------------------------------------------------------------
# JIT / pure-Python equivalence
# ---------------------------------------------------------------------------

def test_jit_matches_pure_python(monkeypatch):
    """The Numba JIT run() and the pure-Python fallback give the same discharge.

    Skipped when Numba is not installed (the JIT path is unavailable, so both
    runs would be pure-Python). In the Numba-enabled CI job and in local
    development this guards against the two code paths diverging.
    """
    pytest.importorskip("numba")
    import mnished
    import mnished.mnished as _m

    col = "Specific Discharge (modeled) [mm/day]"
    monkeypatch.chdir(EXAMPLE_DIR)

    # JIT path (Numba available and used by run()).
    b_jit = mnished.Buckets()
    b_jit.initialize(EXAMPLE_CONFIG)
    b_jit.run()
    q_jit = b_jit.hydrodata[col].astype(float).to_numpy()

    # Force the pure-Python fallback for the same configuration.
    monkeypatch.setattr(_m, "_numba_available", False)
    b_py = mnished.Buckets()
    b_py.initialize(EXAMPLE_CONFIG)
    b_py.run()
    q_py = b_py.hydrodata[col].astype(float).to_numpy()

    np.testing.assert_allclose(q_jit, q_py, rtol=1e-8, atol=1e-10, equal_nan=True)


def test_jit_matches_pure_python_advanced(tmp_path, monkeypatch):
    """JIT and pure-Python agree with the advanced reservoir mechanics on.

    The basic equivalence test above exercises only linear reservoirs. This
    one enables the v3 branches the JIT reimplements but that the cannon
    config never touches: power-law recession, threshold and leakance
    junctions, a tile drain, and a multipath drain — all at once — to guard
    those JIT code paths against diverging from the pure-Python loop.
    """
    pytest.importorskip("numba")
    import yaml
    import mnished
    import mnished.mnished as _m

    cfg = {
        'timeseries': {'datafile': os.path.join(EXAMPLE_DIR, "CannonTestInput.csv")},
        'initial_conditions': {
            'water_reservoir_effective_depths__mm': [15, 40, 500],
            'snowpack__mm_SWE': 0,
        },
        'catchment': {
            'drainage_basin_area__km2': 3800,
            'evapotranspiration_method': 'datafile',
            'water_year_start_month': 10,
        },
        'general': {'spin_up_cycles': 1, 'direct_runoff_fraction': 0.1,
                    'et_alpha': 0.6},
        'reservoirs': {
            'recession_timescales': [14, 40, 500],
            'exfiltration_fractions': [0.19, 0.76, 1.0],
            'maximum_effective_depths__mm': [18.0, float('inf'), float('inf')],
            'recession_exponents': [2.0, 1.5, 1.0],          # power-law on 0,1
            'junction_types': ['threshold', 'leakance', 'fraction'],
            'leakance_R__days': [None, 50.0, None],
            'H_threshold__mm': [3.0, 0.0, 0.0],
            'tile_fractions': [0.3, 0.0, 0.0],
            'tile_residence_times__days': [3.0, None, None],
            'multipath_thresholds__mm': [None, 20.0, None],
            'multipath_timescales__days': [None, 5.0, None],
        },
        'snowmelt': {
            'PDD_melt_factor': 1.0,
            'fgi_decay_coeff': 0.97,
            'snow_insulation_k': 0.0,
        },
        'modules': {
            'snowpack': True,
            'frozen_ground': True,
            'rain_on_snow': True,
            'direct_runoff': True,
            'et_reservoir_draw': True,
        },
    }
    cfg_path = tmp_path / "advanced.yml"
    cfg_path.write_text(yaml.safe_dump(cfg))

    q_col = "Specific Discharge (modeled) [mm/day]"
    s_col = "Subsurface storage (modeled total) [mm]"

    # wp_soil / wp_soil_sigma are run_and_score-only (not YAML); set directly
    # to exercise the soil wilting-point ET-draw branch.
    b_jit = mnished.Buckets()
    b_jit.initialize(str(cfg_path))
    b_jit.wp_soil = 25.0
    b_jit.wp_soil_sigma = 8.0
    b_jit.run()
    q_jit = b_jit.hydrodata[q_col].astype(float).to_numpy()
    s_jit = b_jit.hydrodata[s_col].astype(float).to_numpy()

    monkeypatch.setattr(_m, "_numba_available", False)
    b_py = mnished.Buckets()
    b_py.initialize(str(cfg_path))
    b_py.wp_soil = 25.0
    b_py.wp_soil_sigma = 8.0
    b_py.run()
    q_py = b_py.hydrodata[q_col].astype(float).to_numpy()
    s_py = b_py.hydrodata[s_col].astype(float).to_numpy()

    np.testing.assert_allclose(q_jit, q_py, rtol=1e-7, atol=1e-9, equal_nan=True)
    np.testing.assert_allclose(s_jit, s_py, rtol=1e-7, atol=1e-9, equal_nan=True)


# ---------------------------------------------------------------------------
# ET reservoir draw: condensation respects Hmax
# ---------------------------------------------------------------------------

def _draw_harness(Hmax, H0):
    """A Buckets with one reservoir, configured for a direct ET-draw call."""
    from mnished import Buckets, Reservoir
    b = Buckets()
    b.reservoirs = [Reservoir(recession_coeff=10.0, Hmax=Hmax, H0=H0)]
    b.et_alpha = 1.0          # all of ET_pot drawn from / added to reservoir 0
    b.wp_soil = 0.0
    b.wp_soil_sigma = 0.0
    return b


def test_et_draw_condensation_above_hmax_runs_off():
    """Negative ET (condensation) above Hmax sheds to runoff, not stored."""
    b = _draw_harness(Hmax=20.0, H0=20.0)
    excess = b._draw_et_from_reservoirs(-5.0)        # 5 mm condensation onto a full reservoir
    assert b.reservoirs[0].Hwater == pytest.approx(20.0)   # capped at Hmax, not 25
    assert excess == pytest.approx(5.0)                    # surplus shed to runoff


def test_et_draw_condensation_below_hmax_is_stored():
    """Condensation that stays under Hmax is stored, with no runoff."""
    b = _draw_harness(Hmax=20.0, H0=15.0)
    excess = b._draw_et_from_reservoirs(-3.0)        # 3 mm condensation, room to spare
    assert b.reservoirs[0].Hwater == pytest.approx(18.0)
    assert excess == pytest.approx(0.0)


def test_et_draw_positive_et_returns_no_excess():
    """A normal (positive ET) draw reduces storage and returns zero excess."""
    b = _draw_harness(Hmax=20.0, H0=15.0)
    excess = b._draw_et_from_reservoirs(4.0)
    assert b.reservoirs[0].Hwater == pytest.approx(11.0)
    assert excess == pytest.approx(0.0)
