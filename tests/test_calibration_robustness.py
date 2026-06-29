"""
Robustness tests for run_and_score's handling of chained initial states.

Chaining initial conditions from a partial-data or failed decade can produce
NaN end-states. Passed on unchecked, they propagate silently through the run
(every modelled flow becomes NaN and the score looks merely poor rather than
broken). run_and_score validates the chained `initial_states` /
`post_spinup_states` at the boundary and raises a clear error instead.
"""

import math
import os

import pytest
import yaml

import mnished

EXAMPLE_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "examples", "cannon_forward")
)
CANNON_CSV = os.path.join(EXAMPLE_DIR, "CannonTestInput.csv")

NAN = float("nan")


def _legacy_cfg():
    """A single-cascade (K=1) config: two reservoirs, Cannon forcing."""
    return {
        'timeseries': {'datafile': CANNON_CSV},
        'catchment': {
            'drainage_basin_area__km2': 3800,
            'evapotranspiration_method': 'datafile',
            'water_year_start_month': 10,
        },
        'general': {'spin_up_cycles': 0},
        'reservoirs': {
            'recession_coefficients': [14, 500],
            'exfiltration_fractions': [0.3, 1.0],
            'maximum_effective_depths__mm': [float('inf'), float('inf')],
        },
        'initial_conditions': {
            'water_reservoir_effective_depths__mm': [15, 400],
            'snowpack__mm_SWE': 0,
        },
        'snowmelt': {'PDD_melt_factor': 1.0, 'fgi_decay_coeff': 0.97,
                     'snow_insulation_k': 0.0},
        'modules': {'snowpack': True, 'frozen_ground': True,
                    'rain_on_snow': True, 'direct_runoff': False},
    }


def _write(tmp_path, cfg):
    p = tmp_path / "cfg.yml"
    p.write_text(yaml.safe_dump(cfg))
    return str(p)


# ---------------------------------------------------------------------------
# Non-finite chained states raise at the boundary
# ---------------------------------------------------------------------------

def test_nan_reservoir_in_initial_states_raises(tmp_path):
    cfg_path = _write(tmp_path, _legacy_cfg())
    with pytest.raises(ValueError, match=r"initial_states\['reservoirs'\]\[1\].*non-finite"):
        mnished.run_and_score(
            cfg_path, spin_up_cycles=0, metric='KGE',
            initial_states={'reservoirs': [20.0, NAN], 'snowpack': 0.0,
                            'fgi': 0.0, 'H_deficit_carry': 0.0})


def test_nan_snowpack_in_initial_states_raises(tmp_path):
    cfg_path = _write(tmp_path, _legacy_cfg())
    with pytest.raises(ValueError, match=r"\['snowpack'\].*non-finite"):
        mnished.run_and_score(
            cfg_path, spin_up_cycles=0, metric='KGE',
            initial_states={'reservoirs': [20.0, 300.0], 'snowpack': NAN,
                            'fgi': 0.0})


def test_inf_fgi_in_initial_states_raises(tmp_path):
    cfg_path = _write(tmp_path, _legacy_cfg())
    with pytest.raises(ValueError, match=r"\['fgi'\].*non-finite"):
        mnished.run_and_score(
            cfg_path, spin_up_cycles=0, metric='KGE',
            initial_states={'reservoirs': [20.0, 300.0], 'fgi': float('inf')})


def test_nan_in_nested_initial_states_names_sub_catchment(tmp_path):
    cfg_path = _write(tmp_path, _legacy_cfg())
    with pytest.raises(ValueError,
                       match=r"initial_states\['sub_catchments'\]\[0\]\['reservoirs'\]\[0\]"):
        mnished.run_and_score(
            cfg_path, spin_up_cycles=0, metric='KGE',
            initial_states={'sub_catchments': [
                {'reservoirs': [NAN, 300.0], 'snowpack': 0.0, 'fgi': 0.0,
                 'H_deficit_carry': 0.0}]})


def test_nan_in_post_spinup_states_raises(tmp_path):
    cfg_path = _write(tmp_path, _legacy_cfg())
    with pytest.raises(ValueError, match=r"post_spinup_states\['reservoirs'\]\[1\].*non-finite"):
        mnished.run_and_score(
            cfg_path, start='1993-01-01', end='1994-12-31', spin_up_cycles=1,
            metric='KGE', post_spinup_states={'reservoirs': [None, NAN]})


# ---------------------------------------------------------------------------
# Valid states still run; None reservoir entries are allowed
# ---------------------------------------------------------------------------

def test_finite_initial_states_runs(tmp_path):
    cfg_path = _write(tmp_path, _legacy_cfg())
    result = mnished.run_and_score(
        cfg_path, spin_up_cycles=0, metric='KGE',
        initial_states={'reservoirs': [20.0, 300.0], 'snowpack': 0.0,
                        'fgi': 0.0, 'H_deficit_carry': 0.0})
    assert isinstance(result.final_states, dict)
    assert all(math.isfinite(h) for h in result.final_states['reservoirs'])


def test_none_reservoir_entries_in_post_spinup_allowed(tmp_path):
    """post_spinup_states may use None to keep a reservoir's spin-up value;
    that is not a non-finite value and must not raise."""
    cfg_path = _write(tmp_path, _legacy_cfg())
    result = mnished.run_and_score(
        cfg_path, start='1993-01-01', end='1994-12-31', spin_up_cycles=1,
        metric='KGE', post_spinup_states={'reservoirs': [None, 350.0]})
    assert isinstance(result.final_states, dict)


# --------------------------------------------------------------------------
# Seasonally-weighted objective (issue #37): equal weight per meteorological
# season, so a whole-record-flat fit can't trade a high winter for a high fall.
# --------------------------------------------------------------------------

def _write_cfg(tmp_path, cfg):
    p = tmp_path / "cfg.yml"
    p.write_text(yaml.safe_dump(cfg))
    return str(p)


def test_seasonal_metric_equals_mean_of_per_season_scores(tmp_path):
    import numpy as np
    import pandas as pd
    from mnished.calibration import _kge_logkge, _seasonal_mean

    path = _write_cfg(tmp_path, _legacy_cfg())
    res = mnished.run_and_score(path, enforce_water_balance="global",
                                metric="KGE_logKGE_seasonal")
    assert math.isfinite(res.score)

    # reconstruct the equal-weight seasonal mean by hand from the scored rows
    hd = res.buckets.hydrodata
    q_mod = pd.to_numeric(hd["Specific Discharge (modeled) [mm/day]"],
                          errors="coerce")
    q_obs = hd["Specific Discharge [mm/day]"]
    mask = q_mod.notna() & q_obs.notna()
    m = q_mod[mask].to_numpy(float)
    o = q_obs[mask].to_numpy(float)
    months = pd.DatetimeIndex(hd["Date"][mask]).month.to_numpy()
    assert np.isclose(res.score, _seasonal_mean(_kge_logkge, m, o, months))


def test_seasonal_metric_differs_from_whole_record(tmp_path):
    path = _write_cfg(tmp_path, _legacy_cfg())
    seasonal = mnished.run_and_score(path, enforce_water_balance="global",
                                     metric="KGE_logKGE_seasonal").score
    whole = mnished.run_and_score(path, enforce_water_balance="global",
                                  metric="KGE_logKGE").score
    assert not math.isclose(seasonal, whole)


def test_unknown_metric_lists_seasonal_option(tmp_path):
    path = _write_cfg(tmp_path, _legacy_cfg())
    with pytest.raises(ValueError, match="KGE_logKGE_seasonal"):
        mnished.run_and_score(path, metric="not_a_metric")
