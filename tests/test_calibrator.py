"""
Tests for mnished.calibration: target_kwargs (the declarative parameter
mapping) and Calibrator (the build-once, config-driven calibration problem).

The contract is that a Calibrator scores **bit-identically** to the equivalent
hand-written run_and_score call — the declarative config is just a different way
to spell the same model run.
"""

import os

import yaml

from mnished import Calibrator, run_and_score, target_kwargs

EXAMPLE_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "examples", "cannon_forward")
)
CANNON_CSV = os.path.join(EXAMPLE_DIR, "CannonTestInput.csv")


# --- target_kwargs: the declarative mapping -----------------------------

def test_target_kwargs_flat():
    params = {
        "log__a": {"fixed": 1.0, "target": "recession_coeff[0]"},
        "log__b": {"fixed": 2.0, "target": "recession_coeff[1]"},
        "c":      {"fixed": 0.4, "target": "f_to_discharge[0]"},
        "d":      {"fixed": 1.5, "target": "melt_factor"},
        "off":    {"fixed": 9.9},                       # no target -> ignored
    }
    kw = target_kwargs(params, {"log__a": 1.3})         # log__a overridden
    assert kw["recession_coeff"] == [10 ** 1.3, 10 ** 2.0]
    assert kw["f_to_discharge"] == [0.4]
    assert kw["melt_factor"] == 1.5
    assert "off" not in kw and "sub_catchments" not in kw


def test_target_kwargs_nested_shared_and_lake():
    params = {
        "log__soil": {"fixed": 1.3, "target": "sub_catchments[0,1].recession_coeff[0]"},
        "log__gw":   {"fixed": 3.5, "target": "sub_catchments[0,1].recession_coeff[1]"},
        "f_exf":     {"fixed": 0.5, "target": "sub_catchments[0,1].f_to_discharge[0]"},
        "log__lake": {"fixed": 2.5, "target": "sub_catchments[2].recession_coeff[0]"},
        "f_route":   {"fixed": 0.6, "target": "sub_catchments[2].f_route_lake"},
    }
    sub = target_kwargs(params, {}, n_sub=3)["sub_catchments"]
    assert sub[0] == sub[1]                              # shared land zones
    assert sub[0]["recession_coeff"] == [10 ** 1.3, 10 ** 3.5]
    assert sub[0]["f_to_discharge"] == [0.5]            # short list keeps cfg [1]
    assert sub[2]["f_route_lake"] == 0.6
    assert sub[2]["recession_coeff"] == [10 ** 2.5]


# --- Calibrator: bit-identical to run_and_score -------------------------

def _flat_model(tmp_path):
    cfg = {
        "timeseries": {"datafile": CANNON_CSV},
        "catchment": {"drainage_basin_area__km2": 3800,
                      "evapotranspiration_method": "datafile",
                      "water_year_start_month": 10},
        "general": {"spin_up_cycles": 0},
        "reservoirs": {"recession_coefficients": [14, 500],
                       "exfiltration_fractions": [0.3, 1.0],
                       "maximum_effective_depths__mm": [float("inf"),
                                                        float("inf")]},
        "initial_conditions": {
            "water_reservoir_effective_depths__mm": [15, 400],
            "snowpack__mm_SWE": 0},
        "snowmelt": {"PDD_melt_factor": 1.0},
        "modules": {"snowpack": True, "frozen_ground": False,
                    "rain_on_snow": True, "direct_runoff": False},
    }
    p = tmp_path / "model.yml"
    p.write_text(yaml.safe_dump(cfg))
    return str(p), cfg["modules"]


def _flat_params(model_path):
    return {
        "parameters": {
            "log__rec_shallow": {"lower": 0.3, "upper": 2.0, "initial": 1.2,
                                 "fixed": 1.2, "active": True,
                                 "target": "recession_coeff[0]"},
            "log__rec_deep": {"lower": 2.0, "upper": 3.0, "initial": 2.6,
                              "fixed": 2.6, "active": True,
                              "target": "recession_coeff[1]"},
            "f_exf": {"lower": 0.01, "upper": 0.99, "initial": 0.3,
                      "fixed": 0.3, "active": True,
                      "target": "f_to_discharge[0]"},
            "melt": {"lower": 0.1, "upper": 5.0, "initial": 1.0, "fixed": 1.0,
                     "active": True, "target": "melt_factor"},
        },
        "driver": {"config_template": model_path, "metric": "KGE",
                   "spin_up_cycles": 0, "routing_N": 2},
    }


def test_calibrator_bit_identical(tmp_path):
    model_path, modules = _flat_model(tmp_path)
    cfg = _flat_params(model_path)
    cal = Calibrator(cfg["parameters"], cfg["driver"], modules)
    theta = {"log__rec_shallow": 1.4, "log__rec_deep": 2.7, "f_exf": 0.4,
             "melt": 1.6}
    hand = run_and_score(model_path, modules=modules, routing_N=2,
                         spin_up_cycles=0, metric="KGE",
                         recession_coeff=[10 ** 1.4, 10 ** 2.7],
                         f_to_discharge=[0.4], melt_factor=1.6)
    assert cal.score(theta).score == hand.score


def test_calibrator_vector_equals_dict_and_from_yaml(tmp_path):
    model_path, modules = _flat_model(tmp_path)
    cfg = _flat_params(model_path)
    cfg["modules"] = modules
    pyml = tmp_path / "params.yml"
    pyml.write_text(yaml.safe_dump(cfg))
    cal = Calibrator.from_yaml(str(pyml))
    assert set(cal.names) == {"log__rec_shallow", "log__rec_deep", "f_exf",
                              "melt"}
    theta = {"log__rec_shallow": 1.4, "log__rec_deep": 2.7, "f_exf": 0.4,
             "melt": 1.6}
    vec = [theta[n] for n in cal.names]
    assert cal.score(vec).score == cal.score(theta).score


def test_calibrator_nested_bit_identical(tmp_path):
    """A two-land-zone + lake config calibrated declaratively matches the
    hand-written sub_catchments override."""
    cfg = {
        "timeseries": {"datafile": CANNON_CSV},
        "catchment": {"drainage_basin_area__km2": 3800,
                      "evapotranspiration_method": "datafile",
                      "water_year_start_month": 10},
        "general": {"spin_up_cycles": 0, "enforce_water_balance": "none"},
        "snowmelt": {"PDD_melt_factor": 1.0},
        "modules": {"snowpack": False, "frozen_ground": False,
                    "rain_on_snow": False, "direct_runoff": False},
        "sub_catchments": [
            {"name": "direct_land", "area_fraction": 0.3,
             "reservoirs": {"recession_coefficients": [14, 500],
                            "exfiltration_fractions": [0.3, 1.0],
                            "maximum_effective_depths__mm": [20.0,
                                                             float("inf")]},
             "initial_conditions": {
                 "water_reservoir_effective_depths__mm": [8, 350]}},
            {"name": "lake_basin_land", "area_fraction": 0.4,
             "reservoirs": {"recession_coefficients": [14, 500],
                            "exfiltration_fractions": [0.3, 1.0],
                            "maximum_effective_depths__mm": [20.0,
                                                             float("inf")]},
             "initial_conditions": {
                 "water_reservoir_effective_depths__mm": [8, 350]}},
            {"name": "lake", "kind": "lake", "area_fraction": 0.3,
             "lake": {"outflow_coefficient": 0.05, "sill_storage__mm": 180.0,
                      "gw_partner": "lake_basin_land", "f_route_lake": 0.0},
             "initial_conditions": {"lake_storage__mm": 260.0}}],
    }
    model_path = str(tmp_path / "lake.yml")
    open(model_path, "w").write(yaml.safe_dump(cfg))
    params = {
        "log__soil": {"lower": 0.0, "upper": 3.0, "initial": 1.3, "fixed": 1.3,
                      "active": True,
                      "target": "sub_catchments[0,1].recession_coeff[0]"},
        "log__gw": {"lower": 1.0, "upper": 5.0, "initial": 3.5, "fixed": 3.5,
                    "active": True,
                    "target": "sub_catchments[0,1].recession_coeff[1]"},
        "f_exf": {"lower": 0.01, "upper": 0.99, "initial": 0.5, "fixed": 0.5,
                  "active": True, "target": "sub_catchments[0,1].f_to_discharge[0]"},
        "log__lake": {"lower": 0.3, "upper": 7.0, "initial": 2.5, "fixed": 2.5,
                      "active": True,
                      "target": "sub_catchments[2].recession_coeff[0]"},
        "f_route": {"lower": 0.0, "upper": 1.0, "initial": 0.6, "fixed": 0.6,
                    "active": True, "target": "sub_catchments[2].f_route_lake"},
    }
    driver = {"config_template": model_path, "metric": "KGE",
              "spin_up_cycles": 0, "routing_N": 2,
              "enforce_water_balance": "none"}
    cal = Calibrator(params, driver, cfg["modules"])
    land = {"recession_coeff": [10 ** 1.3, 10 ** 3.5], "f_to_discharge": [0.5, 1.0]}
    lake = {"recession_coeff": [10 ** 2.5], "f_route_lake": 0.6}
    hand = run_and_score(model_path, modules=cfg["modules"], routing_N=2,
                         spin_up_cycles=0, metric="KGE",
                         enforce_water_balance="none",
                         sub_catchments=[land, land, lake])
    assert cal.score({n: params[n]["initial"] for n in cal.names}).score == \
        hand.score


# --- Calibrator.score_windows: multi-window objective -------------------

def test_score_windows_per_window(tmp_path):
    """score_windows returns one CalibResult per driver window, each equal to the
    single-window score() for that span."""
    model_path, modules = _flat_model(tmp_path)
    cfg = _flat_params(model_path)
    windows = [{"start": "1993-01-01", "end": "1993-12-31"},
               {"start": "1994-01-01", "end": "1994-12-31"}]
    cfg["driver"]["decades"] = windows
    cal = Calibrator(cfg["parameters"], cfg["driver"], modules)
    theta = {"log__rec_shallow": 1.4, "log__rec_deep": 2.7, "f_exf": 0.4,
             "melt": 1.6}
    results = cal.score_windows(theta)
    assert len(results) == 2
    for w, r in zip(windows, results):
        assert r.score == cal.score(theta, start=w["start"], end=w["end"]).score
    assert results[0].score != results[1].score          # disjoint years differ


def test_score_windows_single_default_and_override(tmp_path):
    """Without a decades: list, score_windows is one full-record window == score();
    an explicit windows= overrides the driver default."""
    model_path, modules = _flat_model(tmp_path)
    cfg = _flat_params(model_path)
    cal = Calibrator(cfg["parameters"], cfg["driver"], modules)
    theta = {"log__rec_shallow": 1.4, "log__rec_deep": 2.7, "f_exf": 0.4,
             "melt": 1.6}
    default = cal.score_windows(theta)
    assert len(default) == 1
    assert default[0].score == cal.score(theta).score
    override = cal.score_windows(
        theta, windows=[{"start": "1993-01-01", "end": "1993-12-31"},
                        {"start": "1994-01-01", "end": "1994-12-31"}])
    assert len(override) == 2
    assert override[0].score == cal.score(theta, start="1993-01-01",
                                          end="1993-12-31").score
