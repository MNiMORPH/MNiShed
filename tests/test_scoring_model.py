"""
Tests for mnished.calibration.ScoringModel.

ScoringModel is a build-once / score-many wrapper around run_and_score (it
reads the forcing and constructs the model once, then reuses a fresh copy per
evaluation). The contract is that it is **bit-identical** to the equivalent
run_and_score call — that equivalence is the acceptance test here.
"""

import os

import pytest
import yaml

import mnished
from mnished.calibration import ScoringModel

EXAMPLE_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "examples", "cannon_forward")
)
CANNON_CSV = os.path.join(EXAMPLE_DIR, "CannonTestInput.csv")


def _cfg(et="datafile"):
    return {
        "timeseries": {"datafile": CANNON_CSV},
        "catchment": {"drainage_basin_area__km2": 3800,
                      "evapotranspiration_method": et,
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
        "modules": {"snowpack": True, "frozen_ground": True,
                    "rain_on_snow": True, "direct_runoff": False},
    }


def _write(tmp_path, cfg):
    p = tmp_path / "cfg.yml"
    p.write_text(yaml.safe_dump(cfg))
    return str(p)


@pytest.mark.parametrize("et", ["datafile", "ThornthwaiteChang2019"])
def test_scoring_model_bit_identical(tmp_path, et):
    """score() reproduces run_and_score exactly, for both ET modes (the
    Thornthwaite case also exercises the cached raw-ET path)."""
    cfg = _write(tmp_path, _cfg(et))
    kw = dict(recession_coeff=[20.0, 400.0], melt_factor=1.5, metric="KGE")
    direct = mnished.run_and_score(cfg, enforce_water_balance="water-year",
                                   **kw)
    sm = ScoringModel(cfg, enforce_water_balance="water-year")
    reuse = sm.score(**kw)
    assert reuse.score == direct.score
    assert reuse.aic == direct.aic


def test_scoring_model_et_scale_bit_identical(tmp_path):
    """The cached ET base must still respond exactly to a per-eval et_scale."""
    cfg = _write(tmp_path, _cfg("ThornthwaiteChang2019"))
    sm = ScoringModel(cfg, enforce_water_balance="water-year")
    for es in (0.7, 1.0, 1.3):
        direct = mnished.run_and_score(cfg, enforce_water_balance="water-year",
                                       et_scale=es, metric="KGE")
        assert sm.score(et_scale=es, metric="KGE").score == direct.score


def test_scoring_model_reuse_is_stateless(tmp_path):
    """Repeated calls do not leak state: the same parameters give the same
    score, and different parameters give different scores."""
    sm = ScoringModel(_write(tmp_path, _cfg()),
                      enforce_water_balance="water-year")
    a = sm.score(et_scale=0.7, metric="KGE").score
    b = sm.score(et_scale=1.2, metric="KGE").score
    a_again = sm.score(et_scale=0.7, metric="KGE").score
    assert a != b            # responds to parameters
    assert a_again == a      # no state leak between evaluations
