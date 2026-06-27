"""
Tests for the GDD vegetation-phenology coefficient (Buckets.phenology_Kc).

The phenology factor reshapes the Thornthwaite ET demand to follow thermal-time
leaf-out instead of temperature: it suppresses early-spring ET (so a snowmelt
basin's freshet is not evaporated before leaf-out) and lets the water-balance
correction re-close the annual total. The contract tested here:

* disabled (the default) is a no-op — bit-identical to a run without it;
* enabled lowers spring ET and, under global closure, preserves the annual total;
* the coefficient stays within ``[dormant_Kc, full_Kc]`` and leafs out in spring;
* with measured (``datafile``) ET it is ignored, with a warning.
"""

import os
import warnings

import numpy as np
import pytest
import yaml

import mnished
from mnished import Buckets

EXAMPLE_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "examples", "cannon_forward")
)
CANNON_CSV = os.path.join(EXAMPLE_DIR, "CannonTestInput.csv")


def _cfg(et="ThornthwaiteChang2019", phenology=None):
    cfg = {
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
        "modules": {"snowpack": True, "frozen_ground": False,
                    "rain_on_snow": True, "direct_runoff": False},
    }
    if phenology is not None:
        cfg["phenology"] = phenology
    return cfg


def _et_series(tmp_path, cfg, wb="global"):
    """The model ET column and its months for a built+run Buckets."""
    p = tmp_path / "cfg.yml"
    p.write_text(yaml.safe_dump(cfg))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Buckets()
        b.initialize(str(p), enforce_water_balance=wb)
        b.run()
    et = b.hydrodata["ET for model [mm/day]"].to_numpy()
    months = b.hydrodata["Date"].dt.month.to_numpy()
    return et, months


def test_phenology_disabled_is_noop(tmp_path):
    """Absent phenology block == explicit enabled:false, bit-identical."""
    et_absent, _ = _et_series(tmp_path, _cfg())
    et_false, _ = _et_series(tmp_path, _cfg(phenology={"enabled": False}))
    assert np.allclose(et_absent, et_false, equal_nan=True)


def test_phenology_reshapes_and_preserves_annual_total(tmp_path):
    """Enabled suppresses spring ET, raises summer ET, and keeps the annual
    total fixed under global closure (the freshet-preserving reshape)."""
    et_off, mo = _et_series(tmp_path, _cfg())
    et_on, _ = _et_series(tmp_path, _cfg(phenology={"enabled": True}))

    spring = np.isin(mo, [3, 4, 5])
    summer = np.isin(mo, [6, 7, 8])
    assert np.nanmean(et_on[spring]) < np.nanmean(et_off[spring])   # freshet kept
    assert np.nanmean(et_on[summer]) > np.nanmean(et_off[summer])   # compensates
    # April (the freshet month) is strongly suppressed
    apr = mo == 4
    assert np.nanmean(et_on[apr]) < 0.6 * np.nanmean(et_off[apr])
    # global closure preserves the long-term total
    assert np.isclose(np.nansum(et_on), np.nansum(et_off), rtol=1e-6)


def test_phenology_Kc_bounds_and_spring_leafout(tmp_path):
    """Kc stays in [dormant_Kc, full_Kc] and rises from April to June."""
    cfg = _cfg(phenology={"enabled": True, "dormant_Kc": 0.4, "full_Kc": 1.0})
    p = tmp_path / "cfg.yml"
    p.write_text(yaml.safe_dump(cfg))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Buckets()
        b.initialize(str(p), enforce_water_balance="global")
    Kc = np.asarray(b.phenology_Kc())
    mo = b.hydrodata["Date"].dt.month.to_numpy()
    assert Kc.min() >= 0.4 - 1e-9 and Kc.max() <= 1.0 + 1e-9
    # dormant in early spring, full by midsummer: April < June
    assert Kc[mo == 4].mean() < Kc[mo == 6].mean()
    assert np.isclose(Kc[mo == 7].mean(), 1.0, atol=1e-6)   # full canopy summer


def test_phenology_ignored_for_datafile(tmp_path):
    """With measured ET, phenology is ignored (score unchanged) and warns."""
    base = _cfg(et="datafile")
    s0 = mnished.run_and_score(_write(tmp_path, base, "a"),
                               enforce_water_balance="global", metric="KGE").score
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        s1 = mnished.run_and_score(
            _write(tmp_path, _cfg(et="datafile", phenology={"enabled": True}), "b"),
            enforce_water_balance="global", metric="KGE").score
    assert s0 == s1
    assert any("phenology Kc" in str(w.message) for w in caught)


def _write(tmp_path, cfg, name):
    p = tmp_path / f"cfg_{name}.yml"
    p.write_text(yaml.safe_dump(cfg))
    return str(p)


@pytest.mark.parametrize("wb", ["global", "none"])
def test_phenology_runs_under_scalar_closure(tmp_path, wb):
    """Phenology composes with global and none/et_scale closure without error."""
    et, _ = _et_series(tmp_path, _cfg(phenology={"enabled": True}), wb=wb)
    assert np.isfinite(et).any()
