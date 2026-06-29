"""
Tests for mnished.diagnostics.SeasonalMassBalance and the store_fluxes
per-source discharge partition.

The core invariant is that the recorded fast / slow / lake discharge sums exactly
to the modeled specific discharge — both per time step and after seasonal
aggregation — so the decomposition is an exact partition of the hydrograph, not
an approximation. Also covered: the lake-routing case (f_route_lake > 0), and the
clear error when a run did not record the partition.
"""

import os
import warnings

import numpy as np
import pandas as pd
import pytest
import yaml

from mnished import Buckets, SeasonalMassBalance, run_and_score

EXAMPLE_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "examples", "cannon_forward")
)
CANNON_CSV = os.path.join(EXAMPLE_DIR, "CannonTestInput.csv")

FLUX = ('Discharge: fast [mm/day]', 'Discharge: slow [mm/day]',
        'Discharge: lake [mm/day]')


def _flat_cfg():
    return {
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


def _lake_cfg():
    return {
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
                      "gw_partner": "lake_basin_land", "f_route_lake": 0.5},
             "initial_conditions": {"lake_storage__mm": 260.0}}],
    }


def _run(cfg, wb):
    p = "/tmp/_smb_cfg.yml"
    with open(p, "w") as f:
        yaml.safe_dump(cfg, f)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Buckets()
        b.initialize(p, enforce_water_balance=wb)
        b.run(store_fluxes=True)
    return b


def test_partition_sums_to_modeled_per_step(tmp_path):
    b = _run(_flat_cfg(), "water-year")
    hd = b.hydrodata
    mod = hd["Specific Discharge (modeled) [mm/day]"].to_numpy(dtype=float)
    parts = sum(hd[c].to_numpy(dtype=float) for c in FLUX)
    m = np.isfinite(mod)
    assert np.allclose(parts[m], mod[m], atol=1e-9)


def test_partition_sums_with_lake_routing(tmp_path):
    b = _run(_lake_cfg(), "none")
    hd = b.hydrodata
    mod = hd["Specific Discharge (modeled) [mm/day]"].to_numpy(dtype=float)
    lake = hd["Discharge: lake [mm/day]"].to_numpy(dtype=float)
    parts = sum(hd[c].to_numpy(dtype=float) for c in FLUX)
    m = np.isfinite(mod)
    assert np.allclose(parts[m], mod[m], atol=1e-9)
    assert np.nansum(lake[m]) > 0          # lake contributes flow


def test_seasonal_table_is_an_exact_partition(tmp_path):
    b = _run(_flat_cfg(), "water-year")
    st = SeasonalMassBalance(b).seasonal_table()
    assert list(st.index) == ["DJF", "MAM", "JJA", "SON"]
    assert np.allclose((st["fast"] + st["slow"] + st["lake"]).to_numpy(),
                       st["mod"].to_numpy(), atol=1e-9)
    # mod/obs is the ratio of the mod and obs columns
    assert np.allclose((st["mod"] / st["obs"]).to_numpy(),
                       st["mod/obs"].to_numpy(), equal_nan=True)


def test_monthly_table_twelve_rows(tmp_path):
    smb = SeasonalMassBalance(_run(_flat_cfg(), "water-year"))
    mt = smb.monthly_table()
    assert list(mt.index) == list(range(1, 13))
    assert {"SWE", "P", "ET", "obs", "mod"} <= set(mt.columns)


def test_window_restriction(tmp_path):
    b = _run(_flat_cfg(), "water-year")
    full = SeasonalMassBalance(b)
    win = SeasonalMassBalance(b, start="1993-01-01", end="1994-12-31")
    assert len(win.df) < len(full.df)


def test_requires_store_fluxes(tmp_path):
    """A run without store_fluxes gives a clear, actionable error."""
    p = tmp_path / "cfg.yml"
    p.write_text(yaml.safe_dump(_flat_cfg()))
    res = run_and_score(str(p), enforce_water_balance="water-year", metric="KGE")
    with pytest.raises(ValueError, match="store_fluxes=True"):
        SeasonalMassBalance(res.buckets)


def test_run_and_score_passthrough(tmp_path):
    """run_and_score(store_fluxes=True) records the partition on .buckets."""
    p = tmp_path / "cfg.yml"
    p.write_text(yaml.safe_dump(_flat_cfg()))
    res = run_and_score(str(p), enforce_water_balance="water-year",
                        metric="KGE", store_fluxes=True)
    smb = SeasonalMassBalance(res.buckets)
    assert "DJF" in smb.report()


def test_partition_closes_with_phenology(tmp_path):
    """The fast/slow/lake partition is an exact partition of modeled Q even when
    phenology reshapes the ET demand (the split is ET-method-independent)."""
    cfg = _flat_cfg()
    cfg["catchment"]["evapotranspiration_method"] = "ThornthwaiteChang2019"
    cfg["phenology"] = {"enabled": True}
    b = _run(cfg, "global")
    hd = b.hydrodata
    mod = hd["Specific Discharge (modeled) [mm/day]"].to_numpy(dtype=float)
    parts = sum(hd[c].to_numpy(dtype=float) for c in FLUX)
    m = np.isfinite(mod)
    assert np.allclose(parts[m], mod[m], atol=1e-9)


def test_seasonal_et_includes_lake_evaporation(tmp_path):
    """The diagnostic's ET term is area-weighted to include lake open-water
    evaporation, so for a phenology + lake basin it differs from the land ET."""
    cfg = _lake_cfg()
    cfg["catchment"]["evapotranspiration_method"] = "ThornthwaiteChang2019"
    cfg["phenology"] = {"enabled": True}
    b = _run(cfg, "none")
    smb = SeasonalMassBalance(b)
    hd = b.hydrodata
    land = pd.to_numeric(hd["ET for model [mm/day]"], errors="coerce").to_numpy()
    ow = pd.to_numeric(hd["ET for model (open water) [mm/day]"],
                       errors="coerce").to_numpy()
    lake_frac = sum(sc.area_fraction for sc in b.sub_catchments
                    if sc.kind == "lake")
    expect = (1 - lake_frac) * land + lake_frac * ow
    # the diagnostic ET (over its scored rows) matches the area-weighted basin ET
    et_diag = smb.df["ET"].to_numpy()
    expect_scored = expect[hd["Specific Discharge (modeled) [mm/day]"].notna()]
    assert np.allclose(et_diag, expect_scored, equal_nan=True)
    assert not np.allclose(et_diag, land[hd["Specific Discharge (modeled) [mm/day]"].notna()],
                           equal_nan=True)   # lake E actually shifts it


def test_store_fluxes_partition_exact_under_routing(tmp_path):
    """fast/slow/lake stays an exact decomposition of the *routed* discharge: the
    Nash cascade is linear, so routing each source through the same cascade keeps
    the partition consistent with 'Specific Discharge (modeled)'. Regression for
    MNiMORPH/MNiShed#36.2."""
    p = tmp_path / "cfg.yml"
    p.write_text(yaml.safe_dump(_flat_cfg()))
    res = run_and_score(str(p), enforce_water_balance="water-year", metric="KGE",
                        store_fluxes=True, routing_K=10.0)
    hd = res.buckets.hydrodata
    mod = pd.to_numeric(hd["Specific Discharge (modeled) [mm/day]"],
                        errors="coerce").to_numpy()
    parts = sum(pd.to_numeric(hd[c], errors="coerce").to_numpy() for c in FLUX)
    m = np.isfinite(mod)
    assert np.allclose(parts[m], mod[m], atol=1e-9)         # exact after routing
    # and routing actually changed the series (so the test is meaningful)
    unrouted = run_and_score(str(p), enforce_water_balance="water-year",
                             metric="KGE", store_fluxes=True)
    mod0 = pd.to_numeric(
        unrouted.buckets.hydrodata["Specific Discharge (modeled) [mm/day]"],
        errors="coerce").to_numpy()
    both = np.isfinite(mod) & np.isfinite(mod0)
    assert not np.allclose(mod[both], mod0[both], atol=1e-6)


def test_store_fluxes_partition_exact_with_baseflow(tmp_path):
    """baseflow_Q is folded into the slow source so the partition still sums to
    the final discharge (MNiMORPH/MNiShed#36.2)."""
    p = tmp_path / "cfg.yml"
    p.write_text(yaml.safe_dump(_flat_cfg()))
    res = run_and_score(str(p), enforce_water_balance="water-year", metric="KGE",
                        store_fluxes=True, baseflow_Q=0.2)
    hd = res.buckets.hydrodata
    mod = pd.to_numeric(hd["Specific Discharge (modeled) [mm/day]"],
                        errors="coerce").to_numpy()
    parts = sum(pd.to_numeric(hd[c], errors="coerce").to_numpy() for c in FLUX)
    m = np.isfinite(mod)
    assert np.allclose(parts[m], mod[m], atol=1e-9)
