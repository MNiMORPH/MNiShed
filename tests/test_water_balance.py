"""
Water-balance closure tests for the ET scaling.

``enforce_water_balance`` is supposed to scale ET so that P - Q - ET = 0 over
each water year ('water-year') or over the full record ('global'). The scaling
must normalise against the *same* ET demand that is actually applied — measured
ET in datafile mode, Thornthwaite × phenology in Thornthwaite mode. These tests
assert that invariant directly for every et_method × mode (with and without
phenology); its absence let a regression slip in where the water-year multiplier
divided by the raw input ET column instead of the applied Thornthwaite demand,
so Thornthwaite + water-year did not close (off by ~2×).
"""

import os
import warnings

import numpy as np
import pytest
import yaml

from mnished import Buckets

EXAMPLE_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "examples", "cannon_forward")
)
CANNON_CSV = os.path.join(EXAMPLE_DIR, "CannonTestInput.csv")


def _cfg(et, phenology=False):
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
    if phenology:
        cfg["phenology"] = {"enabled": True}
    return cfg


def _build(tmp_path, cfg, wb):
    p = tmp_path / "cfg.yml"
    p.write_text(yaml.safe_dump(cfg))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Buckets()
        b.initialize(str(p), enforce_water_balance=wb)
        b.run()
    return b


def _terms(b):
    et = b.hydrodata["ET for model [mm/day]"].to_numpy()
    P = b.hydrodata["Precipitation [mm/day]"].to_numpy()
    Q = b.hydrodata["Specific Discharge [mm/day]"].to_numpy()
    wy = b.hydrodata["Water Year"].to_numpy()
    finite = np.isfinite(P) & np.isfinite(Q)
    return et, P, Q, wy, finite


@pytest.mark.parametrize("et", ["datafile", "ThornthwaiteChang2019"])
@pytest.mark.parametrize("phenology", [False, True])
def test_water_year_closure(tmp_path, et, phenology):
    """Each water year's scaled ET equals its P - Q (the documented contract)."""
    if phenology and et == "datafile":
        pytest.skip("phenology is ignored for measured ET")
    b = _build(tmp_path, _cfg(et, phenology), wb="water-year")
    et_mod, P, Q, wy, finite = _terms(b)
    for w in np.unique(wy):
        m = (wy == w) & finite
        if m.sum() < 300:                      # skip partial first/last years
            continue
        assert np.isclose(np.nansum(et_mod[m]), np.nansum((P - Q)[m]), rtol=1e-6)


@pytest.mark.parametrize("et", ["datafile", "ThornthwaiteChang2019"])
@pytest.mark.parametrize("phenology", [False, True])
def test_global_closure(tmp_path, et, phenology):
    """Global mode closes the full record (not each year individually)."""
    if phenology and et == "datafile":
        pytest.skip("phenology is ignored for measured ET")
    b = _build(tmp_path, _cfg(et, phenology), wb="global")
    et_mod, P, Q, _, finite = _terms(b)
    assert np.isclose(np.nansum(et_mod[finite]), np.nansum((P - Q)[finite]),
                      rtol=1e-6)


def test_datafile_demand_is_the_input_column(tmp_path):
    """In datafile mode the demand == the input ET column, so routing the
    multipliers through _demand_ET() leaves datafile results bit-identical."""
    b = _build(tmp_path, _cfg("datafile"), wb="water-year")
    assert np.allclose(b._demand_ET(),
                       b.hydrodata["Evapotranspiration [mm/day]"].to_numpy(),
                       equal_nan=True)


def test_water_year_closure_under_ragged_gaps(tmp_path):
    """With P, Q, and the ET demand each missing on *different* days, the
    per-water-year closure must still hold over the days where all three are
    present — the multiplier is built on that common finite mask, not on
    per-column NaN-skipping means (which would close the balance only
    approximately under ragged gaps). Regression for MNiMORPH/MNiShed#36.1."""
    import pandas as pd
    df = pd.read_csv(CANNON_CSV, parse_dates=["Date"])
    i = np.arange(len(df))
    df.loc[i % 9 == 0, "Discharge [m^3/s]"] = np.nan       # Q gaps ...
    df.loc[i % 7 == 0, "Precipitation [mm/day]"] = np.nan  # ... on different days than P
    ragged = tmp_path / "ragged.csv"
    df.to_csv(ragged, index=False)

    cfg = _cfg("datafile")
    cfg["timeseries"]["datafile"] = str(ragged)
    b = _build(tmp_path, cfg, "water-year")

    et, P, Q, wy, _ = _terms(b)
    demand = np.asarray(b._demand_ET(), dtype=float)
    common = np.isfinite(P) & np.isfinite(Q) & np.isfinite(demand)
    for y in np.unique(wy[common]):
        m = common & (wy == y)
        if m.sum() < 5:
            continue
        residual = P[m].mean() - Q[m].mean() - et[m].mean()
        assert abs(residual) < 1e-9, f"WY {y} does not close: residual={residual}"


def test_water_year_multiplier_no_inf_on_zero_demand(tmp_path):
    """A water year whose finite days carry zero ET demand divides by zero; the
    multiplier is mapped to NaN (then handled as raw ET) rather than propagating
    inf into the modelled discharge. Regression for the #36 pre-release review."""
    import pandas as pd
    df = pd.read_csv(CANNON_CSV, parse_dates=["Date"])
    wy93 = (df["Date"] >= "1992-10-01") & (df["Date"] <= "1993-09-30")
    df.loc[wy93, "Evapotranspiration [mm/day]"] = 0.0        # zero demand that WY
    zero = tmp_path / "zero_et.csv"
    df.to_csv(zero, index=False)
    cfg = _cfg("datafile")
    cfg["timeseries"]["datafile"] = str(zero)
    b = _build(tmp_path, cfg, "water-year")
    mult = b.hydrodata_WY_means["ET multiplier"].to_numpy()
    assert not np.isinf(mult).any()                         # guarded: no infinities
    et = b.hydrodata["ET for model [mm/day]"].to_numpy()
    assert not np.isinf(et).any()                           # no inf propagation
