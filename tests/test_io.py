"""
Tests for mnished.io — the input-contract spec and pre-flight validator.

The contract: the real example inputs validate clean; a config/forcing pair is
checked against the documented spec, with a column reported as an *error* when
its absence breaks MNiShed and a *warning* when MNiShed silently degrades.
"""

import os

import numpy as np
import pandas as pd
import pytest

import mnished
from mnished import io

EXAMPLES = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "examples"))
CANNON_CFG = os.path.join(EXAMPLES, "cannon_forward", "cannon_cfg.yml")
CROW_WING_CFG = os.path.join(EXAMPLES, "crow_wing", "crow_wing_config.yml")


# --------------------------------------------------------------------------
# The real examples are the ground truth: they must validate clean.
# --------------------------------------------------------------------------

@pytest.mark.parametrize("cfg", [CANNON_CFG, CROW_WING_CFG])
def test_real_examples_validate_clean(cfg):
    report = mnished.validate_inputs(cfg)
    assert report.ok, str(report)
    assert not report.errors


# --------------------------------------------------------------------------
# Config validation
# --------------------------------------------------------------------------

def _good_config(**over):
    cfg = {
        "timeseries": {"datafile": "f.csv"},
        "catchment": {"drainage_basin_area__km2": 100.0,
                      "evapotranspiration_method": "datafile",
                      "water_year_start_month": 10},
        "general": {"spin_up_cycles": 1},
        "reservoirs": {"recession_coefficients": [14, 500],
                       "exfiltration_fractions": [0.3, 1.0],
                       "maximum_effective_depths__mm": [float("inf"), float("inf")]},
        "initial_conditions": {"water_reservoir_effective_depths__mm": [15, 400]},
    }
    cfg.update(over)
    return cfg


def test_good_config_passes():
    assert mnished.validate_config(_good_config()).ok


def test_missing_required_section_and_key():
    cfg = _good_config()
    del cfg["catchment"]
    cfg["general"] = {}
    r = mnished.validate_config(cfg)
    assert not r.ok
    assert any("catchment" in e for e in r.errors)
    assert any("spin_up_cycles" in e for e in r.errors)


def test_unknown_et_method_rejected():
    cfg = _good_config()
    cfg["catchment"]["evapotranspiration_method"] = "Penman"
    r = mnished.validate_config(cfg)
    assert any("evapotranspiration_method" in e for e in r.errors)


def test_exactly_one_cascade_form():
    # both reservoirs and sub_catchments -> error
    cfg = _good_config(sub_catchments=[{"name": "a"}])
    assert any("exactly one" in e for e in mnished.validate_config(cfg).errors)
    # neither -> error
    cfg = _good_config()
    del cfg["reservoirs"]
    del cfg["initial_conditions"]
    assert any("exactly one" in e for e in mnished.validate_config(cfg).errors)


def test_reservoir_list_length_mismatch():
    cfg = _good_config()
    cfg["reservoirs"]["exfiltration_fractions"] = [0.3]      # wrong length
    assert any("same length" in e for e in mnished.validate_config(cfg).errors)


def test_initial_depths_length_mismatch():
    cfg = _good_config()
    cfg["initial_conditions"]["water_reservoir_effective_depths__mm"] = [15]
    assert any("water_reservoir_effective_depths__mm" in e
               for e in mnished.validate_config(cfg).errors)


def test_sub_catchments_form_passes():
    cfg = _good_config()
    del cfg["reservoirs"]
    del cfg["initial_conditions"]
    cfg["sub_catchments"] = [{"name": "land", "area_fraction": 1.0,
                              "reservoirs": {}}]
    assert mnished.validate_config(cfg).ok


# --------------------------------------------------------------------------
# Forcing validation
# --------------------------------------------------------------------------

def _forcing(cols, n=5):
    base = {"Date": pd.date_range("2000-01-01", periods=n)}
    for c in cols:
        base[c] = np.arange(n, dtype=float) + 1.0
    return pd.DataFrame(base)


def test_forcing_always_required_columns():
    df = _forcing([])                       # only Date
    r = mnished.validate_forcing(df)
    assert any(io.PRECIP in e for e in r.errors)
    assert any(io.DISCHARGE in e for e in r.errors)


def test_datafile_requires_et_column():
    df = _forcing([io.PRECIP, io.DISCHARGE])
    cfg = {"catchment": {"evapotranspiration_method": "datafile"}}
    assert any(io.ET in e for e in mnished.validate_forcing(df, cfg).errors)


def test_thornthwaite_requires_photoperiod():
    df = _forcing([io.PRECIP, io.DISCHARGE])
    cfg = {"catchment": {"evapotranspiration_method": "ThornthwaiteChang2019"}}
    assert any(io.PHOTOPERIOD in e for e in mnished.validate_forcing(df, cfg).errors)


def test_fdd_threshold_requires_mean_temperature():
    df = _forcing([io.PRECIP, io.DISCHARGE])
    cfg = {"catchment": {"evapotranspiration_method": "datafile"},
           "snowmelt": {"fdd_threshold": 15.0}}
    r = mnished.validate_forcing(df, cfg)
    assert any(io.TMEAN in e for e in r.errors)


def test_snowpack_missing_mean_temp_is_warning_not_error():
    df = _forcing([io.PRECIP, io.DISCHARGE, io.ET])   # no Mean Temperature
    cfg = {"catchment": {"evapotranspiration_method": "datafile"},
           "modules": {"snowpack": True, "dtr_fgi_decay": False}}
    r = mnished.validate_forcing(df, cfg)
    assert r.ok                                    # no hard error
    assert any(io.TMEAN in w for w in r.warnings)


def test_dtr_decay_missing_minmax_is_warning():
    df = _forcing([io.PRECIP, io.DISCHARGE, io.TMEAN])
    cfg = {"catchment": {"evapotranspiration_method": "datafile"},
           "modules": {"snowpack": True, "dtr_fgi_decay": True}}
    r = mnished.validate_forcing(df, cfg)
    assert any(io.TMIN in w for w in r.warnings)
    assert any(io.TMAX in w for w in r.warnings)


def test_non_daily_series_is_error():
    df = pd.DataFrame({
        "Date": pd.to_datetime(["2000-01-01", "2000-01-02", "2000-01-05"]),
        io.PRECIP: [1.0, 2.0, 3.0], io.DISCHARGE: [1.0, 2.0, 3.0]})
    assert any("daily" in e for e in mnished.validate_forcing(df).errors)


def test_all_nan_required_column_is_error():
    df = _forcing([io.PRECIP, io.DISCHARGE])
    df[io.PRECIP] = np.nan
    assert any("empty" in e for e in mnished.validate_forcing(df).errors)


def test_missing_forcing_file():
    r = mnished.validate_forcing("/no/such/file.csv")
    assert any("not found" in e for e in r.errors)


# --------------------------------------------------------------------------
# ValidationReport behaviour
# --------------------------------------------------------------------------

def test_report_raise_and_str():
    r = io.ValidationReport(errors=["boom"], warnings=["meh"])
    assert not r.ok
    with pytest.raises(ValueError, match="boom"):
        r.raise_if_errors()
    assert "ERROR" in str(r) and "warning" in str(r)
    assert io.ValidationReport().raise_if_errors().ok       # no-op when clean


def test_spec_is_exposed_for_documentation():
    names = [c.name for c in io.FORCING_COLUMNS]
    assert io.PRECIP in names and io.DISCHARGE in names
    assert "timeseries" in io.CONFIG_SECTIONS


# --------------------------------------------------------------------------
# Pre-release review fixes (#1-#3): temperature disjunction, photoperiod
# false-positive, non-numeric fdd_threshold.
# --------------------------------------------------------------------------

def test_thornthwaite_accepts_minmax_for_temperature():
    """Thornthwaite needs a mean temperature, but Min+Max suffice (the model
    derives the mean) — so a min/max-only forcing validates clean."""
    df = _forcing([io.PRECIP, io.DISCHARGE, io.PHOTOPERIOD, io.TMIN, io.TMAX])
    cfg = {"catchment": {"evapotranspiration_method": "ThornthwaiteChang2019"}}
    assert mnished.validate_forcing(df, cfg).ok


def test_thornthwaite_without_any_temperature_errors():
    df = _forcing([io.PRECIP, io.DISCHARGE, io.PHOTOPERIOD])
    cfg = {"catchment": {"evapotranspiration_method": "ThornthwaiteChang2019"}}
    r = mnished.validate_forcing(df, cfg)
    assert any(io.TMIN in e for e in r.errors)      # Chang ET needs Min and Max
    assert any(io.TMAX in e for e in r.errors)


def test_datafile_phenology_photoperiod_does_not_require_photoperiod():
    """datafile mode ignores phenology, so phenology senescence_method=photoperiod
    must NOT make Photoperiod required (was a false ERROR)."""
    df = _forcing([io.PRECIP, io.DISCHARGE, io.ET, io.TMEAN])
    cfg = {"catchment": {"evapotranspiration_method": "datafile"},
           "phenology": {"enabled": True, "senescence_method": "photoperiod"},
           "modules": {"snowpack": False, "dtr_fgi_decay": False}}
    r = mnished.validate_forcing(df, cfg)
    assert r.ok and not any(io.PHOTOPERIOD in e for e in r.errors)


def test_nonnumeric_fdd_threshold_reports_cleanly():
    cfg = _good_config()
    cfg["snowmelt"] = {"fdd_threshold": "auto"}          # not a number
    r = mnished.validate_config(cfg)
    assert any("fdd_threshold must be a number" in e for e in r.errors)


def test_inf_fdd_threshold_is_accepted():
    cfg = _good_config()
    cfg["snowmelt"] = {"fdd_threshold": float("inf")}     # the 'never frozen' default
    assert mnished.validate_config(cfg).ok


def test_thornthwaite_mean_only_is_insufficient():
    """Thornthwaite-Chang's effective temperature uses Min/Max (the diurnal
    range), so a Mean-only forcing (no Min/Max) is an error even though it has a
    mean temperature."""
    df = _forcing([io.PRECIP, io.DISCHARGE, io.PHOTOPERIOD, io.TMEAN])
    cfg = {"catchment": {"evapotranspiration_method": "ThornthwaiteChang2019"}}
    r = mnished.validate_forcing(df, cfg)
    assert any(io.TMIN in e or io.TMAX in e for e in r.errors)


def test_forcing_catalog_matches_the_validator_functions():
    """The FORCING_COLUMNS catalog and the validator functions stay in sync:
    the catalog lists exactly the known columns, and every column the functions
    can require/recommend is catalogued. Guards against the spec/enforcement
    drift of MNiMORPH/MNiShed#7."""
    catalog = {c.name for c in io.FORCING_COLUMNS}
    known = {io.DATE, io.PRECIP, io.DISCHARGE, io.TMEAN, io.TMIN, io.TMAX,
             io.PHOTOPERIOD, io.ET}
    assert catalog == known
    referenced = set()
    for et in ("datafile", "ThornthwaiteChang2019"):
        cfg = {"catchment": {"evapotranspiration_method": et},
               "modules": {"snowpack": True, "dtr_fgi_decay": True},
               "snowmelt": {"fdd_threshold": 10.0}}
        referenced |= set(io.required_forcing_columns(cfg))
        referenced |= set(io.recommended_forcing_columns(cfg))
    referenced |= {io.TMEAN, io.TMIN, io.TMAX}     # the temperature-disjunction columns
    assert referenced <= catalog                   # no rule references an uncatalogued column
