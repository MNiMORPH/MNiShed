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
import pandas as pd
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


def test_dormant_Kc_default_is_half(tmp_path):
    """The Crow-Wing-tuned default avoids winter over-production."""
    cfg = _cfg(et="ThornthwaiteChang2019", phenology={"enabled": True})
    p = tmp_path / "cfg.yml"
    p.write_text(yaml.safe_dump(cfg))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Buckets()
        b.initialize(str(p), enforce_water_balance="global")
    assert b.phenology_params["dormant_Kc"] == 0.5


def test_leafout_GDD_is_a_calibration_target(tmp_path):
    """leafout_GDD is overridable per evaluation (the one calibratable knob),
    and ScoringModel's build-once path matches run_and_score bit-for-bit — so
    the Thornthwaite cache / Kc-recompute split is correct."""
    from mnished import ScoringModel
    cfg = _cfg(et="ThornthwaiteChang2019", phenology={"enabled": True})
    path = _write(tmp_path, cfg, "ph")
    early = mnished.run_and_score(path, enforce_water_balance="global",
                                  metric="KGE", leafout_GDD=80).score
    late = mnished.run_and_score(path, enforce_water_balance="global",
                                 metric="KGE", leafout_GDD=320).score
    assert early != late                                  # the knob moves the fit
    sm = ScoringModel(path, enforce_water_balance="global")
    assert sm.score(metric="KGE", leafout_GDD=80).score == early
    assert sm.score(metric="KGE", leafout_GDD=320).score == late


def test_leafout_GDD_ignored_without_phenology(tmp_path):
    """leafout_GDD with phenology disabled warns and has no effect."""
    path = _write(tmp_path, _cfg(et="ThornthwaiteChang2019"), "noph")
    base = mnished.run_and_score(path, enforce_water_balance="global",
                                 metric="KGE").score
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        got = mnished.run_and_score(path, enforce_water_balance="global",
                                    metric="KGE", leafout_GDD=200).score
    assert got == base
    assert any("leafout_GDD" in str(w.message) for w in caught)


def test_phenology_with_et_reservoir_draw(tmp_path):
    """Phenology composes with et_reservoir_draw + none/et_scale closure — the
    mode the Crow Wing calibration uses — and leafout_GDD still moves the fit."""
    cfg = _cfg(et="ThornthwaiteChang2019", phenology={"enabled": True})
    cfg["modules"]["et_reservoir_draw"] = True
    et, _ = _et_series(tmp_path, cfg, wb="none")
    assert np.isfinite(et).any()
    path = _write(tmp_path, cfg, "rdraw")
    a = mnished.run_and_score(path, enforce_water_balance="none", metric="KGE",
                              leafout_GDD=70).score
    b = mnished.run_and_score(path, enforce_water_balance="none", metric="KGE",
                              leafout_GDD=260).score
    assert a != b


# --------------------------------------------------------------------------
# Photoperiod-driven autumn senescence (issue #35, part 1)
# --------------------------------------------------------------------------

def _Kc(tmp_path, phenology, name="ph"):
    """Phenology Kc array and its day-of-year for a built Buckets."""
    path = _write(tmp_path, _cfg(phenology=phenology), name)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Buckets()
        b.initialize(path, enforce_water_balance="none")
    Kc = np.asarray(b.phenology_Kc())
    doy = b.hydrodata["Date"].dt.dayofyear.to_numpy()
    return Kc, doy


def test_senescence_method_defaults_to_doy(tmp_path):
    """Default senescence_method is 'doy' — backward-compatible with v3.2.0."""
    path = _write(tmp_path, _cfg(phenology={"enabled": True}), "def")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Buckets()
        b.initialize(path, enforce_water_balance="none")
    assert b.phenology_params["senescence_method"] == "doy"


def test_photoperiod_and_doy_agree_before_autumn(tmp_path):
    """The two senescence forms differ only in autumn: identical green-up and
    summer plateau (senescence inactive before the photoperiod/DOY trigger)."""
    kc_doy, doy = _Kc(tmp_path, {"enabled": True}, "doy")
    kc_pp, _ = _Kc(tmp_path, {"enabled": True,
                              "senescence_method": "photoperiod"}, "pp")
    pre_autumn = doy < 230          # green-up + summer, before either trigger
    assert np.allclose(kc_doy[pre_autumn], kc_pp[pre_autumn])
    # they do diverge somewhere in the autumn brown-down
    assert not np.allclose(kc_doy, kc_pp)


def test_photoperiod_no_spring_browndown(tmp_path):
    """The post-solstice gate keeps spring's equally short days from triggering
    senescence: Kc still leafs out to full canopy by midsummer."""
    kc, doy = _Kc(tmp_path, {"enabled": True,
                             "senescence_method": "photoperiod"})
    # full canopy reached in July despite short photoperiods earlier in spring
    jul = (doy >= 182) & (doy <= 212)
    assert np.isclose(kc[jul].max(), 1.0, atol=1e-6)
    # spring rises April -> June (leaf-out, not senescing)
    apr = (doy >= 91) & (doy <= 120)
    jun = (doy >= 152) & (doy <= 181)
    assert kc[apr].mean() < kc[jun].mean()


def test_photoperiod_senesces_to_dormant_in_late_autumn(tmp_path):
    """As day length falls well below the critical value, Kc returns toward
    dormant_Kc — the autumn ET draw-down the freshet/fall lever depends on."""
    kc, doy = _Kc(tmp_path, {"enabled": True, "dormant_Kc": 0.5,
                             "senescence_method": "photoperiod"})
    nov = (doy >= 320) & (doy <= 350)        # short days, fully senesced
    assert np.isclose(kc[nov].min(), 0.5, atol=1e-6)


def test_bad_senescence_method_raises(tmp_path):
    """An unknown senescence_method is rejected at build time."""
    path = _write(tmp_path, _cfg(phenology={"enabled": True,
                                            "senescence_method": "bogus"}), "bad")
    with pytest.raises(ValueError, match="senescence_method"):
        b = Buckets()
        b.initialize(path, enforce_water_balance="none")


def test_photoperiod_method_needs_photoperiod_column(tmp_path):
    """Photoperiod senescence raises a clear error when the forcing lacks the
    'Photoperiod [hr]' column (rather than failing obscurely downstream)."""
    df = pd.read_csv(CANNON_CSV).drop(columns=["Photoperiod [hr]"])
    csv = tmp_path / "no_photoperiod.csv"
    df.to_csv(csv, index=False)
    cfg = _cfg(phenology={"enabled": True, "senescence_method": "photoperiod"})
    cfg["timeseries"]["datafile"] = str(csv)
    path = _write(tmp_path, cfg, "nopp")
    with pytest.raises(ValueError, match="Photoperiod"):
        b = Buckets()
        b.initialize(path, enforce_water_balance="none")


# --------------------------------------------------------------------------
# Photoperiod-senescence latitude TRANSFER (issue #35) — a deterministic
# property of the mechanism, testable without any calibration: feed identical
# forcings differing only in latitude (the Photoperiod column) and check where
# the autumn brown-down lands.
# --------------------------------------------------------------------------

def _day_length_hours(latitude_deg, doy):
    """Astronomical day length [hr] (same model as the docs senescence figure)."""
    decl = np.radians(23.44) * np.sin(2 * np.pi * (doy - 81) / 365.0)
    h = np.clip(-np.tan(np.radians(latitude_deg)) * np.tan(decl), -1.0, 1.0)
    return (24.0 / np.pi) * np.arccos(h)


def _synthetic_forcing(tmp_path, latitude, name):
    """A multi-year daily forcing whose ONLY latitude dependence is Photoperiod;
    temperature is a fixed warm-summer sinusoid so green-up is identical."""
    dates = pd.date_range("2000-01-01", "2003-12-31", freq="D")
    doy = dates.dayofyear.to_numpy()
    tmean = 10.0 + 18.0 * np.sin(2 * np.pi * (doy - 110) / 365.0)   # warm summer
    df = pd.DataFrame({
        "Date": dates,
        "Precipitation [mm/day]": 2.0,
        "Discharge [m^3/s]": 5.0,
        "Mean Temperature [C]": tmean,
        "Minimum Temperature [C]": tmean - 5.0,
        "Maximum Temperature [C]": tmean + 5.0,
        "Photoperiod [hr]": _day_length_hours(latitude, doy),
    })
    p = tmp_path / f"forcing_{name}.csv"
    df.to_csv(p, index=False)
    return str(p)


def _browndown_doy(tmp_path, latitude, phenology, name, level="midpoint"):
    """Day-of-year (year 2001) where autumn Kc falls to a target level.

    ``level='midpoint'`` uses the full/dormant midpoint (most latitude spread
    when the critical day length is well below the equinox); ``level='onset'``
    uses the first decline below full canopy (Kc crosses near N = the critical
    day length), which resolves the transfer direction even near the equinox.
    """
    cfg = _cfg(et="ThornthwaiteChang2019", phenology=phenology)
    cfg["timeseries"]["datafile"] = _synthetic_forcing(tmp_path, latitude, name)
    path = _write(tmp_path, cfg, name)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Buckets()
        b.initialize(path, enforce_water_balance="none")
    Kc = np.asarray(b.phenology_Kc())
    dates = b.hydrodata["Date"]
    doy = dates.dt.dayofyear.to_numpy()
    yr = (dates.dt.year == 2001).to_numpy()
    p = b.phenology_params
    if level == "onset":
        target = p["full_Kc"] - 0.01 * (p["full_Kc"] - p["dormant_Kc"])
    else:
        target = 0.5 * (p["full_Kc"] + p["dormant_Kc"])
    autumn = yr & (doy >= 200)
    d, k = doy[autumn], Kc[autumn]
    order = np.argsort(d)
    d, k = d[order], k[order]
    below = np.where(k <= target)[0]
    return int(d[below[0]]) if len(below) else None


def test_photoperiod_senescence_transfers_north_earlier(tmp_path):
    """The headline claim: with a SUB-equinox critical day length (the 11.0-hr
    default), a higher-latitude basin browns down EARLIER — and the same two
    numbers do it, with no per-basin re-tuning."""
    phen = {"enabled": True, "senescence_method": "photoperiod"}   # defaults 11.0/1.8
    north = _browndown_doy(tmp_path, 47.5, phen, "n")
    south = _browndown_doy(tmp_path, 41.5, phen, "s")
    assert north is not None and south is not None
    assert north < south - 2          # north meaningfully earlier (not just noise)


def test_doy_senescence_does_not_transfer(tmp_path):
    """The contrast: the fixed calendar window places brown-down on the same date
    regardless of latitude (zero transfer)."""
    phen = {"enabled": True, "senescence_method": "doy"}
    north = _browndown_doy(tmp_path, 47.5, phen, "dn")
    south = _browndown_doy(tmp_path, 41.5, phen, "ds")
    assert north == south


def test_photoperiod_transfer_reverses_above_the_equinox(tmp_path):
    """The equinox-pivot finding, locked in: a critical day length ABOVE the
    ~12-hr equinox is crossed BEFORE the equinox, where higher-latitude days are
    still longer — so the transfer reverses (north later). This is why the
    default sits below the equinox."""
    below = {"enabled": True, "senescence_method": "photoperiod",
             "senescence_photoperiod__hr": 11.0}
    above = {"enabled": True, "senescence_method": "photoperiod",
             "senescence_photoperiod__hr": 13.0}
    n_below = _browndown_doy(tmp_path, 47.5, below, "nb", level="onset")
    s_below = _browndown_doy(tmp_path, 41.5, below, "sb", level="onset")
    n_above = _browndown_doy(tmp_path, 47.5, above, "na", level="onset")
    s_above = _browndown_doy(tmp_path, 41.5, above, "sa", level="onset")
    assert n_below < s_below          # sub-equinox: north earlier (correct)
    assert n_above > s_above          # above-equinox: north later (reversed)


def _kc_from_df(tmp_path, df, name):
    """Build a Thornthwaite+phenology model from a forcing DataFrame; return
    (buckets, Kc array)."""
    csv = tmp_path / f"{name}.csv"
    df.to_csv(csv, index=False)
    cfg = _cfg(phenology={"enabled": True})
    cfg["timeseries"]["datafile"] = str(csv)
    path = _write(tmp_path, cfg, name)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Buckets()
        b.initialize(path, enforce_water_balance="global")
    return b, np.asarray(b.phenology_Kc())


def test_phenology_falls_back_on_minmax_mean(tmp_path):
    """phenology GDD uses 'Mean Temperature [C]'; with that column absent the
    model synthesizes it from (Tmin+Tmax)/2 at load, so phenology runs from a
    min/max-only record and gives the SAME Kc as when the (midpoint) Mean column
    is present. Confirms the fallback (Cannon's Mean already equals the
    midpoint, so dropping it must not change Kc)."""
    df = pd.read_csv(CANNON_CSV, parse_dates=["Date"])
    b_with, kc_with = _kc_from_df(tmp_path, df, "with_mean")
    b_no, kc_no = _kc_from_df(
        tmp_path, df.drop(columns=["Mean Temperature [C]"]), "no_mean")
    assert "Mean Temperature [C]" in b_no.hydrodata.columns   # synthesized at load
    assert np.allclose(kc_with, kc_no)                        # identical GDD -> Kc
    assert np.isfinite(kc_no).all()
