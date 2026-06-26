"""
Unit tests for mnished.identifiability — post-fit parameter-identifiability
diagnostics.

All tests use cheap analytic objectives (no model runs), so the suite is
fast and deterministic.  Three reference objectives are used throughout:

* ``quad`` — a separable objective with one sharply-constrained parameter,
  one nearly-flat parameter, and one whose optimum sits on a bound.
* ``ridge_obj`` — a deliberately degenerate objective, sharp along the
  combination ``x0 - x1`` and flat (sloppy) along ``x0 + x1``.
* ``tilted_opt`` — used to check that a supplied point which is *not* a
  true maximum is reported honestly (positive curvature).
"""

import matplotlib
matplotlib.use("Agg")  # headless: plotting must not require a display

import numpy as np
import pytest

from mnished.identifiability import (
    Parameter,
    ParameterSet,
    eigenspectrum,
    profile,
    profile_all,
    ridge,
    _CachedObjective,
)


# --- reference analytic objectives --------------------------------------

def quad(t):
    """Sharp in x0, ~flat in x1, optimum of x2 at its upper bound (2)."""
    return -(t["x0"]) ** 2 + 0.001 * t["x1"] - 0.5 * (t["x2"] - 2.0) ** 2


def ridge_obj(t):
    """Stiff along (x0 - x1); sloppy along (x0 + x1)."""
    return -10.0 * (t["x0"] - t["x1"]) ** 2 - 0.01 * (t["x0"] + t["x1"]) ** 2


def quad_pset():
    return ParameterSet([
        Parameter("x0", value=0.0, lower=-1, upper=1, log=False),
        Parameter("x1", value=0.0, lower=-1, upper=1, log=False),
        Parameter("x2", value=2.0, lower=0,  upper=2, log=False),
    ])


def ridge_pset():
    return ParameterSet([
        Parameter("x0", value=0.0, lower=-1, upper=1, log=False),
        Parameter("x1", value=0.0, lower=-1, upper=1, log=False),
    ])


# --- Parameter / ParameterSet -------------------------------------------

def test_parameter_infers_log_from_name():
    assert Parameter("log__tau", 1.0, 0, 2).log is True
    assert Parameter("f_frac", 0.5, 0, 1).log is False
    # explicit flag overrides the name-based inference
    assert Parameter("log__x", 1.0, 0, 2, log=False).log is False


def test_parameter_rejects_bad_bounds():
    with pytest.raises(ValueError):
        Parameter("x", 0.0, lower=1.0, upper=1.0)


def test_parameter_normalize():
    p = Parameter("x", 0.0, lower=-2, upper=2)
    assert p.normalize(0.0) == pytest.approx(0.5)
    assert p.normalize(-2.0) == pytest.approx(0.0)


def test_parameterset_rejects_duplicates():
    with pytest.raises(ValueError):
        ParameterSet([Parameter("x", 0, 0, 1), Parameter("x", 0, 0, 1)])


def test_from_params_yml_filters_fixed_and_uses_optimum():
    cfg = {
        "log__a": {"lower": 0.0, "upper": 2.0, "initial": 1.0, "active": True},
        "b":      {"lower": 0.0, "upper": 1.0, "initial": 0.5, "active": True},
        "c_off":  {"lower": 0.0, "upper": 1.0, "initial": 0.5, "active": False},
    }
    pset = ParameterSet.from_params_yml(cfg, optimum={"log__a": 1.7})
    assert pset.names == ["log__a", "b"]           # fixed 'c_off' dropped
    assert pset["log__a"].value == pytest.approx(1.7)   # optimum override
    assert pset["b"].value == pytest.approx(0.5)        # initial fallback


# --- profiles -----------------------------------------------------------

def test_profile_sharp_parameter_is_constrained():
    pr = profile(quad, quad_pset(), "x0", n=21)
    assert pr.curvature() < 0            # concave at a real maximum
    assert pr.half_width() < 0.2         # narrow → well constrained
    assert pr.at_bound() is False


def test_profile_flat_parameter_is_unconstrained():
    pr = profile(quad, quad_pset(), "x1", n=21)
    assert pr.half_width() > 0.4         # spans most of the range
    assert pr.at_bound() is False        # flat-but-tilted is NOT pegged


def test_profile_detects_bound_pegging():
    pr = profile(quad, quad_pset(), "x2", n=21)
    assert pr.at_bound() is True


def test_profile_reports_non_maximum_honestly():
    """A held point that isn't the optimum: the profile reveals a better value."""
    # hold x2 at 1.0 (interior) though the true optimum is at the bound (2.0)
    pset = ParameterSet([
        Parameter("x0", 0.0, -1, 1), Parameter("x1", 0.0, -1, 1),
        Parameter("x2", 1.0, 0, 2),
    ])
    pr = profile(quad, pset, "x2", n=21)
    best = int(np.nanargmax(pr.score))
    assert pr.x[best] > pr.x_opt          # a better x2 exists above the held point
    assert pr.score[best] > pr.score_opt  # and it scores higher


def test_profile_carries_model_failures_as_nan():
    def flaky(t):
        return np.nan if t["x0"] > 0.5 else -(t["x0"]) ** 2
    pset = ParameterSet([Parameter("x0", 0.0, -1, 1)])
    pr = profile(flaky, pset, "x0", n=21)
    assert np.isnan(pr.score).any()
    assert np.isfinite(pr.half_width())  # robust to the nan tail


# --- eigenspectrum ------------------------------------------------------

def test_eigenspectrum_recovers_stiff_and_sloppy_axes():
    sp = eigenspectrum(ridge_obj, ridge_pset(), step=0.05)
    stiff, sloppy = sp.eigvecs[:, 0], sp.eigvecs[:, 1]
    # stiff direction is (1, -1); sloppy is (1, 1) (up to sign)
    assert stiff[0] * stiff[1] < 0
    assert sloppy[0] * sloppy[1] > 0
    assert abs(abs(stiff[0]) - abs(stiff[1])) < 0.05


def test_eigenspectrum_condition_number_large_for_degenerate():
    sp = eigenspectrum(ridge_obj, ridge_pset(), step=0.05)
    assert sp.condition_number > 100


def test_spectrum_excludes_pegged_parameter():
    """A parameter at a bound is dropped, not allowed to void the spectrum."""
    # x2 sits exactly at its upper bound -> no FD room -> excluded
    sp = eigenspectrum(quad, quad_pset(), step=0.05)
    assert "x2" in sp.excluded
    assert sp.names == ["x0", "x1"]
    assert np.all(np.isfinite(sp.eigvals))   # remaining spectrum is finite
    assert sp.eigvecs.shape == (2, 2)


def test_flat_directions_names_the_combination():
    sp = eigenspectrum(ridge_obj, ridge_pset(), step=0.05)
    flat = sp.flat_directions(rel_tol=1e-2)
    assert len(flat) == 1
    _, loading = flat[0]
    # both parameters load with comparable magnitude and the same sign
    a, b = loading["x0"], loading["x1"]
    assert a * b > 0 and abs(abs(a) - abs(b)) < 0.05


# --- 2-D ridge ----------------------------------------------------------

def test_ridge_axis_recovers_degenerate_direction():
    rg = ridge(ridge_obj, ridge_pset(), "x0", "x1", n=21)
    v = rg.ridge_axis(delta=0.01)
    assert v is not None
    assert abs(abs(v[0]) - 0.707) < 0.1 and v[0] * v[1] > 0


def test_ridge_axis_matches_sloppy_eigenvector():
    """The Hessian names the ridge; the grid confirms it — they agree."""
    pset = ridge_pset()
    sp = eigenspectrum(ridge_obj, pset, step=0.05)
    sloppy = sp.eigvecs[:, 1]
    v = np.array(ridge(ridge_obj, pset, "x0", "x1", n=21).ridge_axis())
    # compare as undirected axes: |cos angle| ~ 1
    cos = abs(v @ sloppy) / (np.linalg.norm(v) * np.linalg.norm(sloppy))
    assert cos > 0.95


# --- report -------------------------------------------------------------

def test_report_summary_and_seam():
    rep = profile_all(quad, quad_pset(), n=21, aic=-123.4)
    rep.spectrum = eigenspectrum(quad, quad_pset(), step=0.05)
    text = rep.summary()
    assert "x0" in text and "x1" in text and "x2" in text
    assert "AIC" in text or "-123" in text
    assert "eigenspectrum" in text
    # the a-priori reconciliation seam is empty until the forward tools exist
    assert rep.predicted is None
    assert "not available" in text


# --- caching ------------------------------------------------------------

def test_cached_objective_deduplicates():
    calls = {"n": 0}

    def counting(t):
        calls["n"] += 1
        return -(t["x0"]) ** 2
    cobj = _CachedObjective(counting, ["x0"])
    cobj({"x0": 0.3})
    cobj({"x0": 0.3})   # identical → served from cache
    assert calls["n"] == 1
    assert cobj.n_calls == 1


# --- plotting smoke tests (Agg backend) ---------------------------------

def test_plot_methods_do_not_raise():
    import matplotlib.pyplot as plt
    rep = profile_all(quad, quad_pset(), n=11)
    rep.spectrum = eigenspectrum(quad, quad_pset(), step=0.05)
    rep.profiles["x0"].plot()
    rep.spectrum.plot()
    ridge(quad, quad_pset(), "x0", "x1", n=7).plot()
    rep.plot()
    plt.close("all")
