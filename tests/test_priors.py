"""
Unit tests for suggest_priors and the Priors class.
"""

import warnings

import numpy as np
import pytest
from mnished import suggest_priors, Priors, BrutsaertNieber, HydrographSeparation


def _recession_Q(n=600, tau1=8.0, tau2=400.0, Q1=5.0, Q2=2.0):
    """Two-reservoir exponential recession; BrutsaertNieber finds many pairs in it."""
    t = np.arange(n, dtype=float)
    return Q1 * np.exp(-t / tau1) + Q2 * np.exp(-t / tau2)


# ---------------------------------------------------------------------------
# Return type and basic structure
# ---------------------------------------------------------------------------

def test_suggest_priors_returns_priors_instance():
    Q = _recession_Q()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pr = suggest_priors(Q, n_reservoirs=2)
    assert isinstance(pr, Priors)


def test_suggest_priors_attached_objects():
    """pr.bn and pr.hs are the fitted analysis objects of the right types."""
    Q = _recession_Q()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pr = suggest_priors(Q, n_reservoirs=2)
    assert isinstance(pr.bn, BrutsaertNieber)
    assert isinstance(pr.hs, HydrographSeparation)
    assert pr.n_reservoirs == 2


# ---------------------------------------------------------------------------
# recession_exponents structure (pure logic; does not depend on fit quality)
# ---------------------------------------------------------------------------

def test_recession_exponents_n1():
    """n_reservoirs=1 → [b_fast]."""
    Q = _recession_Q()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pr = suggest_priors(Q, n_reservoirs=1)
    assert len(pr.recession_exponents) == 1


def test_recession_exponents_n2():
    """n_reservoirs=2 → [b_fast, 2.203]."""
    Q = _recession_Q()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pr = suggest_priors(Q, n_reservoirs=2)
    assert len(pr.recession_exponents) == 2
    assert pr.recession_exponents[-1] == pytest.approx(2.203)


def test_recession_exponents_n3():
    """n_reservoirs=3 → [b_fast, 2.203, 1.0]."""
    Q = _recession_Q()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pr = suggest_priors(Q, n_reservoirs=3)
    assert len(pr.recession_exponents) == 3
    assert pr.recession_exponents[1] == pytest.approx(2.203)
    assert pr.recession_exponents[-1] == pytest.approx(1.0)


def test_list_lengths_match_n_reservoirs():
    """recession_coeff, recession_exponents, and initial_depths all have length n_reservoirs."""
    Q = _recession_Q()
    for n in (1, 2, 3):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pr = suggest_priors(Q, n_reservoirs=n)
        assert len(pr.recession_exponents) == n
        assert len(pr.recession_coeff) == n
        assert len(pr.initial_depths) == n


# ---------------------------------------------------------------------------
# Fallback behaviour when BrutsaertNieber cannot fit
# ---------------------------------------------------------------------------

@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_bn_failure_defaults_b_fast_to_2():
    """Flat Q has no recession pairs; BN fit fails → b_fast defaults to 2.0.

    RuntimeWarnings are suppressed: numpy's Welch/nanmean emits them when
    computing the PSD of a zero-variance series, which is expected for this
    degenerate input and handled by suggest_priors' except-Exception fallback.
    """
    Q = np.ones(100)
    with pytest.warns(UserWarning, match="BrutsaertNieber fit failed"):
        pr = suggest_priors(Q, n_reservoirs=2)
    assert pr.recession_exponents[0] == pytest.approx(2.0)


# ---------------------------------------------------------------------------
# initial_depths physical checks
# ---------------------------------------------------------------------------

def test_initial_depths_non_negative():
    """All estimable initial storage depths are non-negative."""
    Q = _recession_Q()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pr = suggest_priors(Q, n_reservoirs=2)
    assert all(h >= 0 for h in pr.initial_depths if h is not None)


# ---------------------------------------------------------------------------
# to_yaml_snippet
# ---------------------------------------------------------------------------

def test_to_yaml_snippet_contains_required_sections():
    """YAML snippet contains all required reservoir and initial-condition keys."""
    Q = _recession_Q()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pr = suggest_priors(Q, n_reservoirs=2)
    snippet = pr.to_yaml_snippet()
    assert 'reservoirs:' in snippet
    assert 'recession_coefficients:' in snippet
    assert 'recession_exponents:' in snippet
    assert 'initial_conditions:' in snippet
    assert 'water_reservoir_effective_depths__mm:' in snippet


def test_to_yaml_snippet_line_count_matches_n_reservoirs():
    """YAML snippet contains one entry per reservoir for each parameter block."""
    Q = _recession_Q()
    n = 3
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pr = suggest_priors(Q, n_reservoirs=n)
    snippet = pr.to_yaml_snippet()
    # Count lines starting with '        - ' (8-space indent = parameter list items)
    item_lines = [l for l in snippet.splitlines() if l.startswith('        - ')]
    # 5 parameter blocks × n reservoirs
    assert len(item_lines) == 5 * n


# ---------------------------------------------------------------------------
# log_recession_coeff_bounds structure
# ---------------------------------------------------------------------------

def test_log_recession_coeff_bounds_structure():
    """Non-None bounds satisfy lower < initial < upper."""
    Q = _recession_Q()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pr = suggest_priors(Q, n_reservoirs=2)
    for val in pr.log_recession_coeff_bounds.values():
        if val is not None:
            assert val['lower'] < val['initial'] < val['upper']


# ---------------------------------------------------------------------------
# leafout_GDD_from_date: green-up prior from a regional leaf-out date (#35)
# ---------------------------------------------------------------------------

import pandas as pd  # noqa: E402
from mnished import leafout_GDD_from_date  # noqa: E402


def _synthetic_forcing(years=range(2000, 2005), amp=15.0, mean=8.0):
    """Daily Tmax/Tmin from a sinusoidal annual cycle (warmest ~ day 200)."""
    dates = pd.date_range(f"{min(years)}-01-01", f"{max(years)}-12-31", freq="D")
    doy = dates.dayofyear.to_numpy()
    tmean = mean + amp * np.sin(2 * np.pi * (doy - 100) / 365.0)
    return pd.DataFrame({
        "Date": dates,
        "Maximum Temperature [C]": tmean + 4.0,
        "Minimum Temperature [C]": tmean - 4.0,
    })


def test_leafout_GDD_from_date_monotonic_in_date():
    """A later leaf-out date integrates more thermal time on the same basin."""
    df = _synthetic_forcing()
    g_early = leafout_GDD_from_date(df, 5, 1)
    g_late = leafout_GDD_from_date(df, 6, 1)
    assert g_late > g_early > 0.0


def test_leafout_GDD_from_date_matches_hand_accumulation():
    """The returned value equals an explicit per-year GDD accumulation."""
    df = _synthetic_forcing(years=range(2000, 2003))
    base = 5.0
    tmean = 0.5 * (df["Maximum Temperature [C]"] + df["Minimum Temperature [C]"])
    gdd = np.maximum(tmean - base, 0.0)
    totals = []
    for y in (2000, 2001, 2002):
        target = pd.Timestamp(year=y, month=5, day=20)
        mask = (df["Date"] >= pd.Timestamp(y, 1, 1)) & (df["Date"] <= target)
        totals.append(gdd[mask].sum())
    expected = float(np.mean(totals))
    assert np.isclose(leafout_GDD_from_date(df, 5, 20, base_temperature__C=base),
                      expected)


def test_leafout_GDD_from_date_returns_per_year():
    df = _synthetic_forcing(years=range(2000, 2004))
    mean_gdd, per_year = leafout_GDD_from_date(df, 5, 20, return_years=True)
    assert set(per_year) == {2000, 2001, 2002, 2003}
    assert np.isclose(mean_gdd, np.mean(list(per_year.values())))


def test_leafout_GDD_from_date_skips_incomplete_years():
    """A trailing year whose record stops before the leaf-out date is skipped
    (with a warning), not silently under-counted."""
    df = _synthetic_forcing(years=range(2000, 2003))
    # truncate the last year at March 1 — it cannot reach a May leaf-out
    df = df[df["Date"] <= pd.Timestamp(2002, 3, 1)]
    with pytest.warns(UserWarning, match="does not reach"):
        _, per_year = leafout_GDD_from_date(df, 5, 20, return_years=True)
    assert 2002 not in per_year and 2000 in per_year


def test_leafout_GDD_from_date_base_temperature_scales():
    """A higher base temperature yields a smaller GDD accumulation."""
    df = _synthetic_forcing()
    assert (leafout_GDD_from_date(df, 5, 20, base_temperature__C=10.0)
            < leafout_GDD_from_date(df, 5, 20, base_temperature__C=0.0))


def test_leafout_GDD_from_date_missing_column_raises():
    df = _synthetic_forcing().drop(columns=["Minimum Temperature [C]"])
    with pytest.raises(KeyError, match="Minimum Temperature"):
        leafout_GDD_from_date(df, 5, 20)


def test_leafout_GDD_from_date_no_year_reaches_date_raises():
    df = _synthetic_forcing(years=range(2000, 2001))
    df = df[df["Date"] <= pd.Timestamp(2000, 2, 1)]
    with pytest.raises(ValueError, match="no year"):
        leafout_GDD_from_date(df, 5, 20)


def test_leafout_GDD_from_date_prefers_mean_column():
    """When a 'Mean Temperature [C]' column is present, the helper uses it
    (matching the model) rather than the min/max midpoint (MNiMORPH/MNiShed#8)."""
    df = _synthetic_forcing(years=range(2000, 2003))
    midpoint = 0.5 * (df["Maximum Temperature [C]"] + df["Minimum Temperature [C]"])
    df["Mean Temperature [C]"] = midpoint + 3.0          # a true mean, 3 C warmer
    with_mean = leafout_GDD_from_date(df, 5, 20)
    without = leafout_GDD_from_date(
        df.drop(columns=["Mean Temperature [C]"]), 5, 20)
    assert with_mean > without                            # warmer mean -> more GDD
