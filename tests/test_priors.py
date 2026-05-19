"""
Unit tests for suggest_priors and the Priors class.
"""

import warnings

import numpy as np
import pytest
from hydroravens import suggest_priors, Priors, BrutsaertNieber, HydrographSeparation


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
    """t_efold, recession_exponents, and initial_depths all have length n_reservoirs."""
    Q = _recession_Q()
    for n in (1, 2, 3):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pr = suggest_priors(Q, n_reservoirs=n)
        assert len(pr.recession_exponents) == n
        assert len(pr.t_efold) == n
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
    assert 'e_folding_residence_times__days:' in snippet
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
# log_t_efold_bounds structure
# ---------------------------------------------------------------------------

def test_log_t_efold_bounds_structure():
    """Non-None bounds satisfy lower < initial < upper."""
    Q = _recession_Q()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pr = suggest_priors(Q, n_reservoirs=2)
    for val in pr.log_t_efold_bounds.values():
        if val is not None:
            assert val['lower'] < val['initial'] < val['upper']
