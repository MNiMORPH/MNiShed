"""
Unit tests for HydrographSeparation.

Most tests supply tau_fast and tau_deep directly to bypass spectral and
recession fitting; this makes them fast and deterministic.  One smoke test
runs the full pipeline on synthetic data.
"""

import warnings

import numpy as np
import pytest
from mnished import HydrographSeparation


def _two_res_Q(n=600, tau1=8.0, tau2=400.0, Q1=5.0, Q2=2.0):
    """Pure two-reservoir exponential recession Q(t) = Q1·e^{-t/τ1} + Q2·e^{-t/τ2}."""
    t = np.arange(n, dtype=float)
    return Q1 * np.exp(-t / tau1) + Q2 * np.exp(-t / tau2)


# ---------------------------------------------------------------------------
# Low-pass filter
# ---------------------------------------------------------------------------

def test_apply_lowpass_step_response():
    """Filter converges to 1 after a step input of sufficient length."""
    # Step from 0 to 1 at t=100; filter with tau=5 days.
    # After 5*tau = 25 steps past the step the filter is >99% converged.
    signal = np.concatenate([np.zeros(100), np.ones(200)])
    hs = HydrographSeparation(signal, tau_fast=[5.0], tau_deep=100.0)
    result = hs._apply_lowpass(signal, tau=5.0)
    assert result[0] == pytest.approx(0.0)
    assert result[-1] == pytest.approx(1.0, abs=0.01)


def test_apply_lowpass_constant_input():
    """Filter is a fixed point for a constant input."""
    signal = np.full(100, 3.7)
    hs = HydrographSeparation(signal, tau_fast=[5.0], tau_deep=100.0)
    result = hs._apply_lowpass(signal, tau=10.0)
    np.testing.assert_allclose(result, 3.7)


def test_apply_lowpass_large_tau_is_identity():
    """Very large tau → filter barely moves; output ≈ first value everywhere."""
    signal = np.random.default_rng(0).uniform(1, 10, 200)
    hs = HydrographSeparation(signal, tau_fast=[5.0], tau_deep=100.0)
    result = hs._apply_lowpass(signal, tau=1e8)
    assert result[-1] == pytest.approx(signal[0], rel=1e-3)


# ---------------------------------------------------------------------------
# Bypass constructors (tau_fast / tau_deep supplied)
# ---------------------------------------------------------------------------

def test_tau_deep_bypass_sets_tau_karst():
    """Supplying tau_deep skips recession fitting and uses the value directly."""
    Q = _two_res_Q()
    hs = HydrographSeparation(Q, tau_fast=[8.0], tau_deep=400.0)
    hs.fit()
    assert hs.tau_karst == pytest.approx(400.0)


def test_tau_fast_bypass_skips_spectral():
    """Supplying tau_fast sets aic_scores to None (spectral fitting skipped)."""
    Q = _two_res_Q()
    hs = HydrographSeparation(Q, tau_fast=[8.0], tau_deep=400.0)
    hs.fit()
    assert hs.aic_scores is None


def test_n_reservoirs_fitted_matches_tau_fast_length():
    """n_reservoirs_fitted = len(tau_fast) + 1 (for the deep reservoir)."""
    Q = _two_res_Q()
    hs = HydrographSeparation(Q, tau_fast=[5.0, 30.0], tau_deep=400.0)
    hs.fit()
    assert hs.n_reservoirs_fitted == 3


# ---------------------------------------------------------------------------
# fit() return value and state
# ---------------------------------------------------------------------------

def test_fit_returns_self():
    """fit() supports method chaining."""
    Q = _two_res_Q()
    hs = HydrographSeparation(Q, tau_fast=[8.0], tau_deep=400.0)
    assert hs.fit() is hs


def test_tau_and_h0_set_after_fit():
    """tau and H0 arrays are populated after fit()."""
    Q = _two_res_Q()
    hs = HydrographSeparation(Q, tau_fast=[8.0], tau_deep=400.0)
    hs.fit()
    assert hs.tau is not None
    assert hs.H0 is not None
    assert len(hs.tau) == hs.n_reservoirs_fitted
    assert len(hs.H0) == hs.n_reservoirs_fitted


def test_h0_non_negative():
    """All initial storage depths are non-negative."""
    Q = _two_res_Q()
    hs = HydrographSeparation(Q, tau_fast=[8.0], tau_deep=400.0)
    hs.fit()
    assert (hs.H0 >= 0).all()


def test_h0_karst_formula():
    """H0_karst (slowest, index 0) equals Q_residual[0] * tau_karst."""
    Q = _two_res_Q()
    hs = HydrographSeparation(Q, tau_fast=[8.0], tau_deep=400.0)
    hs.fit()
    expected = float(hs._Q_residual[0]) * hs.tau_karst
    assert hs.H0[0] == pytest.approx(expected, rel=1e-9)


# ---------------------------------------------------------------------------
# get_initial_conditions
# ---------------------------------------------------------------------------

def test_get_initial_conditions_structure():
    """Returns a dict with 'H0' key; length matches n_reservoirs."""
    Q = _two_res_Q()
    hs = HydrographSeparation(Q, n_reservoirs=2, tau_fast=[8.0], tau_deep=400.0)
    hs.fit()
    ic = hs.get_initial_conditions()
    assert 'H0' in ic
    assert len(ic['H0']) == 2


def test_get_initial_conditions_order():
    """get_initial_conditions() returns fastest-first (reversed from internal)."""
    Q = _two_res_Q()
    hs = HydrographSeparation(Q, tau_fast=[8.0], tau_deep=400.0)
    hs.fit()
    ic = hs.get_initial_conditions()
    # Internal self.H0[0] is karst (slowest); ic['H0'][0] is fast (shallowest).
    # Fastest reservoir has shorter tau and lower H0 for this synthetic series.
    assert ic['H0'][0] <= ic['H0'][-1]


def test_get_initial_conditions_raises_before_fit():
    """get_initial_conditions() raises RuntimeError if fit() not called."""
    hs = HydrographSeparation(np.ones(100), tau_fast=[5.0], tau_deep=100.0)
    with pytest.raises(RuntimeError):
        hs.get_initial_conditions()


# ---------------------------------------------------------------------------
# get_parameter_priors
# ---------------------------------------------------------------------------

def test_get_parameter_priors_structure():
    """Each non-None prior has initial, lower, upper with lower < initial < upper."""
    Q = _two_res_Q()
    hs = HydrographSeparation(Q, n_reservoirs=2, tau_fast=[8.0], tau_deep=400.0)
    hs.fit()
    priors = hs.get_parameter_priors()
    for val in priors.values():
        if val is not None:
            assert 'initial' in val and 'lower' in val and 'upper' in val
            assert val['lower'] < val['initial'] < val['upper']


def test_get_parameter_priors_raises_before_fit():
    """get_parameter_priors() raises RuntimeError if fit() not called."""
    hs = HydrographSeparation(np.ones(100), tau_fast=[5.0], tau_deep=100.0)
    with pytest.raises(RuntimeError):
        hs.get_parameter_priors()


def test_get_parameter_priors_log_units():
    """Prior bounds are in log10(days): initial ≈ log10(tau_deep) for the deep slot."""
    Q = _two_res_Q()
    tau_deep = 400.0
    hs = HydrographSeparation(Q, n_reservoirs=2, tau_fast=[8.0], tau_deep=tau_deep)
    hs.fit()
    priors = hs.get_parameter_priors()
    # The last (karst) key corresponds to tau_deep = 400 d → log10(400) ≈ 2.602
    last_key = list(priors.keys())[-1]
    if priors[last_key] is not None:
        assert priors[last_key]['initial'] == pytest.approx(np.log10(tau_deep), abs=0.01)


# ---------------------------------------------------------------------------
# Full pipeline smoke test (spectral + recession fitting, no bypasses)
# ---------------------------------------------------------------------------

def test_full_pipeline_smoke():
    """fit() completes on a long synthetic recession without raising."""
    Q = _two_res_Q(n=600)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        hs = HydrographSeparation(Q, n_reservoirs=2).fit()
    assert hs.tau is not None
    assert hs.tau_karst > 0
    assert (hs.H0 >= 0).all()
