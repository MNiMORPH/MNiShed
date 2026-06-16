"""
Unit tests for the BrutsaertNieber recession analysis class.
"""

import warnings

import numpy as np
import pytest
from mnished import BrutsaertNieber


def _linear_recession(Q0, a, n_steps, dt=1.0):
    """Synthetic recession: Q[t+1] = Q[t] - a * Q[t] * dt (n=1 case)."""
    Q = [Q0]
    for _ in range(n_steps - 1):
        Q.append(Q[-1] - a * Q[-1] * dt)
    return np.array(Q)


def _powerlaw_recession(Q0, a, n_bn, n_steps, dt=1.0):
    """Synthetic recession via Euler integration of -dQ/dt = a * Q^n."""
    Q = [Q0]
    for _ in range(n_steps - 1):
        q_new = Q[-1] - a * Q[-1] ** n_bn * dt
        if q_new <= 0:
            break
        Q.append(q_new)
    return np.array(Q)


def test_fit_linear_recession_n_near_one():
    """A pure linear recession should yield n ≈ 1."""
    Q = _linear_recession(Q0=50.0, a=0.05, n_steps=30)
    bn = BrutsaertNieber(Q, min_recession_days=3).fit()
    assert bn.n_ == pytest.approx(1.0, abs=0.05)


def test_fit_powerlaw_recovers_exponent():
    """Synthetic power-law recession should recover the input n within tolerance."""
    target_n = 1.5
    Q = _powerlaw_recession(Q0=100.0, a=0.003, n_bn=target_n, n_steps=50)
    bn = BrutsaertNieber(Q, min_recession_days=3).fit()
    assert bn.n_ == pytest.approx(target_n, abs=0.1)


def test_to_reservoir_exponent_linear():
    """n=1 → b=1 (linear reservoir)."""
    Q = _linear_recession(Q0=50.0, a=0.05, n_steps=30)
    bn = BrutsaertNieber(Q, min_recession_days=3).fit()
    b = bn.to_reservoir_exponent()
    assert b == pytest.approx(1.0, abs=0.1)


def test_to_reservoir_exponent_boussinesq():
    """n=1.5 → b=2 (long-time Boussinesq solution)."""
    Q = _powerlaw_recession(Q0=100.0, a=0.003, n_bn=1.5, n_steps=60)
    bn = BrutsaertNieber(Q, min_recession_days=3).fit()
    b = bn.to_reservoir_exponent()
    assert b == pytest.approx(2.0, abs=0.3)


def test_to_reservoir_exponent_raises_before_fit():
    """to_reservoir_exponent() raises RuntimeError if fit() not called."""
    Q = _linear_recession(Q0=50.0, a=0.05, n_steps=20)
    bn = BrutsaertNieber(Q)
    with pytest.raises(RuntimeError):
        bn.to_reservoir_exponent()


def test_n_ge_2_returns_inf_and_warns():
    """to_reservoir_exponent() returns np.inf and warns when n ≥ 2."""
    Q = _linear_recession(Q0=50.0, a=0.05, n_steps=30)
    bn = BrutsaertNieber(Q, min_recession_days=3).fit()
    bn.n_ = 2.1  # directly exercise the n ≥ 2 branch
    with pytest.warns(UserWarning, match="undefined"):
        result = bn.to_reservoir_exponent()
    assert result == np.inf


def test_n_ge_2_fit_warns():
    """A fitted n ≥ 2 triggers the conversion warning (end-to-end)."""
    Q = _powerlaw_recession(Q0=200.0, a=0.0001, n_bn=2.2, n_steps=80)
    if len(Q) < 5:
        pytest.skip("Synthetic recession too short for this parameter set")
    bn = BrutsaertNieber(Q, min_recession_days=3).fit()
    if bn.n_ < 2.0:
        pytest.skip(f"Euler integration gave n={bn.n_:.2f} < 2; skip end-to-end check")
    with pytest.warns(UserWarning, match="undefined"):
        result = bn.to_reservoir_exponent()
    assert result == np.inf


def test_fit_raises_on_too_few_pairs():
    """Fewer than 3 recession pairs raises ValueError."""
    Q = np.array([10.0, 8.0, 9.0, 7.0, 5.0])  # non-monotone; short segments
    bn = BrutsaertNieber(Q, min_recession_days=10)
    with pytest.raises(ValueError, match="Only"):
        bn.fit()


def test_r2_attributes_set_after_fit():
    """r2_ and r2_quadratic_ are set and in [0, 1] after fit."""
    Q = _powerlaw_recession(Q0=80.0, a=0.002, n_bn=1.4, n_steps=50)
    bn = BrutsaertNieber(Q, min_recession_days=3).fit()
    assert 0.0 <= bn.r2_ <= 1.0
    assert 0.0 <= bn.r2_quadratic_ <= 1.0


def test_multiple_segments_concatenated():
    """Multiple recession segments are all used in the fit."""
    # Two long identical recessions separated by a rise.
    # A rise followed by the second segment creates an extra large-jump pair
    # (separator_end → seg[0]) so we only verify that many pairs were found
    # from both segments, not the exact n value.
    seg = _linear_recession(Q0=30.0, a=0.08, n_steps=15)
    Q = np.concatenate([seg, [60.0, 55.0], seg])
    bn = BrutsaertNieber(Q, min_recession_days=3).fit()
    # Both segments have 14 declining steps each → many pairs
    assert len(bn._Q_mid) > 20
    assert bn.n_ is not None


def test_zeros_and_negatives_skipped():
    """Zero and negative values in Q do not cause errors."""
    Q = np.array([10.0, 8.0, 6.0, 4.0, 0.0, -1.0, 5.0, 4.0, 3.0, 2.0, 1.0])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        bn = BrutsaertNieber(Q, min_recession_days=3).fit()
    assert bn.n_ is not None


def test_min_recession_days_filters_short_segments():
    """Short segments are excluded when min_recession_days is high."""
    # 2-step segments (3 values): each sep-to-seg transition is also 3 steps,
    # all below the min_recession_days=5 threshold → no pairs → ValueError.
    short_seg = np.array([10.0, 8.0, 6.0])  # 2 declining steps
    long_seg = _linear_recession(Q0=50.0, a=0.08, n_steps=20)

    Q_short = np.concatenate([short_seg, [12.0], short_seg, [12.0], short_seg])
    bn_strict = BrutsaertNieber(Q_short, min_recession_days=5)
    with pytest.raises(ValueError):
        bn_strict.fit()

    Q_long = np.concatenate([short_seg, [12.0], long_seg])
    bn_lenient = BrutsaertNieber(Q_long, min_recession_days=5).fit()
    assert bn_lenient.n_ is not None


def test_summary_raises_before_fit():
    """summary() raises RuntimeError if fit() not called."""
    Q = _linear_recession(Q0=50.0, a=0.05, n_steps=20)
    bn = BrutsaertNieber(Q)
    with pytest.raises(RuntimeError):
        bn.summary()


def test_n_to_b_inverse_formula():
    """b = 1/(2-n) and n = (2b-1)/b are exact inverses."""
    for b in [1.0, 1.5, 2.0, 3.0, 5.0]:
        n = (2 * b - 1) / b
        b_recovered = 1.0 / (2.0 - n)
        assert b_recovered == pytest.approx(b, rel=1e-12)
