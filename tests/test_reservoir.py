"""
Unit tests for the Reservoir class.
"""

import numpy as np
import pytest
from hydroravens import Reservoir


def test_recharge_adds_water():
    res = Reservoir(t_efold=10.0)
    res.recharge(20.0)
    assert res.Hwater == pytest.approx(20.0)


def test_recharge_negative_creates_deficit():
    res = Reservoir(t_efold=10.0, H0=5.0)
    res.recharge(-10.0)
    assert res.Hwater == pytest.approx(0.0)
    assert res.H_deficit == pytest.approx(-5.0)


def test_recharge_overflow_hmax():
    res = Reservoir(t_efold=10.0, Hmax=50.0, H0=40.0)
    res.recharge(20.0)
    assert res.Hwater == pytest.approx(50.0)
    assert res.H_excess == pytest.approx(10.0)


def test_discharge_linear_exponential_decay():
    """Linear reservoir drains by exactly 1 - exp(-dt/tau) each step."""
    tau = 10.0
    H0 = 100.0
    res = Reservoir(t_efold=tau, H0=H0)
    res.discharge(dt=1.0)
    expected_dH = H0 * (1 - np.exp(-1.0 / tau))
    assert res.H_exfiltrated == pytest.approx(expected_dH, rel=1e-9)
    assert res.Hwater == pytest.approx(H0 - expected_dH, rel=1e-9)


def test_discharge_water_balance():
    """H_before == H_after + H_exfiltrated + H_infiltrated."""
    res = Reservoir(t_efold=10.0, f_to_discharge=0.6, H0=80.0)
    H_before = res.Hwater
    res.discharge(dt=1.0)
    assert H_before == pytest.approx(
        res.Hwater + res.H_exfiltrated + res.H_infiltrated, rel=1e-9
    )


def test_discharge_partitioning():
    """f_to_discharge splits drainage correctly between stream and infiltration."""
    f = 0.3
    res = Reservoir(t_efold=10.0, f_to_discharge=f, H0=100.0)
    res.discharge(dt=1.0)
    total = res.H_exfiltrated + res.H_infiltrated
    assert res.H_exfiltrated == pytest.approx(f * total, rel=1e-9)
    assert res.H_infiltrated == pytest.approx((1 - f) * total, rel=1e-9)


def test_discharge_empty_reservoir_produces_nothing():
    res = Reservoir(t_efold=10.0, H0=0.0)
    res.discharge(dt=1.0)
    assert res.H_exfiltrated == pytest.approx(0.0)
    assert res.Hwater == pytest.approx(0.0)


def test_mrt_linear_equals_t_efold():
    """MRT of a linear reservoir equals t_efold regardless of Q_ref."""
    res = Reservoir(t_efold=42.0)
    assert res.mean_residence_time(Q_ref=1.0) == pytest.approx(42.0)
    assert res.mean_residence_time(Q_ref=5.0) == pytest.approx(42.0)


def test_mrt_nonlinear_formula():
    """MRT = tau^(1/b) / Q_ref^(1 - 1/b) for b > 1."""
    tau = 10.0
    b = 2.0
    Q_ref = 2.0
    res = Reservoir(t_efold=tau)
    res.recession_exponent = b
    expected = tau ** (1.0 / b) / Q_ref ** (1.0 - 1.0 / b)
    assert res.mean_residence_time(Q_ref=Q_ref) == pytest.approx(expected)


def test_mrt_raises_on_nonpositive_q():
    res = Reservoir(t_efold=10.0)
    with pytest.raises(ValueError):
        res.mean_residence_time(Q_ref=0.0)


def test_invalid_t_efold_raises():
    with pytest.raises(ValueError):
        Reservoir(t_efold=0.0)


def test_invalid_f_to_discharge_raises():
    with pytest.raises(ValueError):
        Reservoir(t_efold=10.0, f_to_discharge=1.5)


def test_tile_drain_increases_discharge():
    """With tile drainage, total discharge exceeds non-tile case."""
    res_plain = Reservoir(t_efold=10.0, f_to_discharge=0.5, H0=100.0)
    res_plain.discharge(dt=1.0)

    res_tile = Reservoir(t_efold=10.0, f_to_discharge=0.5,
                         f_tile=0.5, tau_tile=5.0, H0=100.0)
    res_tile.discharge(dt=1.0)

    assert res_tile.H_discharge > res_plain.H_discharge
