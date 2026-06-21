"""
Unit tests for the Reservoir class.
"""

import numpy as np
import pytest
from mnished import Reservoir


def test_recharge_adds_water():
    res = Reservoir(recession_coeff=10.0)
    res.recharge(20.0)
    assert res.Hwater == pytest.approx(20.0)


def test_recharge_negative_creates_deficit():
    res = Reservoir(recession_coeff=10.0, H0=5.0)
    res.recharge(-10.0)
    assert res.Hwater == pytest.approx(0.0)
    assert res.H_deficit == pytest.approx(-5.0)


def test_recharge_overflow_hmax():
    res = Reservoir(recession_coeff=10.0, Hmax=50.0, H0=40.0)
    res.recharge(20.0)
    assert res.Hwater == pytest.approx(50.0)
    assert res.H_excess == pytest.approx(10.0)


def test_discharge_linear_exponential_decay():
    """Linear reservoir drains by exactly 1 - exp(-dt/tau) each step."""
    tau = 10.0
    H0 = 100.0
    res = Reservoir(recession_coeff=tau, H0=H0)
    res.discharge(dt=1.0)
    expected_dH = H0 * (1 - np.exp(-1.0 / tau))
    assert res.H_exfiltrated == pytest.approx(expected_dH, rel=1e-9)
    assert res.Hwater == pytest.approx(H0 - expected_dH, rel=1e-9)


def test_discharge_water_balance():
    """H_before == H_after + H_exfiltrated + H_to_next."""
    res = Reservoir(recession_coeff=10.0, f_to_discharge=0.6, H0=80.0)
    H_before = res.Hwater
    res.discharge(dt=1.0)
    assert H_before == pytest.approx(
        res.Hwater + res.H_exfiltrated + res.H_to_next, rel=1e-9
    )


def test_discharge_partitioning():
    """f_to_discharge splits drainage correctly between stream and infiltration."""
    f = 0.3
    res = Reservoir(recession_coeff=10.0, f_to_discharge=f, H0=100.0)
    res.discharge(dt=1.0)
    total = res.H_exfiltrated + res.H_to_next
    assert res.H_exfiltrated == pytest.approx(f * total, rel=1e-9)
    assert res.H_to_next == pytest.approx((1 - f) * total, rel=1e-9)


def test_discharge_empty_reservoir_produces_nothing():
    res = Reservoir(recession_coeff=10.0, H0=0.0)
    res.discharge(dt=1.0)
    assert res.H_exfiltrated == pytest.approx(0.0)
    assert res.Hwater == pytest.approx(0.0)


def test_mrt_linear_equals_recession_coeff():
    """MRT of a linear reservoir equals recession_coeff regardless of Q_ref."""
    res = Reservoir(recession_coeff=42.0)
    assert res.mean_residence_time(Q_ref=1.0) == pytest.approx(42.0)
    assert res.mean_residence_time(Q_ref=5.0) == pytest.approx(42.0)


def test_mrt_nonlinear_formula():
    """MRT = tau^(1/b) / Q_ref^(1 - 1/b) for b > 1."""
    tau = 10.0
    b = 2.0
    Q_ref = 2.0
    res = Reservoir(recession_coeff=tau)
    res.recession_exponent = b
    expected = tau ** (1.0 / b) / Q_ref ** (1.0 - 1.0 / b)
    assert res.mean_residence_time(Q_ref=Q_ref) == pytest.approx(expected)


def test_mrt_raises_on_nonpositive_q():
    res = Reservoir(recession_coeff=10.0)
    with pytest.raises(ValueError):
        res.mean_residence_time(Q_ref=0.0)


def test_invalid_recession_coeff_raises():
    with pytest.raises(ValueError):
        Reservoir(recession_coeff=0.0)


def test_old_t_efold_kwarg_raises():
    """t_efold is no longer a valid keyword; passing it must raise TypeError."""
    with pytest.raises(TypeError):
        Reservoir(t_efold=10.0)


def test_invalid_f_to_discharge_raises():
    with pytest.raises(ValueError):
        Reservoir(recession_coeff=10.0, f_to_discharge=1.5)


def test_nonlinear_discharge_exceeds_linear():
    """b > 1 drains faster than b = 1 from the same initial storage."""
    H0 = 4.0
    tau = 10.0
    res_linear = Reservoir(recession_coeff=tau, H0=H0)
    res_linear.discharge(dt=1.0)

    res_nonlinear = Reservoir(recession_coeff=tau, H0=H0)
    res_nonlinear.recession_exponent = 2.0
    res_nonlinear.discharge(dt=1.0)

    assert res_nonlinear.H_exfiltrated > res_linear.H_exfiltrated


def test_nonlinear_discharge_exact_formula():
    """Nonlinear (b=2) discharge matches the exact analytical solution."""
    tau = 10.0
    b = 2.0
    H0 = 4.0
    dt = 1.0
    # H_ref = 1.0 (default); tau_eff = tau * H_ref^(b-1) = tau
    # H_new = [H0^(1-b) + (b-1)*dt/tau_eff]^(1/(1-b))
    tau_eff = tau  # * 1.0^(b-1)
    H_new = (H0 ** (1 - b) + (b - 1) * dt / tau_eff) ** (1 / (1 - b))
    expected_dH = H0 - H_new

    res = Reservoir(recession_coeff=tau, H0=H0)
    res.recession_exponent = b
    res.discharge(dt=dt)

    assert res.H_exfiltrated == pytest.approx(expected_dH, rel=1e-9)
    assert res.Hwater == pytest.approx(H0 - expected_dH, rel=1e-9)


def test_nonlinear_discharge_water_balance():
    """Water balance holds for b > 1: H_before == H_after + dH."""
    res = Reservoir(recession_coeff=10.0, f_to_discharge=0.6, H0=4.0)
    res.recession_exponent = 2.0
    H_before = res.Hwater
    res.discharge(dt=1.0)
    assert H_before == pytest.approx(
        res.Hwater + res.H_exfiltrated + res.H_to_next, rel=1e-9
    )


def test_tile_drain_increases_discharge():
    """With tile drainage, total discharge exceeds non-tile case."""
    res_plain = Reservoir(recession_coeff=10.0, f_to_discharge=0.5, H0=100.0)
    res_plain.discharge(dt=1.0)

    res_tile = Reservoir(recession_coeff=10.0, f_to_discharge=0.5,
                         f_tile=0.5, tau_tile=5.0, H0=100.0)
    res_tile.discharge(dt=1.0)

    assert res_tile.H_discharge > res_plain.H_discharge


# --- Junction type tests ---

def test_leakance_invalid_without_R():
    """leakance junction requires leakance_R."""
    with pytest.raises(ValueError):
        Reservoir(recession_coeff=10.0, junction_type='leakance')


def test_invalid_junction_type_raises():
    with pytest.raises(ValueError):
        Reservoir(recession_coeff=10.0, junction_type='magic')


def test_leakance_all_to_stream_when_reservoirs_equal():
    """When H_this == H_next, Q_leak = 0: all drainage goes to stream."""
    res = Reservoir(recession_coeff=10.0, leakance_R=100.0,
                    junction_type='leakance', H0=50.0)
    res.discharge(dt=1.0, H_next=50.0)
    assert res.H_to_next == pytest.approx(0.0, abs=1e-12)
    assert res.H_exfiltrated == pytest.approx(res.H_discharge, abs=1e-12)


def test_leakance_all_infiltrates_when_head_difference_large():
    """When head difference / R >> dH, essentially all drainage is leakance."""
    H0 = 50.0
    res = Reservoir(recession_coeff=1000.0, leakance_R=0.01,
                    junction_type='leakance', H0=H0)
    res.discharge(dt=1.0, H_next=0.0)
    # Q_leak = (50 - 0) / 0.01 = 5000 >> dH; should be capped at dH
    assert res.H_exfiltrated == pytest.approx(0.0, abs=1e-9)
    assert res.H_to_next == pytest.approx(
        res.H_to_next + res.H_exfiltrated, rel=1e-9)


def test_leakance_water_balance():
    """H_before == H_after + H_exfiltrated + H_to_next for leakance."""
    res = Reservoir(recession_coeff=10.0, leakance_R=50.0,
                    junction_type='leakance', H0=80.0)
    H_before = res.Hwater
    res.discharge(dt=1.0, H_next=20.0)
    assert H_before == pytest.approx(
        res.Hwater + res.H_exfiltrated + res.H_to_next, rel=1e-9)


def test_leakance_more_to_next_when_head_difference_larger():
    """Larger head difference → more leakance to next reservoir."""
    res_small = Reservoir(recession_coeff=10.0, leakance_R=100.0,
                          junction_type='leakance', H0=80.0)
    res_small.discharge(dt=1.0, H_next=70.0)

    res_large = Reservoir(recession_coeff=10.0, leakance_R=100.0,
                          junction_type='leakance', H0=80.0)
    res_large.discharge(dt=1.0, H_next=10.0)

    assert res_large.H_to_next > res_small.H_to_next


def test_threshold_no_drainage_below_threshold():
    """Below H_threshold, no water drains."""
    res = Reservoir(recession_coeff=10.0, junction_type='threshold',
                    H_threshold=100.0, H0=50.0)
    res.discharge(dt=1.0)
    assert res.H_exfiltrated == pytest.approx(0.0, abs=1e-12)
    assert res.Hwater == pytest.approx(50.0, abs=1e-12)


def test_threshold_normal_drainage_above_threshold():
    """Above H_threshold, only the excess drains (dead storage preserved)."""
    H_threshold = 20.0
    H0 = 70.0
    tau = 10.0
    res = Reservoir(recession_coeff=tau, junction_type='threshold',
                    H_threshold=H_threshold, H0=H0)
    res.discharge(dt=1.0)
    H_eff = H0 - H_threshold
    expected_dH = H_eff * (1 - np.exp(-1.0 / tau))
    assert res.H_exfiltrated == pytest.approx(expected_dH, rel=1e-9)
    assert res.Hwater == pytest.approx(H0 - expected_dH, rel=1e-9)


def test_threshold_water_balance():
    """Water balance holds for threshold junction."""
    res = Reservoir(recession_coeff=10.0, junction_type='threshold',
                    H_threshold=30.0, H0=80.0)
    H_before = res.Hwater
    res.discharge(dt=1.0)
    assert H_before == pytest.approx(
        res.Hwater + res.H_exfiltrated + res.H_to_next, rel=1e-9)


# --- Fast-path discharge components (recorded for the BMI flux partition) ---

def test_tile_component_recorded():
    """H_tile records the tile-drain contribution; H_multipath stays 0."""
    res = Reservoir(recession_coeff=10.0, f_to_discharge=0.5,
                    f_tile=0.4, tau_tile=3.0, H0=100.0)
    res.discharge(dt=1.0)
    assert res.H_tile > 0.0
    assert res.H_multipath == 0.0


def test_multipath_component_recorded():
    """H_multipath records the multipath contribution; H_tile stays 0."""
    res = Reservoir(recession_coeff=50.0, multipath_threshold=20.0,
                    multipath_timescale=3.0, H0=100.0)
    res.discharge(dt=1.0)
    assert res.H_multipath > 0.0
    assert res.H_tile == 0.0


def test_no_fast_path_components_zero():
    """With no tile or multipath, both component records are exactly 0."""
    res = Reservoir(recession_coeff=10.0, f_to_discharge=1.0, H0=100.0)
    res.discharge(dt=1.0)
    assert res.H_tile == 0.0
    assert res.H_multipath == 0.0
