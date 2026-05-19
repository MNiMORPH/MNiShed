"""
Unit tests for the Snowpack class.
"""

import pytest
from hydroravens.hydroravens import Snowpack


def test_accumulation_below_freezing():
    """Precipitation accumulates as SWE when T <= 0."""
    sp = Snowpack(melt_factor=2.0)
    sp.set_temperature(-5.0)
    sp.recharge(10.0)
    assert sp.Hwater == pytest.approx(10.0)
    assert sp.H_infiltrated == pytest.approx(0.0)


def test_precipitation_bypasses_at_positive_temperature():
    """Rainfall passes through directly when T > 0."""
    sp = Snowpack(melt_factor=2.0)
    sp.set_temperature(5.0)
    sp.recharge(10.0)
    assert sp.Hwater == pytest.approx(0.0)
    assert sp.H_infiltrated == pytest.approx(10.0)


def test_melt_reduces_swe():
    """Positive temperature melts SWE at the PDD rate."""
    sp = Snowpack(melt_factor=2.0)
    sp.Hwater = 20.0
    sp.set_temperature(3.0)
    sp.recharge(0.0)
    sp.melt(dt=1.0)
    expected_melt = 2.0 * 3.0 * 1.0  # melt_factor * T * dt = 6 mm
    assert sp.Hwater == pytest.approx(20.0 - expected_melt)
    assert sp.H_infiltrated == pytest.approx(expected_melt)


def test_melt_exhausts_snowpack():
    """Melt cannot reduce SWE below zero; excess energy is returned."""
    sp = Snowpack(melt_factor=2.0)
    sp.Hwater = 3.0
    sp.set_temperature(5.0)
    sp.recharge(0.0)
    excess_dd = sp.melt(dt=1.0)
    assert sp.Hwater == pytest.approx(0.0)
    assert sp.H_infiltrated == pytest.approx(3.0)
    # 2.0 * 5.0 * 1.0 = 10 mm available, only 3 mm SWE: excess = 7 mm / 2.0 mf
    assert excess_dd == pytest.approx(7.0 / 2.0)


def test_no_melt_below_freezing():
    """No melt when T <= 0."""
    sp = Snowpack(melt_factor=2.0)
    sp.Hwater = 50.0
    sp.set_temperature(-3.0)
    sp.recharge(0.0)
    excess_dd = sp.melt(dt=1.0)
    assert sp.Hwater == pytest.approx(50.0)
    assert excess_dd == pytest.approx(0.0)


def test_sublimation_removes_swe():
    """Negative recharge removes water from snowpack as sublimation."""
    sp = Snowpack(melt_factor=2.0)
    sp.Hwater = 10.0
    sp.set_temperature(-5.0)
    sp.recharge(-3.0)
    assert sp.Hwater == pytest.approx(7.0)
    assert sp.H_deficit == pytest.approx(0.0)


def test_sublimation_deficit_when_insufficient_swe():
    """Negative recharge exceeding SWE produces a deficit."""
    sp = Snowpack(melt_factor=2.0)
    sp.Hwater = 2.0
    sp.set_temperature(-5.0)
    sp.recharge(-5.0)
    assert sp.Hwater == pytest.approx(0.0)
    assert sp.H_deficit == pytest.approx(-3.0)
