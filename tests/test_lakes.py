"""
Tests for lake (open-water) sub-catchments.

A ``kind: lake`` sub-catchment is a single threshold power-law reservoir fed by
direct ``P - E`` (open-water evaporation = the model ET column, no soil-moisture
water stress) and discharging over a sill, ``Q_out = a*max(H - H_sill, 0)^b``
(b = 5/3 by default). v1 holds the lake hydrologically disconnected from river
inflow (``f_route_lake = 0``); the bidirectional groundwater exchange ``Q_gw``
to a land sub-catchment is added in a later stage. See ``DESIGN_lakes.md``.

Forcing uses the Cannon River forward example (examples/cannon_forward/).
"""

import os

import numpy as np
import pandas as pd
import pytest
import yaml

import mnished

EXAMPLE_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "examples", "cannon_forward")
)
CANNON_CSV = os.path.join(EXAMPLE_DIR, "CannonTestInput.csv")

Q_COL = "Specific Discharge (modeled) [mm/day]"


def _base_cfg():
    """Cannon-forced skeleton with snowpack off and no water-balance closure,
    so the mass-balance test is a real check of the bookkeeping (not an
    ET auto-tune that trivially forces P - Q - ET = 0)."""
    return {
        'timeseries': {'datafile': CANNON_CSV},
        'catchment': {'drainage_basin_area__km2': 3800,
                      'evapotranspiration_method': 'datafile',
                      'water_year_start_month': 10},
        'general': {'spin_up_cycles': 0, 'enforce_water_balance': 'none'},
        'snowmelt': {'PDD_melt_factor': 1.0},
        'modules': {'snowpack': False, 'frozen_ground': False,
                    'rain_on_snow': False, 'direct_runoff': False},
    }


def _land(area, recession=(14, 500), exfil=(0.3, 1.0),
          hmax=(20.0, float('inf')), h0=(5, 300), name='upland'):
    return {'name': name, 'area_fraction': area,
            'reservoirs': {'recession_coefficients': list(recession),
                           'exfiltration_fractions': list(exfil),
                           'maximum_effective_depths__mm': list(hmax)},
            'initial_conditions': {
                'water_reservoir_effective_depths__mm': list(h0)}}


def _lake(area, a=0.05, sill=200.0, b=5.0 / 3.0, h0=250.0, name='lake',
          **extra):
    lk = {'outflow_coefficient': a, 'sill_storage__mm': sill,
          'outflow_exponent': b}
    lk.update(extra)
    return {'name': name, 'kind': 'lake', 'area_fraction': area,
            'lake': lk, 'initial_conditions': {'lake_storage__mm': h0}}


def _write(tmp_path, cfg, name='cfg.yml'):
    p = tmp_path / name
    p.write_text(yaml.safe_dump(cfg))
    return str(p)


def _num(hd, col):
    return pd.to_numeric(hd[col], errors='coerce').to_numpy()


def _init(tmp_path, sub_catchments):
    cfg = _base_cfg()
    cfg['sub_catchments'] = sub_catchments
    b = mnished.Buckets()
    b.initialize(_write(tmp_path, cfg))
    return b


# ---------------------------------------------------------------------------
# Structure
# ---------------------------------------------------------------------------

def test_lake_config_builds(tmp_path):
    b = _init(tmp_path, [_land(0.7), _lake(0.3, a=0.05, sill=200.0, b=1.6)])
    assert b.has_lake
    assert [sc.kind for sc in b.sub_catchments] == ['land', 'lake']
    lake = b.sub_catchments[1]
    assert len(lake.reservoirs) == 1
    res = lake.reservoirs[0]
    assert res.junction_type == 'threshold'
    assert res.H_threshold == 200.0
    assert res.recession_exponent == 1.6
    assert res.recession_coeff == pytest.approx(1.0 / 0.05)   # a = 1/recession_coeff
    assert lake.snowpack is None
    assert lake.f_route_lake == 0.0
    # Flattened reservoirs include the lake's single store.
    assert len(b.reservoirs) == 3


def test_lake_default_exponent_is_five_thirds(tmp_path):
    b = _init(tmp_path, [_land(0.6), {'name': 'lk', 'kind': 'lake',
              'area_fraction': 0.4,
              'lake': {'outflow_coefficient': 0.05, 'sill_storage__mm': 100.0},
              'initial_conditions': {'lake_storage__mm': 150.0}}])
    assert b.sub_catchments[1].reservoirs[0].recession_exponent == \
        pytest.approx(5.0 / 3.0)


def test_lake_forces_python_path(tmp_path):
    """A lake present means has_lake is True; run() must use the pure-Python
    loop (the JIT does not yet cover lakes — issue #19)."""
    b = _init(tmp_path, [_land(0.7), _lake(0.3)])
    assert b.has_lake is True
    b.run()   # must not raise regardless of numba availability
    assert np.isfinite(_num(b.hydrodata, Q_COL)).any()


# ---------------------------------------------------------------------------
# Water balance and outflow physics
# ---------------------------------------------------------------------------

def test_lake_mass_balance_closes(tmp_path):
    b = _init(tmp_path, [_land(0.7), _lake(0.3, sill=200.0, h0=250.0)])

    def storage():
        return sum(sc.area_fraction * sum(r.Hwater for r in sc.reservoirs)
                   for sc in b.sub_catchments)
    s0 = storage()
    b.run()
    s1 = storage()
    hd = b.hydrodata
    P, ET, Qm = (_num(hd, 'Precipitation [mm/day]'),
                 _num(hd, 'ET for model [mm/day]'), _num(hd, Q_COL))
    mask = np.isfinite(Qm)
    residual = np.nansum((P - ET)[mask]) - np.nansum(Qm[mask]) - (s1 - s0)
    assert residual == pytest.approx(0.0, abs=1e-6)


def test_lake_dead_pool_holds_below_sill(tmp_path):
    """With the sill unreachable, the lake never discharges above it. Two land
    zones leave the lake partner ambiguous (no Q_gw), isolating the outlet."""
    b = _init(tmp_path, [
        _land(0.6, recession=(14,), exfil=(1.0,), hmax=(float('inf'),),
              h0=(10,), name='a'),
        _land(0.3, recession=(14,), exfil=(1.0,), hmax=(float('inf'),),
              h0=(10,), name='b'),
        _lake(0.1, sill=1e6, h0=100.0)])
    assert b.sub_catchments[2].gw_partner is None
    b.run()
    lake_res = b.sub_catchments[2].reservoirs[0]
    # Never rose above the (unreachable) sill, so no over-sill discharge path.
    assert lake_res.Hwater < 1e6
    assert lake_res.H_discharge == pytest.approx(0.0)


def test_lake_discharges_toward_sill(tmp_path):
    """A lake started above its sill drains toward it."""
    b = _init(tmp_path, [_land(0.7), _lake(0.3, sill=200.0, h0=400.0)])
    b.run()
    assert b.sub_catchments[1].reservoirs[0].Hwater < 400.0


# ---------------------------------------------------------------------------
# Validation / guards
# ---------------------------------------------------------------------------

def test_f_route_lake_nonzero_not_implemented(tmp_path):
    with pytest.raises(NotImplementedError):
        _init(tmp_path, [_land(0.7), _lake(0.3, f_route_lake=0.2)])


def test_lake_requires_positive_outflow_coefficient(tmp_path):
    with pytest.raises(ValueError):
        _init(tmp_path, [_land(0.7), {'name': 'lk', 'kind': 'lake',
              'area_fraction': 0.3,
              'lake': {'outflow_coefficient': 0.0, 'sill_storage__mm': 100.0}}])


def test_lake_area_counts_toward_basin_sum(tmp_path):
    """Lake area is part of the basin; land + lake must sum to 1."""
    with pytest.raises(ValueError):
        _init(tmp_path, [_land(0.7), _lake(0.5)])   # sums to 1.2


# ---------------------------------------------------------------------------
# Groundwater exchange Q_gw (bidirectional, lake <-> land subsurface)
# ---------------------------------------------------------------------------

def _isolated_exchange(tmp_path, h_s, H_lake, a_land=0.7, a_lake=0.3):
    """A land + lake basin with the heads set and the sill unreachable (no
    outlet), so a single _apply_lake_gw_exchange() isolates Q_gw."""
    b = _init(tmp_path, [
        _land(a_land, recession=(500,), exfil=(1.0,), hmax=(float('inf'),),
              h0=(h_s,), name='L'),
        _lake(a_lake, sill=1e9, h0=H_lake, name='lk')])
    return b


def test_gw_partner_auto_resolves_single_land(tmp_path):
    b = _init(tmp_path, [_land(0.7), _lake(0.3)])
    lake = b.sub_catchments[1]
    assert lake.gw_partner is b.sub_catchments[0]


def test_gw_partner_named(tmp_path):
    b = _init(tmp_path, [
        _land(0.4, name='till'), _land(0.3, name='clay'),
        _lake(0.3, gw_partner='clay')])
    assert b.sub_catchments[2].gw_partner.name == 'clay'


def test_gw_partner_missing_raises(tmp_path):
    with pytest.raises(ValueError):
        _init(tmp_path, [_land(0.7), _lake(0.3, gw_partner='nope')])


def test_gw_partner_cannot_be_a_lake(tmp_path):
    with pytest.raises(ValueError):
        _init(tmp_path, [_land(0.4), _lake(0.3, name='lk1'),
                         _lake(0.3, name='lk2', gw_partner='lk1')])


def test_gw_no_partner_when_ambiguous(tmp_path):
    """Several land zones and no name given -> no exchange (partner None)."""
    b = _init(tmp_path, [_land(0.4, name='a'), _land(0.3, name='b'),
                         _lake(0.3)])
    assert b.sub_catchments[2].gw_partner is None


def test_gw_flows_into_lake_when_aquifer_higher(tmp_path):
    b = _isolated_exchange(tmp_path, h_s=300.0, H_lake=100.0)
    land, lake = b.sub_catchments[0].reservoirs[-1], b.sub_catchments[1].reservoirs[0]
    b._apply_lake_gw_exchange()
    assert land.Hwater < 300.0   # aquifer loses
    assert lake.Hwater > 100.0   # lake gains


def test_gw_flows_out_of_lake_when_lake_higher(tmp_path):
    b = _isolated_exchange(tmp_path, h_s=100.0, H_lake=300.0)
    land, lake = b.sub_catchments[0].reservoirs[-1], b.sub_catchments[1].reservoirs[0]
    b._apply_lake_gw_exchange()
    assert land.Hwater > 100.0   # aquifer gains
    assert lake.Hwater < 300.0   # lake loses (sign flip from one term)


@pytest.mark.parametrize("h_s,H_lake", [(300.0, 100.0), (100.0, 300.0)])
def test_gw_conserves_volume(tmp_path, h_s, H_lake):
    a_land, a_lake = 0.7, 0.3
    b = _isolated_exchange(tmp_path, h_s, H_lake, a_land, a_lake)
    land, lake = b.sub_catchments[0].reservoirs[-1], b.sub_catchments[1].reservoirs[0]
    v0 = a_land * land.Hwater + a_lake * lake.Hwater
    b._apply_lake_gw_exchange()
    v1 = a_land * land.Hwater + a_lake * lake.Hwater
    assert v1 == pytest.approx(v0, abs=1e-9)


def test_gw_mass_balance_closes(tmp_path):
    """Q_gw is internal; the basin mass balance still closes with it active."""
    b = _init(tmp_path, [_land(0.7), _lake(0.3, sill=200.0, h0=250.0)])
    assert b.sub_catchments[1].gw_partner is not None

    def storage():
        return sum(sc.area_fraction * sum(r.Hwater for r in sc.reservoirs)
                   for sc in b.sub_catchments)
    s0 = storage()
    b.run()
    s1 = storage()
    hd = b.hydrodata
    P, ET, Qm = (_num(hd, 'Precipitation [mm/day]'),
                 _num(hd, 'ET for model [mm/day]'), _num(hd, Q_COL))
    mask = np.isfinite(Qm)
    residual = np.nansum((P - ET)[mask]) - np.nansum(Qm[mask]) - (s1 - s0)
    assert residual == pytest.approx(0.0, abs=1e-6)


def test_gw_changes_lake_trajectory(tmp_path):
    """Turning the exchange off (no partner) gives a different lake state than
    leaving it on, with a large head contrast to drive exchange."""
    common = dict(sill=1e9, h0=50.0)   # sill unreachable -> only P-E and Q_gw act
    on = _init(tmp_path, [
        _land(0.5, recession=(500,), exfil=(1.0,), hmax=(float('inf'),),
              h0=(800,), name='a'),
        _lake(0.5, name='lk', **common)])
    off = _init(tmp_path, [
        _land(0.4, recession=(500,), exfil=(1.0,), hmax=(float('inf'),),
              h0=(800,), name='a'),
        _land(0.1, recession=(500,), exfil=(1.0,), hmax=(float('inf'),),
              h0=(0,), name='b'),     # second land zone -> lake partner ambiguous -> off
        _lake(0.5, name='lk', **common)])
    assert on.sub_catchments[1].gw_partner is not None
    assert off.sub_catchments[2].gw_partner is None
    on.run()
    off.run()
    h_on = on.sub_catchments[1].reservoirs[0].Hwater
    h_off = off.sub_catchments[2].reservoirs[0].Hwater
    assert h_on != pytest.approx(h_off)


# ---------------------------------------------------------------------------
# Calibration, mass balance, snowpack interaction
# ---------------------------------------------------------------------------

def _calib_cfg(tmp_path):
    """Land + lake config with snowpack ON and observed discharge (Cannon),
    written to disk for run_and_score / initialize()."""
    cfg = {
        'timeseries': {'datafile': CANNON_CSV},
        'catchment': {'drainage_basin_area__km2': 3800,
                      'evapotranspiration_method': 'datafile',
                      'water_year_start_month': 10},
        'general': {'spin_up_cycles': 0},
        'snowmelt': {'PDD_melt_factor': 1.0, 'fgi_decay_coeff': 0.97,
                     'snow_insulation_k': 0.0},
        'modules': {'snowpack': True, 'frozen_ground': True,
                    'rain_on_snow': True, 'direct_runoff': False},
        'sub_catchments': [
            {'name': 'upland', 'area_fraction': 0.7,
             'reservoirs': {'recession_coefficients': [14, 500],
                            'exfiltration_fractions': [0.3, 1.0],
                            'maximum_effective_depths__mm':
                                [float('inf'), float('inf')]},
             'initial_conditions': {
                 'water_reservoir_effective_depths__mm': [15, 400]}},
            {'name': 'lake', 'kind': 'lake', 'area_fraction': 0.3,
             'lake': {'outflow_coefficient': 0.05, 'sill_storage__mm': 200.0},
             'initial_conditions': {'lake_storage__mm': 250.0}}],
    }
    return _write(tmp_path, cfg)


def test_lake_initializes_with_snowpack(tmp_path):
    """A lake config with snowpack enabled must initialize (the lake carries no
    snowpack; regression for the per-sub-catchment SWE-init loop)."""
    b = mnished.Buckets()
    b.initialize(_calib_cfg(tmp_path))
    assert b.has_snowpack
    assert b.sub_catchments[1].snowpack is None
    b.run()   # must not raise
    assert b.check_mass_balance() == pytest.approx(0.0, abs=1e-6)


def test_lake_calibrates_via_sub_catchment_overrides(tmp_path):
    """The lake outflow coefficient (a -> recession_coeff = 1/a) and sill
    (H_threshold) calibrate through the existing per-sub-catchment override
    path; the exponent stays fixed at 5/3."""
    from mnished.calibration import _apply_reservoir_overrides
    cfg_path = _calib_cfg(tmp_path)
    res = mnished.run_and_score(
        cfg_path, spin_up_cycles=1, metric='KGE',
        sub_catchments=[{'recession_coeff': [12.0, 480.0]},
                        {'recession_coeff': [1.0 / 0.08],
                         'H_threshold': [180.0]}])
    assert np.isfinite(res.score)

    # The override reaches the lake reservoir and leaves b = 5/3 untouched.
    b = mnished.Buckets()
    b.initialize(cfg_path)
    lake_res = b.sub_catchments[1].reservoirs[0]
    _apply_reservoir_overrides(b.sub_catchments[1].reservoirs,
                               {'recession_coeff': [1.0 / 0.08],
                                'H_threshold': [180.0]})
    assert lake_res.recession_coeff == pytest.approx(1.0 / 0.08)
    assert lake_res.H_threshold == 180.0
    assert lake_res.recession_exponent == pytest.approx(5.0 / 3.0)
