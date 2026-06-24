"""
Tests for parallel sub-catchments.

A SubCatchment is a parallel hydraulic compartment of the basin with its own
reservoir cascade; basin discharge is the area-weighted mean of the
sub-catchments. A single sub-catchment of area_fraction 1.0 reproduces the
original (non-partitioned) behaviour exactly.

Forcing uses the Cannon River forward example (examples/cannon_forward/).
"""

import os

import numpy as np
import pytest
import yaml

EXAMPLE_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "examples", "cannon_forward")
)
CANNON_CSV = os.path.join(EXAMPLE_DIR, "CannonTestInput.csv")

Q_COL = "Specific Discharge (modeled) [mm/day]"
S_COL = "Subsurface storage (modeled total) [mm]"


def _base_cfg(spin_up_cycles=0):
    """Cannon-forced config skeleton without a reservoir structure."""
    return {
        'timeseries': {'datafile': CANNON_CSV},
        'catchment': {
            'drainage_basin_area__km2': 3800,
            'evapotranspiration_method': 'datafile',
            'water_year_start_month': 10,
        },
        'general': {'spin_up_cycles': spin_up_cycles},
        'snowmelt': {'PDD_melt_factor': 1.0, 'fgi_decay_coeff': 0.97,
                     'snow_insulation_k': 0.0},
        'modules': {'snowpack': True, 'frozen_ground': True,
                    'rain_on_snow': True, 'direct_runoff': False},
    }


def _legacy_cfg(recession, exfil, hmax, h0):
    cfg = _base_cfg()
    cfg['reservoirs'] = {
        'recession_coefficients': recession,
        'exfiltration_fractions': exfil,
        'maximum_effective_depths__mm': hmax,
    }
    cfg['initial_conditions'] = {
        'water_reservoir_effective_depths__mm': h0,
        'snowpack__mm_SWE': 0,
    }
    return cfg


def _write(tmp_path, cfg, name="cfg.yml"):
    p = tmp_path / name
    p.write_text(yaml.safe_dump(cfg))
    return str(p)


def _run(cfg_path, cycles=3):
    import mnished
    b = mnished.Buckets()
    b.initialize(cfg_path)
    for _ in range(cycles):
        b.run()
    b.run()
    return b


# ---------------------------------------------------------------------------
# 1. Legacy single-cascade configs are one sub-catchment of area 1.0
# ---------------------------------------------------------------------------

def test_legacy_yaml_is_single_sub_catchment(tmp_path):
    cfg = _legacy_cfg([14, 40, 500], [0.19, 0.76, 1.0],
                      [18.0, float('inf'), float('inf')], [15, 40, 500])
    import mnished
    b = mnished.Buckets()
    b.initialize(_write(tmp_path, cfg))
    assert b.n_sub_catchments == 1
    assert b.sub_catchments[0].area_fraction == 1.0
    assert len(b.reservoirs) == 3
    assert b.sub_catchments[0].reservoirs is not None
    assert len(b.sub_catchments[0].reservoirs) == 3


# ---------------------------------------------------------------------------
# 2. A two-sub-catchment config loads with the right structure
# ---------------------------------------------------------------------------

def _two_sc_cfg(area_a=0.55, area_b=0.45):
    cfg = _base_cfg()
    cfg['sub_catchments'] = [
        {'name': 'till', 'area_fraction': area_a,
         'reservoirs': {'recession_coefficients': [14, 500],
                        'exfiltration_fractions': [0.3, 1.0],
                        'maximum_effective_depths__mm': [20.0, float('inf')]},
         'initial_conditions': {'water_reservoir_effective_depths__mm': [5, 300]}},
        {'name': 'clay', 'area_fraction': area_b,
         'reservoirs': {'recession_coefficients': [1500],
                        'exfiltration_fractions': [1.0],
                        'maximum_effective_depths__mm': [float('inf')]}},
    ]
    return cfg


def test_two_sub_catchments_load(tmp_path):
    import mnished
    b = mnished.Buckets()
    b.initialize(_write(tmp_path, _two_sc_cfg()))
    assert b.n_sub_catchments == 2
    assert [sc.name for sc in b.sub_catchments] == ['till', 'clay']
    assert [sc.area_fraction for sc in b.sub_catchments] == [0.55, 0.45]
    assert [len(sc.reservoirs) for sc in b.sub_catchments] == [2, 1]
    assert len(b.reservoirs) == 3            # flattened across sub-catchments
    # Per-sub-catchment initial depths; the second defaults to zero.
    assert [r.Hwater for r in b.sub_catchments[0].reservoirs] == [5, 300]
    assert [r.Hwater for r in b.sub_catchments[1].reservoirs] == [0.0]


# ---------------------------------------------------------------------------
# 3. Validation
# ---------------------------------------------------------------------------

def test_area_fractions_must_sum_to_one(tmp_path):
    import mnished
    cfg = _two_sc_cfg(area_a=0.55, area_b=0.55)
    with pytest.raises(ValueError, match="sum to 1"):
        mnished.Buckets().initialize(_write(tmp_path, cfg))


def test_duplicate_names_raise(tmp_path):
    import mnished
    cfg = _two_sc_cfg()
    cfg['sub_catchments'][1]['name'] = 'till'
    with pytest.raises(ValueError, match="unique"):
        mnished.Buckets().initialize(_write(tmp_path, cfg))


def test_empty_reservoirs_raise(tmp_path):
    import mnished
    cfg = _two_sc_cfg()
    cfg['sub_catchments'][1]['reservoirs'] = {
        'recession_coefficients': [], 'exfiltration_fractions': [],
        'maximum_effective_depths__mm': []}
    with pytest.raises(ValueError, match="at least one reservoir"):
        mnished.Buckets().initialize(_write(tmp_path, cfg))


# ---------------------------------------------------------------------------
# 4. Mass balance is conserved across sub-catchments
# ---------------------------------------------------------------------------

def test_mass_balance_conserved(tmp_path):
    """A two-sub-catchment run conserves mass (excess ≈ 0 vs total flux)."""
    b = _run(_write(tmp_path, _two_sc_cfg()))
    excess = b.check_mass_balance()
    total_P = b.hydrodata['Precipitation [mm/day]'].sum()
    assert abs(excess) < 1e-6 * total_P


# ---------------------------------------------------------------------------
# 5. Two identical sub-catchments reproduce one cascade
# ---------------------------------------------------------------------------

def test_two_identical_sub_catchments_equal_one(tmp_path):
    legacy = _legacy_cfg([14, 500], [0.3, 1.0], [20.0, float('inf')], [5, 300])
    b_one = _run(_write(tmp_path, legacy, "one.yml"))
    q_one = b_one.hydrodata[Q_COL].astype(float).to_numpy()

    cfg = _base_cfg()

    def _sc():
        return {'reservoirs': {'recession_coefficients': [14, 500],
                               'exfiltration_fractions': [0.3, 1.0],
                               'maximum_effective_depths__mm': [20.0, float('inf')]},
                'initial_conditions': {
                    'water_reservoir_effective_depths__mm': [5, 300]}}
    cfg['sub_catchments'] = [dict(name='a', area_fraction=0.6, **_sc()),
                             dict(name='b', area_fraction=0.4, **_sc())]
    b_two = _run(_write(tmp_path, cfg, "two.yml"))
    q_two = b_two.hydrodata[Q_COL].astype(float).to_numpy()

    np.testing.assert_allclose(q_two, q_one, rtol=1e-9, atol=1e-9, equal_nan=True)


# ---------------------------------------------------------------------------
# 6. Parallel sub-catchments differ from a serial cascade
# ---------------------------------------------------------------------------

def test_parallel_differs_from_serial(tmp_path):
    """Two reservoirs run as parallel sub-catchments are not the same model as
    the two reservoirs in a vertical cascade."""
    serial = _legacy_cfg([14, 500], [0.3, 1.0], [float('inf'), float('inf')],
                         [50, 50])
    b_serial = _run(_write(tmp_path, serial, "serial.yml"))
    q_serial = b_serial.hydrodata[Q_COL].astype(float).to_numpy()

    cfg = _base_cfg()
    cfg['sub_catchments'] = [
        {'name': 'fast', 'area_fraction': 0.5,
         'reservoirs': {'recession_coefficients': [14],
                        'exfiltration_fractions': [1.0],
                        'maximum_effective_depths__mm': [float('inf')]},
         'initial_conditions': {'water_reservoir_effective_depths__mm': [50]}},
        {'name': 'slow', 'area_fraction': 0.5,
         'reservoirs': {'recession_coefficients': [500],
                        'exfiltration_fractions': [1.0],
                        'maximum_effective_depths__mm': [float('inf')]},
         'initial_conditions': {'water_reservoir_effective_depths__mm': [50]}},
    ]
    b_parallel = _run(_write(tmp_path, cfg, "parallel.yml"))
    q_parallel = b_parallel.hydrodata[Q_COL].astype(float).to_numpy()

    assert np.nanmax(np.abs(q_serial - q_parallel)) > 1e-3


# ---------------------------------------------------------------------------
# 7. Multipath drainage works inside a sub-catchment
# ---------------------------------------------------------------------------

def test_multipath_in_sub_catchment(tmp_path):
    """A sub-catchment with a multipath drain responds differently from one
    without — the per-sub-catchment cascade honours multipath parameters."""
    def cfg(multipath):
        c = _base_cfg()
        sc1 = {'name': 'a', 'area_fraction': 0.5,
               'reservoirs': {'recession_coefficients': [200],
                              'exfiltration_fractions': [1.0],
                              'maximum_effective_depths__mm': [float('inf')]},
               'initial_conditions': {'water_reservoir_effective_depths__mm': [100]}}
        sc2 = {'name': 'b', 'area_fraction': 0.5,
               'reservoirs': {'recession_coefficients': [200],
                              'exfiltration_fractions': [1.0],
                              'maximum_effective_depths__mm': [float('inf')]},
               'initial_conditions': {'water_reservoir_effective_depths__mm': [100]}}
        if multipath:
            sc2['reservoirs']['multipath_thresholds__mm'] = [20.0]
            sc2['reservoirs']['multipath_timescales__days'] = [3.0]
        c['sub_catchments'] = [sc1, sc2]
        return c

    q_plain = _run(_write(tmp_path, cfg(False), "plain.yml")
                   ).hydrodata[Q_COL].astype(float).to_numpy()
    q_mp = _run(_write(tmp_path, cfg(True), "mp.yml")
                ).hydrodata[Q_COL].astype(float).to_numpy()
    assert np.nanmax(np.abs(q_plain - q_mp)) > 1e-3


# ---------------------------------------------------------------------------
# 8. Per-sub-catchment forcing is reserved but not yet implemented
# ---------------------------------------------------------------------------

def test_per_sub_catchment_forcing_not_implemented(tmp_path):
    import mnished
    cfg = _two_sc_cfg()
    cfg['sub_catchments'][0]['forcing'] = {'datafile': 'till_forcing.csv'}
    with pytest.raises(NotImplementedError):
        mnished.Buckets().initialize(_write(tmp_path, cfg))


# ---------------------------------------------------------------------------
# 9. JIT and pure-Python agree for multiple sub-catchments
# ---------------------------------------------------------------------------

def test_jit_matches_pure_python_sub_catchments(tmp_path, monkeypatch):
    """The Numba and pure-Python loops agree for two sub-catchments with
    differing cascades and the advanced reservoir mechanics enabled."""
    pytest.importorskip("numba", exc_type=ImportError)
    import mnished
    import mnished.mnished as _m

    cfg = _base_cfg()
    cfg['general']['et_alpha'] = 0.6
    cfg['modules']['et_reservoir_draw'] = True
    cfg['sub_catchments'] = [
        {'name': 'till', 'area_fraction': 0.55,
         'reservoirs': {
             'recession_coefficients': [14, 500],
             'exfiltration_fractions': [0.3, 1.0],
             'maximum_effective_depths__mm': [20.0, float('inf')],
             'recession_exponents': [2.0, 1.0],
             'junction_types': ['threshold', 'fraction'],
             'H_threshold__mm': [3.0, 0.0],
             'tile_fractions': [0.3, 0.0],
             'tile_residence_times__days': [3.0, None],
             'multipath_thresholds__mm': [10.0, None],
             'multipath_timescales__days': [5.0, None]},
         'initial_conditions': {
             'water_reservoir_effective_depths__mm': [5, 300],
             'snowpack__mm_SWE': 2.0}},
        {'name': 'clay', 'area_fraction': 0.45,
         'reservoirs': {
             'recession_coefficients': [40, 1500],
             'exfiltration_fractions': [0.5, 1.0],
             'maximum_effective_depths__mm': [float('inf'), float('inf')],
             'junction_types': ['leakance', 'fraction'],
             'leakance_R__days': [50.0, None]},
         'initial_conditions': {
             'water_reservoir_effective_depths__mm': [30, 400],
             'snowpack__mm_SWE': 1.0}},
    ]
    cfg_path = _write(tmp_path, cfg, "advanced_sc.yml")

    b_jit = mnished.Buckets()
    b_jit.initialize(cfg_path)
    b_jit.run()
    q_jit = b_jit.hydrodata[Q_COL].astype(float).to_numpy()
    s_jit = b_jit.hydrodata[S_COL].astype(float).to_numpy()

    monkeypatch.setattr(_m, "_numba_available", False)
    b_py = mnished.Buckets()
    b_py.initialize(cfg_path)
    b_py.run()
    q_py = b_py.hydrodata[Q_COL].astype(float).to_numpy()
    s_py = b_py.hydrodata[S_COL].astype(float).to_numpy()

    np.testing.assert_allclose(q_jit, q_py, rtol=1e-7, atol=1e-9, equal_nan=True)
    np.testing.assert_allclose(s_jit, s_py, rtol=1e-7, atol=1e-9, equal_nan=True)


# ---------------------------------------------------------------------------
# 10. run_and_score drives a multi-sub-catchment calibration
# ---------------------------------------------------------------------------

def test_run_and_score_sub_catchments(tmp_path):
    import mnished
    cfg = _two_sc_cfg()
    cfg['general']['spin_up_cycles'] = 0
    cfg_path = _write(tmp_path, cfg)

    res = mnished.run_and_score(
        cfg_path, spin_up_cycles=1, metric='KGE',
        sub_catchments=[{'area_fraction': 0.6, 'recession_coeff': [20, 400]},
                        {'area_fraction': 0.4, 'recession_coeff': [1200]}])
    assert np.isfinite(res.score)
    b = res.buckets
    assert [sc.area_fraction for sc in b.sub_catchments] == [0.6, 0.4]
    assert [r.recession_coeff for r in b.sub_catchments[0].reservoirs] == [20, 400]
    assert [r.recession_coeff for r in b.sub_catchments[1].reservoirs] == [1200]


def test_run_and_score_rejects_flat_and_structured(tmp_path):
    import mnished
    cfg = _two_sc_cfg()
    cfg['general']['spin_up_cycles'] = 0
    cfg_path = _write(tmp_path, cfg)
    with pytest.raises(ValueError, match="inside each sub-catchment"):
        mnished.run_and_score(
            cfg_path, spin_up_cycles=0, recession_coeff=[1, 2, 3],
            sub_catchments=[{'recession_coeff': [20, 400]},
                            {'recession_coeff': [1200]}])
