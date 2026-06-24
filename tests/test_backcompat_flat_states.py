"""
Backward-compatibility sentinels for the K=1 *flat* state shape — REMOVE OR
INVERT AT v4.0 (MNiMORPH/MNiShed#18).

Parallel sub-catchments shipped with a back-compat shim: for a single
sub-catchment, the public state representation stays flat/scalar
(``{'reservoirs': [...], 'snowpack': float, 'fgi': float}``) and only becomes
nested (``{'sub_catchments': [...]}``) for several. v4.0 will remove the flat
form and make the representation uniformly per-sub-catchment.

The tests below pin the *current* (v3.x) back-compat behaviour, so they act as
an executable contract for the shim:

* Tests marked ``SENTINEL`` assert the deprecated flat behaviour. They pass
  while the shim is in place and **fail by design once #18 removes it** — at
  which point delete them (or invert to assert the nested-only behaviour).
* The one ``FORWARD`` test asserts the nested contract that should *survive*
  v4.0 unchanged; it should keep passing across the transition.

If a sentinel here fails unexpectedly (i.e. not as part of doing #18), the
back-compat shim regressed.
"""

import os
import warnings

import pytest
import yaml

import mnished

EXAMPLE_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "examples", "cannon_forward")
)
CANNON_CSV = os.path.join(EXAMPLE_DIR, "CannonTestInput.csv")

DEPRECATION_MATCH = "MNiMORPH/MNiShed#18"


def _legacy_cfg():
    """A single-cascade (K=1) config: two reservoirs, Cannon forcing."""
    return {
        'timeseries': {'datafile': CANNON_CSV},
        'catchment': {
            'drainage_basin_area__km2': 3800,
            'evapotranspiration_method': 'datafile',
            'water_year_start_month': 10,
        },
        'general': {'spin_up_cycles': 0},
        'reservoirs': {
            'recession_coefficients': [14, 500],
            'exfiltration_fractions': [0.3, 1.0],
            'maximum_effective_depths__mm': [float('inf'), float('inf')],
        },
        'initial_conditions': {
            'water_reservoir_effective_depths__mm': [15, 400],
            'snowpack__mm_SWE': 0,
        },
        'snowmelt': {'PDD_melt_factor': 1.0, 'fgi_decay_coeff': 0.97,
                     'snow_insulation_k': 0.0},
        'modules': {'snowpack': True, 'frozen_ground': True,
                    'rain_on_snow': True, 'direct_runoff': False},
    }


def _write(tmp_path, cfg, name="cfg.yml"):
    p = tmp_path / name
    p.write_text(yaml.safe_dump(cfg))
    return str(p)


def _flat_states():
    return {'reservoirs': [20.0, 300.0], 'snowpack': 0.0, 'fgi': 0.0,
            'H_deficit_carry': 0.0}


def _nested_states():
    return {'sub_catchments': [
        {'reservoirs': [20.0, 300.0], 'snowpack': 0.0, 'fgi': 0.0,
         'H_deficit_carry': 0.0}]}


def _no_deprecation(record):
    return [w for w in record
            if issubclass(w.category, DeprecationWarning)
            and DEPRECATION_MATCH in str(w.message)]


# ---------------------------------------------------------------------------
# SENTINEL — flat behaviour; FAILS BY DESIGN once #18 removes the shim
# ---------------------------------------------------------------------------

def test_final_states_is_flat_for_single_sub_catchment(tmp_path):
    """SENTINEL: run_and_score returns the flat state shape at K=1.

    v4.0 (#18): final_states becomes nested under 'sub_catchments' always.
    """
    result = mnished.run_and_score(_write(tmp_path, _legacy_cfg()),
                                   spin_up_cycles=1, metric='KGE')
    fs = result.final_states
    assert 'sub_catchments' not in fs
    assert set(fs) == {'reservoirs', 'snowpack', 'fgi', 'H_deficit_carry'}
    assert len(fs['reservoirs']) == 2


def test_flat_initial_states_accepted_and_warns(tmp_path):
    """SENTINEL: a flat initial_states is accepted and emits the v4.0
    DeprecationWarning.

    v4.0 (#18): the flat form is removed (rejected, no warning).
    """
    with pytest.warns(DeprecationWarning, match=DEPRECATION_MATCH):
        result = mnished.run_and_score(
            _write(tmp_path, _legacy_cfg()),
            initial_states=_flat_states(), spin_up_cycles=0, metric='KGE')
    assert isinstance(result.final_states, dict)


def test_flat_post_spinup_states_accepted_and_warns(tmp_path):
    """SENTINEL: a flat post_spinup_states is accepted and warns (decade mode).

    v4.0 (#18): removed.
    """
    with pytest.warns(DeprecationWarning, match=DEPRECATION_MATCH):
        mnished.run_and_score(
            _write(tmp_path, _legacy_cfg()),
            start='1993-01-01', end='1994-12-31',
            post_spinup_states={'reservoirs': [None, 350.0]},
            spin_up_cycles=1, metric='KGE')


def test_buckets_scalar_alias_properties(tmp_path):
    """SENTINEL: Buckets exposes scalar snowpack/_fgi/H_deficit_carry aliases
    over the (single) sub-catchment, with broadcast setters.

    v4.0 (#18): these basin-level aliases may be removed or replaced with
    explicit per-sub-catchment / area-weighted accessors.
    """
    b = mnished.Buckets()
    b.initialize(_write(tmp_path, _legacy_cfg()))
    # Read: scalar aliases of the first (only) sub-catchment.
    assert isinstance(b._fgi, float)
    assert isinstance(b.H_deficit_carry, float)
    assert b.snowpack is b.sub_catchments[0].snowpack
    # Write: broadcasts to every sub-catchment; round-trips at K=1.
    b._fgi = 3.5
    assert b.sub_catchments[0].fgi == 3.5
    b.H_deficit_carry = 2.0
    assert b.H_deficit_carry == 2.0
    assert b.sub_catchments[0].H_deficit_carry == 2.0


def test_store_depths_legacy_columns_at_single_sub_catchment(tmp_path):
    """SENTINEL: store_depths uses the unlabeled 'H_reservoir_{i}' columns at
    K=1.

    v4.0 (#18): columns may be sub-catchment-labeled uniformly.
    """
    b = mnished.Buckets()
    b.initialize(_write(tmp_path, _legacy_cfg()))
    b.run(store_depths=True)
    assert 'H_reservoir_0 (modeled) [mm]' in b.hydrodata.columns
    assert 'H_reservoir_1 (modeled) [mm]' in b.hydrodata.columns


# ---------------------------------------------------------------------------
# FORWARD — nested contract; should SURVIVE v4.0 unchanged
# ---------------------------------------------------------------------------

def test_nested_initial_states_accepted_without_warning(tmp_path):
    """FORWARD: the nested state form works at K=1 and does NOT warn.

    This is the v4.0-and-beyond contract; it should keep passing across the
    transition. (If this one starts failing, the forward path regressed.)
    """
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        result = mnished.run_and_score(
            _write(tmp_path, _legacy_cfg()),
            initial_states=_nested_states(), spin_up_cycles=0, metric='KGE')
    assert not _no_deprecation(rec), "nested initial_states must not deprecate"
    assert isinstance(result.final_states, dict)
