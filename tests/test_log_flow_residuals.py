"""
Tests for mnished.calibration.log_flow_residual_terms.

The load-bearing property is that the helper reproduces run_and_score's
scoring mask *exactly* — otherwise a Dakota Bayesian likelihood built from its
residuals would be fit to a different set of days than the optimisation score.
"""

import os

import numpy as np
import pytest
import yaml

import mnished
from mnished.calibration import log_flow_residual_terms, _METRICS

EXAMPLE_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "examples", "cannon_forward")
)
CANNON_CSV = os.path.join(EXAMPLE_DIR, "CannonTestInput.csv")


def _cfg():
    """A small two-reservoir Cannon config (mirrors the back-compat tests)."""
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


def _write(tmp_path, cfg):
    p = tmp_path / "cfg.yml"
    p.write_text(yaml.safe_dump(cfg))
    return str(p)


def test_terms_reproduce_the_score(tmp_path):
    """KGE recomputed from the helper's (mod, obs) equals run_and_score's."""
    res = mnished.run_and_score(_write(tmp_path, _cfg()), metric='KGE')
    terms = log_flow_residual_terms(res)
    kge = _METRICS['KGE'](terms['mod'].to_numpy(), terms['obs'].to_numpy())
    assert kge == pytest.approx(res.score, abs=1e-12)


def test_terms_are_well_formed(tmp_path):
    res = mnished.run_and_score(_write(tmp_path, _cfg()), metric='KGE')
    terms = log_flow_residual_terms(res)
    assert len(terms) > 0
    assert np.all(np.isfinite(terms['residual']))
    # residual is exactly log_mod - log_obs
    np.testing.assert_array_equal(
        terms['residual'].to_numpy(),
        (terms['log_mod'] - terms['log_obs']).to_numpy())
    # the observed column is independent of the modelled flows
    assert set(terms.columns) >= {'date', 'obs', 'mod',
                                  'log_obs', 'log_mod', 'residual'}


def test_window_restricts_scored_days(tmp_path):
    """A start bound drops earlier days from the residual set."""
    res = mnished.run_and_score(_write(tmp_path, _cfg()), metric='KGE')
    full = log_flow_residual_terms(res)
    mid = full['date'].iloc[len(full) // 2]
    sub = log_flow_residual_terms(res, start=mid)
    assert 0 < len(sub) < len(full)
