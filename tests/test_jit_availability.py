"""
The pure-Python time-loop fallback is not silent when it matters.

Buckets.run() warns once per process when it falls back to the pure-Python loop
because Numba is installed but fails to import (most often a NumPy/Numba version
mismatch). A plain "Numba not installed" stays quiet — pure Python is the
expected default without the ``jit`` extra. (PDM and et_water_stress are
JIT-supported, so they no longer trigger a fallback notice.)

These tests drive the notice by monkeypatching the module-level numba state, so
they are deterministic regardless of whether numba is installed in the test
environment.
"""

import os
import warnings

import pytest
import yaml

import mnished
import mnished.mnished as M

EXAMPLE_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "examples", "cannon_forward")
)
CANNON_CSV = os.path.join(EXAMPLE_DIR, "CannonTestInput.csv")


def _legacy_cfg():
    return {
        'timeseries': {'datafile': CANNON_CSV},
        'catchment': {'drainage_basin_area__km2': 3800,
                      'evapotranspiration_method': 'datafile',
                      'water_year_start_month': 10},
        'general': {'spin_up_cycles': 0},
        'reservoirs': {'recession_coefficients': [14, 500],
                       'exfiltration_fractions': [0.3, 1.0],
                       'maximum_effective_depths__mm': [float('inf'),
                                                        float('inf')]},
        'initial_conditions': {
            'water_reservoir_effective_depths__mm': [15, 400],
            'snowpack__mm_SWE': 0},
        'snowmelt': {'PDD_melt_factor': 1.0, 'fgi_decay_coeff': 0.97,
                     'snow_insulation_k': 0.0},
        'modules': {'snowpack': True, 'frozen_ground': True,
                    'rain_on_snow': True, 'direct_runoff': False},
    }


def _buckets(tmp_path):
    p = tmp_path / "cfg.yml"
    p.write_text(yaml.safe_dump(_legacy_cfg()))
    b = mnished.Buckets()
    b.initialize(str(p))
    return b


def _pure_python_warnings(record):
    return [w for w in record if "pure-Python" in str(w.message)]


def _set_broken_numba(monkeypatch):
    """Simulate numba installed but failing to import (e.g. NumPy too new)."""
    monkeypatch.setattr(M, "_numba_available", False)
    monkeypatch.setattr(M, "_numba_import_error",
                        ImportError("Numba needs NumPy 2.2 or less"))
    monkeypatch.setattr(M, "_jit_unavailable_notified", False)


def test_installed_but_broken_numba_notifies(tmp_path, monkeypatch):
    _set_broken_numba(monkeypatch)
    b = _buckets(tmp_path)
    with pytest.warns(UserWarning, match="failed to import"):
        b.run()


def test_numba_not_installed_is_silent(tmp_path, monkeypatch):
    """numba simply absent -> pure Python is expected -> no notice."""
    monkeypatch.setattr(M, "_numba_available", False)
    monkeypatch.setattr(M, "_numba_import_error", None)
    monkeypatch.setattr(M, "_jit_unavailable_notified", False)
    b = _buckets(tmp_path)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        b.run()
    assert not _pure_python_warnings(rec)


def test_notice_fires_at_most_once_per_process(tmp_path, monkeypatch):
    _set_broken_numba(monkeypatch)
    b = _buckets(tmp_path)
    with pytest.warns(UserWarning, match="pure-Python"):
        b.run()                                   # first run: notifies
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        b.run()                                   # second run: silent
    assert not _pure_python_warnings(rec)
