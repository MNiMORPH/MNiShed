"""
Regression tests for the Buckets end-to-end model.

Uses the Cannon River forward example (examples/cannon_forward/).

Forcing data are derived from Livneh et al. (2015) gridded meteorology,
extracted for the Cannon River catchment by Pierce et al. (2021) via
https://github.com/MNiMORPH/LivnehPierce-hydro-extractor.

References
----------
Livneh, B., et al. (2015). A spatially comprehensive, hydrometeorological
data set for Mexico, the U.S., and Southern Canada 1950–2013. Scientific
Data, 2, 150042. https://doi.org/10.1038/sdata.2015.42

Pierce, A., et al. (2021). LivnehPierce hydro extractor.
https://github.com/MNiMORPH/LivnehPierce-hydro-extractor
"""

import os

import numpy as np
import pytest

EXAMPLE_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "examples", "cannon_forward")
)
EXAMPLE_CONFIG = "cannon_cfg.yml"


@pytest.fixture(scope="module")
def cannon():
    """Run the Cannon River forward example once; share across tests."""
    import hydroravens

    orig_dir = os.getcwd()
    os.chdir(EXAMPLE_DIR)
    try:
        b = hydroravens.Buckets()
        b.initialize(EXAMPLE_CONFIG)
        b.run()
    finally:
        os.chdir(orig_dir)
    return b


# ---------------------------------------------------------------------------
# Smoke / output columns
# ---------------------------------------------------------------------------

def test_run_completes(cannon):
    """Model runs from initialize through run() without error."""
    assert cannon.hydrodata is not None
    assert len(cannon.hydrodata) > 0


def test_output_columns_present(cannon):
    """Expected output columns are present in hydrodata."""
    df = cannon.hydrodata
    assert "Specific Discharge (modeled) [mm/day]" in df.columns
    assert "Snowpack (modeled) [mm SWE]" in df.columns
    assert "Subsurface storage (modeled total) [mm]" in df.columns
    assert "ET for model [mm/day]" in df.columns


def test_discharge_non_negative(cannon):
    """Modeled discharge is non-negative for every timestep."""
    q = cannon.hydrodata["Specific Discharge (modeled) [mm/day]"]
    assert (q >= 0).all()


def test_swe_non_negative(cannon):
    """SWE is non-negative for every timestep."""
    swe = cannon.hydrodata["Snowpack (modeled) [mm SWE]"]
    assert (swe >= 0).all()


def test_storage_non_negative(cannon):
    """Total subsurface storage is non-negative for every timestep."""
    stor = cannon.hydrodata["Subsurface storage (modeled total) [mm]"]
    assert (stor >= 0).all()


# ---------------------------------------------------------------------------
# NSE regression
# ---------------------------------------------------------------------------

def test_nse_regression(cannon):
    """NSE is within 0.01 of the known uncalibrated value (-3.6592)."""
    nse = cannon.compute_NSE(return_nse=True, verbose=False)
    assert nse == pytest.approx(-3.659183, abs=0.01)


# ---------------------------------------------------------------------------
# Water balance
# ---------------------------------------------------------------------------

def test_water_balance_closure(cannon):
    """P - Q - ET - ΔS - ΔSWe < 0.1% of total P."""
    df = cannon.hydrodata
    total_P = df["Precipitation [mm/day]"].sum()
    total_Q = df["Specific Discharge (modeled) [mm/day]"].sum()
    total_ET = df["ET for model [mm/day]"].sum()
    delta_S = (
        df["Subsurface storage (modeled total) [mm]"].iloc[-1]
        - df["Subsurface storage (modeled total) [mm]"].iloc[0]
    )
    delta_swe = (
        df["Snowpack (modeled) [mm SWE]"].iloc[-1]
        - df["Snowpack (modeled) [mm SWE]"].iloc[0]
    )
    wb_error_fraction = abs(total_P - total_Q - total_ET - delta_S - delta_swe) / total_P
    assert wb_error_fraction < 0.001  # < 0.1% of P


# ---------------------------------------------------------------------------
# Basic physical checks
# ---------------------------------------------------------------------------

def test_discharge_finite(cannon):
    """All modeled discharge values are finite."""
    q = cannon.hydrodata["Specific Discharge (modeled) [mm/day]"].astype(float)
    assert np.isfinite(q).all()


def test_et_positive(cannon):
    """ET is positive over the run period."""
    et = cannon.hydrodata["ET for model [mm/day]"]
    assert et.sum() > 0


def test_snowpack_peaks_in_winter(cannon):
    """SWE peaks in winter (Dec–Feb), not in summer (Jun–Aug)."""
    import pandas as pd
    dates = pd.to_datetime(cannon.hydrodata["Date"])
    swe = cannon.hydrodata["Snowpack (modeled) [mm SWE]"].values
    winter_mask = dates.dt.month.isin([12, 1, 2])
    summer_mask = dates.dt.month.isin([6, 7, 8])
    assert swe[winter_mask].mean() > swe[summer_mask].mean()


def test_three_reservoirs_initialized(cannon):
    """Cannon River config uses 3 reservoirs."""
    assert len(cannon.reservoirs) == 3


def test_record_length(cannon):
    """Cannon River example spans 3 years (1992–1994, ~1096 days including leap year)."""
    n = len(cannon.hydrodata)
    assert 1050 <= n <= 1150
