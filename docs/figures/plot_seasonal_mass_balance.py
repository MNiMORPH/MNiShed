"""
Generate the seasonal mass-balance diagnostic figure for the documentation.

Runs a small two-land-zone + lake configuration on the Cannon example forcing
with ``store_fluxes=True`` and plots the per-season discharge split by source
(fast/event, slow/baseflow, lake outlet) as stacked bars against observed Q —
the "which source carries each season" read that ``SeasonalMassBalance`` gives.

Usage (from the repository root)::

    python docs/figures/plot_seasonal_mass_balance.py

Writes ``docs/source/_static/seasonal_mass_balance.png``.
"""

import os
import warnings

import matplotlib
import yaml

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from mnished import SeasonalMassBalance, run_and_score  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
CANNON_CSV = os.path.join(ROOT, "examples", "cannon_forward", "CannonTestInput.csv")
OUT = os.path.join(ROOT, "docs", "source", "_static", "seasonal_mass_balance.png")


def _cfg():
    return {
        "timeseries": {"datafile": CANNON_CSV},
        "catchment": {"drainage_basin_area__km2": 3800,
                      "evapotranspiration_method": "datafile",
                      "water_year_start_month": 10},
        "general": {"spin_up_cycles": 0, "enforce_water_balance": "none"},
        "snowmelt": {"PDD_melt_factor": 1.0},
        "modules": {"snowpack": False, "frozen_ground": False,
                    "rain_on_snow": False, "direct_runoff": False},
        "sub_catchments": [
            {"name": "upland", "area_fraction": 0.3,
             "reservoirs": {"recession_coefficients": [14, 500],
                            "exfiltration_fractions": [0.3, 1.0],
                            "maximum_effective_depths__mm": [20.0, float("inf")]},
             "initial_conditions": {"water_reservoir_effective_depths__mm": [8, 350]}},
            {"name": "lake_basin", "area_fraction": 0.4,
             "reservoirs": {"recession_coefficients": [14, 500],
                            "exfiltration_fractions": [0.3, 1.0],
                            "maximum_effective_depths__mm": [20.0, float("inf")]},
             "initial_conditions": {"water_reservoir_effective_depths__mm": [8, 350]}},
            {"name": "lake", "kind": "lake", "area_fraction": 0.3,
             "lake": {"outflow_coefficient": 0.05, "sill_storage__mm": 180.0,
                      "gw_partner": "lake_basin", "f_route_lake": 0.5},
             "initial_conditions": {"lake_storage__mm": 260.0}},
        ],
    }


def main():
    path = os.path.join(HERE, "_smb_cfg.yml")
    with open(path, "w") as f:
        yaml.safe_dump(_cfg(), f)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = run_and_score(path, enforce_water_balance="none", metric="KGE",
                            store_fluxes=True)
    os.remove(path)
    st = SeasonalMassBalance(res.buckets).seasonal_table()

    seasons = list(st.index)
    x = range(len(seasons))
    fig, ax = plt.subplots(figsize=(5.6, 3.4))
    fast, slow, lake = st["fast"].values, st["slow"].values, st["lake"].values
    ax.bar(x, fast, color="#5aa0d0", label="fast (event)")
    ax.bar(x, slow, bottom=fast, color="#2e6b9e", label="slow (baseflow)")
    ax.bar(x, lake, bottom=fast + slow, color="#7fc97f", label="lake outlet")
    ax.plot(x, st["obs"].values, "o-", color="k", lw=1.5, ms=5, label="observed Q")
    ax.set_xticks(list(x))
    ax.set_xticklabels(seasons)
    ax.set_ylabel("discharge (mm/day)")
    ax.set_title("Seasonal mass balance — modeled Q by source vs. observed",
                 fontsize=10)
    ax.legend(fontsize=7.5, loc="upper right", framealpha=0.9)
    fig.savefig(OUT, dpi=130, bbox_inches="tight")
    print(f"wrote {os.path.relpath(OUT, ROOT)}")


if __name__ == "__main__":
    main()
