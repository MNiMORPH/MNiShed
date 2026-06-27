"""
Generate the lake store-and-release figure for the documentation.

Runs a two-land-zone + lake configuration on the Cannon example and plots the
lake storage above its sill together with the lake-outlet discharge: storage
rises as the lake fills, and the threshold power-law outlet releases it slowly —
the store-and-release buffering that lets the land cascade keep a physical
recession timescale.

Usage (from the repository root)::

    python docs/figures/plot_lake.py

Writes ``docs/source/_static/lake_store_release.png``.
"""

import os
import warnings

import matplotlib
import pandas as pd
import yaml

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from mnished import Buckets  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
CANNON_CSV = os.path.join(ROOT, "examples", "cannon_forward", "CannonTestInput.csv")
OUT = os.path.join(ROOT, "docs", "source", "_static", "lake_store_release.png")
SILL = 180.0


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
            {"name": "upland", "area_fraction": 0.6,
             "reservoirs": {"recession_coefficients": [14, 500],
                            "exfiltration_fractions": [0.3, 1.0],
                            "maximum_effective_depths__mm": [20.0, float("inf")]},
             "initial_conditions": {"water_reservoir_effective_depths__mm": [8, 350]}},
            {"name": "lake", "kind": "lake", "area_fraction": 0.4,
             "lake": {"outflow_coefficient": 0.04, "sill_storage__mm": SILL,
                      "gw_partner": "upland", "f_route_lake": 0.5},
             "initial_conditions": {"lake_storage__mm": SILL + 30.0}},
        ],
    }


def main():
    path = os.path.join(HERE, "_lake_cfg.yml")
    with open(path, "w") as f:
        yaml.safe_dump(_cfg(), f)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Buckets()
        b.initialize(path, enforce_water_balance="none")
        b.run(store_depths=True, store_fluxes=True)
    os.remove(path)

    hd = b.hydrodata
    date = pd.to_datetime(hd["Date"])
    lake_col = [c for c in hd.columns if "lake" in c and "reservoir" in c][0]
    H = pd.to_numeric(hd[lake_col], errors="coerce")
    Qlake = pd.to_numeric(hd["Discharge: lake [mm/day]"], errors="coerce")

    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    ax.fill_between(date, H - SILL, 0, where=(H >= SILL), color="#bde0bd",
                    label="storage above sill")
    ax.axhline(0, color="#3a7", ls=":", lw=1)
    ax.set_ylabel("lake storage above sill (mm)", color="#2a7")
    ax.tick_params(axis="y", labelcolor="#2a7")
    ax2 = ax.twinx()
    ax2.plot(date, Qlake, color="#1f6fb2", lw=1.1, label="lake outlet Q")
    ax2.set_ylabel("lake-outlet discharge (mm/day)", color="#1f6fb2")
    ax2.tick_params(axis="y", labelcolor="#1f6fb2")
    ax.set_title("Lake store-and-release: stage above sill drives the outlet",
                 fontsize=10)
    fig.savefig(OUT, dpi=130, bbox_inches="tight")
    print(f"wrote {os.path.relpath(OUT, ROOT)}")


if __name__ == "__main__":
    main()
