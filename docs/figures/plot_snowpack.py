"""
Generate the snowpack figure for the documentation.

Runs the Cannon example with the snowpack on and plots modeled snow-water
equivalent (SWE) against mean temperature, showing degree-day accumulation when
cold and melt when the temperature rises above freezing.

Usage (from the repository root)::

    python docs/figures/plot_snowpack.py

Writes ``docs/source/_static/snowpack_swe.png``.
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
OUT = os.path.join(ROOT, "docs", "source", "_static", "snowpack_swe.png")


def main():
    cfg = {
        "timeseries": {"datafile": CANNON_CSV},
        "catchment": {"drainage_basin_area__km2": 3800,
                      "evapotranspiration_method": "datafile",
                      "water_year_start_month": 10},
        "general": {"spin_up_cycles": 0},
        "reservoirs": {"recession_coefficients": [14, 500],
                       "exfiltration_fractions": [0.3, 1.0],
                       "maximum_effective_depths__mm": [float("inf"), float("inf")]},
        "initial_conditions": {"water_reservoir_effective_depths__mm": [15, 400],
                               "snowpack__mm_SWE": 0},
        "snowmelt": {"PDD_melt_factor": 2.0},
        "modules": {"snowpack": True, "frozen_ground": False,
                    "rain_on_snow": True, "direct_runoff": False},
    }
    path = os.path.join(HERE, "_snow_cfg.yml")
    with open(path, "w") as f:
        yaml.safe_dump(cfg, f)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Buckets()
        b.initialize(path)
        b.run()
    os.remove(path)

    hd = b.hydrodata
    date = pd.to_datetime(hd["Date"])
    swe = pd.to_numeric(hd["Snowpack (modeled) [mm SWE]"], errors="coerce")
    temp = pd.to_numeric(hd["Mean Temperature [C]"], errors="coerce")

    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    ax.fill_between(date, swe, color="#9ecae1", label="modeled SWE")
    ax.set_ylabel("SWE (mm)", color="#3182bd")
    ax.tick_params(axis="y", labelcolor="#3182bd")
    ax2 = ax.twinx()
    ax2.plot(date, temp, color="#d95f0e", lw=0.7, alpha=0.8, label="mean T")
    ax2.axhline(0, color="gray", ls=":", lw=0.8)
    ax2.set_ylabel("mean temperature (°C)", color="#d95f0e")
    ax2.tick_params(axis="y", labelcolor="#d95f0e")
    ax.set_title("Snowpack: degree-day accumulation and melt (Cannon)",
                 fontsize=10)
    fig.savefig(OUT, dpi=130, bbox_inches="tight")
    print(f"wrote {os.path.relpath(OUT, ROOT)}  (peak SWE {swe.max():.0f} mm)")


if __name__ == "__main__":
    main()
