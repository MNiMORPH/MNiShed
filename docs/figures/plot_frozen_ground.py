"""
Generate the frozen-ground figure for the documentation.

Frozen ground reduces infiltration (a high frozen-ground index routes more
meltwater straight to the stream), sharpening the spring freshet. This runs the
Cannon example with the frozen-ground module off and on and overlays the modeled
discharge so the freshet enhancement is visible.

Usage (from the repository root)::

    python docs/figures/plot_frozen_ground.py

Writes ``docs/source/_static/frozen_ground.png``.
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
OUT = os.path.join(ROOT, "docs", "source", "_static", "frozen_ground.png")


def _discharge(frozen):
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
        "snowmelt": {"PDD_melt_factor": 2.0, "fdd_threshold": 15.0},
        "modules": {"snowpack": True, "frozen_ground": frozen,
                    "rain_on_snow": True, "direct_runoff": False},
    }
    path = os.path.join(HERE, "_fg_cfg.yml")
    with open(path, "w") as f:
        yaml.safe_dump(cfg, f)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Buckets()
        b.initialize(path)
        b.run()
    os.remove(path)
    return (pd.to_datetime(b.hydrodata["Date"]),
            pd.to_numeric(b.hydrodata["Specific Discharge (modeled) [mm/day]"],
                          errors="coerce"))


def main():
    date, q_off = _discharge(False)
    _, q_on = _discharge(True)
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    ax.plot(date, q_off, color="#888", lw=1.1, label="frozen ground off")
    ax.plot(date, q_on, color="#1f6fb2", lw=1.3, label="frozen ground on")
    ax.set_ylabel("modeled discharge (mm/day)")
    ax.set_title("Frozen ground sharpens the spring freshet (Cannon)",
                 fontsize=10)
    ax.legend(fontsize=8)
    fig.savefig(OUT, dpi=130, bbox_inches="tight")
    print(f"wrote {os.path.relpath(OUT, ROOT)}")


if __name__ == "__main__":
    main()
