"""
Generate the parameter-identifiability figure for the documentation.

Builds a small Cannon calibration over three parameters and shows the two
identifiability views around that point: the curvature eigenspectrum (stiffness
per eigen-direction, and the loadings that name each combination) and a 2-D
ridge over a degenerate pair (the soil recession timescale against the shallow
exfiltration fraction), where the fit barely changes along a diagonal valley.

Usage (from the repository root)::

    python docs/figures/plot_identifiability.py

Writes ``docs/source/_static/identifiability.png``.
"""

import os
import warnings

import matplotlib
import yaml

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from mnished import Calibrator, eigenspectrum, ridge  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
CANNON_CSV = os.path.join(ROOT, "examples", "cannon_forward", "CannonTestInput.csv")
OUT = os.path.join(ROOT, "docs", "source", "_static", "identifiability.png")


def _calibrator():
    model = {
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
        "snowmelt": {"PDD_melt_factor": 1.0},
        "modules": {"snowpack": True, "frozen_ground": False,
                    "rain_on_snow": True, "direct_runoff": False},
    }
    model_path = os.path.join(HERE, "_ident_model.yml")
    with open(model_path, "w") as f:
        yaml.safe_dump(model, f)
    params = {
        "log__tau_soil": {"lower": 0.3, "upper": 2.0, "initial": 1.15,
                          "target": "recession_coeff[0]"},
        "log__tau_gw":   {"lower": 2.0, "upper": 3.0, "initial": 2.6,
                          "target": "recession_coeff[1]"},
        "f_exf_shallow": {"lower": 0.05, "upper": 0.95, "initial": 0.45,
                          "target": "f_to_discharge[0]"},
    }
    driver = {"config_template": model_path, "metric": "KGE",
              "spin_up_cycles": 0, "routing_N": 2}
    cal = Calibrator(params, driver, model["modules"])
    return cal, model_path


def main():
    cal, model_path = _calibrator()
    pset = cal.parameter_set

    def objective(theta):                       # minimise 1 - KGE
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return 1.0 - cal.score(theta).score

    spectrum = eigenspectrum(objective, pset)
    rdg = ridge(objective, pset, "log__tau_soil", "f_exf_shallow", n=25)
    os.remove(model_path)

    fig = plt.figure(figsize=(12.5, 3.8))
    ax0 = fig.add_subplot(1, 3, 1)
    ax1 = fig.add_subplot(1, 3, 2)
    ax2 = fig.add_subplot(1, 3, 3)
    spectrum.plot(axes=(ax0, ax1))
    rdg.plot(ax=ax2)
    ax2.set_title("ridge: soil timescale vs. exfiltration", fontsize=9)
    fig.suptitle("Parameter identifiability (Cannon, three parameters)",
                 fontsize=11)
    fig.savefig(OUT, dpi=125, bbox_inches="tight")
    print(f"wrote {os.path.relpath(OUT, ROOT)}")


if __name__ == "__main__":
    main()
