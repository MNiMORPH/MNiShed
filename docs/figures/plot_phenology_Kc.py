"""
Generate the vegetation-phenology figure for the documentation.

Two stacked, x-aligned panels on the Cannon example forcing with the default
``phenology:`` settings: cumulative growing-degree-days (GDD) with the
``leafout``/``full`` thresholds on top, and the resulting vegetation coefficient
``Kc`` below — so the figure reads as cause (GDD crossing a threshold) to effect
(the Kc ramp). Illustrates the equations in ``model_description.rst``.

Usage (from the repository root)::

    python docs/figures/plot_phenology_Kc.py

Writes ``docs/source/_static/phenology_Kc.png``.
"""

import os
import warnings

import matplotlib
import numpy as np
import pandas as pd
import yaml

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402  (after Agg backend)

from mnished import Buckets  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
CANNON_CSV = os.path.join(ROOT, "examples", "cannon_forward", "CannonTestInput.csv")
OUT = os.path.join(ROOT, "docs", "source", "_static", "phenology_Kc.png")


def _model():
    cfg = {
        "timeseries": {"datafile": CANNON_CSV},
        "catchment": {"drainage_basin_area__km2": 3800,
                      "evapotranspiration_method": "ThornthwaiteChang2019",
                      "water_year_start_month": 10},
        "general": {"spin_up_cycles": 0},
        "reservoirs": {"recession_coefficients": [14, 500],
                       "exfiltration_fractions": [0.3, 1.0],
                       "maximum_effective_depths__mm": [float("inf"),
                                                        float("inf")]},
        "initial_conditions": {"water_reservoir_effective_depths__mm": [15, 400],
                               "snowpack__mm_SWE": 0},
        "snowmelt": {"PDD_melt_factor": 1.0},
        "modules": {"snowpack": True, "frozen_ground": False,
                    "rain_on_snow": True, "direct_runoff": False},
        "phenology": {"enabled": True},
    }
    path = os.path.join(HERE, "_phenology_cfg.yml")
    with open(path, "w") as f:
        yaml.safe_dump(cfg, f)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = Buckets()
        b.initialize(path, enforce_water_balance="global")
    os.remove(path)
    return b


def main():
    b = _model()
    p = b.phenology_params
    hd = b.hydrodata
    Tmean = 0.5 * (hd["Maximum Temperature [C]"].to_numpy(float)
                   + hd["Minimum Temperature [C]"].to_numpy(float))
    dates = pd.DatetimeIndex(hd["Date"])
    years, doy = dates.year.to_numpy(), dates.dayofyear.to_numpy()
    gdd_day = np.maximum(Tmean - p["base_temperature__C"], 0.0)
    GDD = np.empty_like(gdd_day)
    for y in np.unique(years):                 # reset each calendar year
        m = years == y
        GDD[m] = np.nancumsum(gdd_day[m])
    gddc = pd.Series(GDD, index=doy).groupby(level=0).mean()       # climatology
    kcc = pd.Series(np.asarray(b.phenology_Kc()),
                    index=doy).groupby(level=0).mean()
    d_leaf = int(gddc.index[gddc.values >= p["leafout_GDD"]][0])
    d_full = int(gddc.index[gddc.values >= p["full_canopy_GDD"]][0])

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.2, 4.6), sharex=True,
                                 gridspec_kw=dict(hspace=0.12))
    a1.plot(gddc.index, gddc.values, color="#c63", lw=2.0)
    for thr, lab in [(p["leafout_GDD"], r"$\mathrm{GDD}_\mathrm{leaf}$"),
                     (p["full_canopy_GDD"], r"$\mathrm{GDD}_\mathrm{full}$")]:
        a1.axhline(thr, ls=":", color="gray", lw=1)
        a1.text(6, thr + 15, lab, fontsize=8, color="gray")
    a1.set_ylabel("cumulative GDD\n(°C·day)")
    a1.set_ylim(0, max(gddc.values) * 1.05)
    a1.set_title("GDD vegetation-phenology coefficient (Cannon defaults)",
                 fontsize=10)

    a2.plot(kcc.index, kcc.values, color="#2a7", lw=2.2)
    a2.axhline(p["dormant_Kc"], ls=":", color="gray", lw=1)
    a2.axhline(p["full_Kc"], ls=":", color="gray", lw=1)
    a2.text(6, p["dormant_Kc"] + 0.015, r"$K_{c,\mathrm{dormant}}$", fontsize=8,
            color="gray")
    a2.text(6, p["full_Kc"] - 0.05, r"$K_{c,\mathrm{full}}$", fontsize=8,
            color="gray")
    a2.axvspan(p["senescence_start_doy"], p["senescence_end_doy"],
               color="#d98", alpha=0.18)
    a2.text((p["senescence_start_doy"] + p["senescence_end_doy"]) / 2, 0.55,
            "senescence\n(day-of-year)", ha="center", fontsize=7.5, color="#a63")
    a2.set_ylabel(r"coefficient $K_c$")
    a2.set_ylim(0.35, 1.05)
    a2.set_xlabel("day of year")

    for d in (d_leaf, d_full):                 # link GDD crossings to the Kc ramp
        for ax in (a1, a2):
            ax.axvline(d, ls="--", color="#555", lw=0.8, alpha=0.6)
    a2.annotate("leaf-out", xy=(d_leaf, p["dormant_Kc"]),
                xytext=(d_leaf + 18, 0.62), fontsize=7.5, color="#262",
                arrowprops=dict(arrowstyle="->", color="#262", lw=1))
    a1.set_xlim(1, 365)
    a2.set_xticks([1, 60, 121, 182, 244, 305, 365])
    a2.set_xticklabels(["Jan", "Mar", "May", "Jul", "Sep", "Nov", ""])
    fig.savefig(OUT, dpi=130, bbox_inches="tight")
    print(f"wrote {os.path.relpath(OUT, ROOT)}  "
          f"(leaf-out day {d_leaf}, full-canopy day {d_full})")


if __name__ == "__main__":
    main()
