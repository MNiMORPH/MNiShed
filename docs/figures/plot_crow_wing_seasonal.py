"""
Generate the Crow Wing seasonal mass-balance figure for the documentation.

Runs the Crow Wing example (``examples/crow_wing/``) at an SCE-UA-calibrated
parameter set over a representative decade (2001-2010) with ``store_fluxes`` and
plots the per-season discharge split by source (fast/event, slow/baseflow, lake
outlet) against observed Q — the real-basin version of the seasonal mass-balance
diagnostic. The parameters below are a calibration result (KGE_logKGE ~0.68 over
the decade); re-run the example's ``calibrate.py`` to refit.

Usage (from the repository root)::

    python docs/figures/plot_crow_wing_seasonal.py

Writes ``docs/source/_static/crow_wing_seasonal.png``.
"""

import os
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from mnished import Calibrator, SeasonalMassBalance, run_and_score  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
EXAMPLE = os.path.join(ROOT, "examples", "crow_wing")
OUT = os.path.join(ROOT, "docs", "source", "_static", "crow_wing_seasonal.png")

# SCE-UA best-fit (calibration coordinates; log__ names are log10). Regenerate
# with `python calibrate.py sceua` in examples/crow_wing/.
THETA = {
    "PDD_melt_factor": 0.8385,
    "et_scale": 0.9313,
    "log__recession_coeff_soil": 0.9527,
    "log__recession_coeff_gw": 3.4310,
    "f_exfil_soil": 0.3074,
    "log__recession_coeff_lake": 3.9846,
    "log__H_sill_lake": 3.0082,
    "f_route_lake": 0.9006,
    "log__fdd_threshold": 1.8869,
    "leafout_GDD": 225.04,
}
START, END = "2001-01-01", "2010-12-31"


def main():
    cwd = os.getcwd()
    os.chdir(EXAMPLE)                       # relative config_template / forcing
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cal = Calibrator.from_yaml("params.yml")
            d = cal.driver
            res = run_and_score(
                d["config_template"], modules=cal.modules,
                routing_N=d["routing_N"], spin_up_cycles=d["spin_up_cycles"],
                metric=d["metric"],
                enforce_water_balance=d.get("enforce_water_balance", "none"),
                start=START, end=END, store_fluxes=True, **cal.run_kwargs(THETA))
    finally:
        os.chdir(cwd)
    st = SeasonalMassBalance(res.buckets, start=START, end=END).seasonal_table()

    seasons = list(st.index)
    x = range(len(seasons))
    fast, slow, lake = st["fast"].values, st["slow"].values, st["lake"].values
    fig, ax = plt.subplots(figsize=(5.8, 3.4))
    ax.bar(x, fast, color="#5aa0d0", label="fast (event)")
    ax.bar(x, slow, bottom=fast, color="#2e6b9e", label="slow (baseflow)")
    ax.bar(x, lake, bottom=fast + slow, color="#7fc97f", label="lake outlet")
    ax.plot(x, st["obs"].values, "o-", color="k", lw=1.5, ms=5, label="observed Q")
    ax.set_xticks(list(x))
    ax.set_xticklabels(seasons)
    ax.set_ylabel("discharge (mm/day)")
    ax.set_title("Crow Wing seasonal mass balance, 2001-2010 (modeled Q by source)",
                 fontsize=9.5)
    ax.legend(fontsize=7.5, loc="upper left", framealpha=0.9)
    fig.savefig(OUT, dpi=130, bbox_inches="tight")
    print(f"wrote {os.path.relpath(OUT, ROOT)}  (decade score {res.score:.3f})")


if __name__ == "__main__":
    main()
