"""
Generate the PDM saturation-excess figure for the documentation.

Plots the probability-distributed-model saturated fraction
``f_sat = 1 - exp(-H / H0)`` against shallow-reservoir storage for a few
characteristic depths ``H0``: that fraction of incoming water is shed directly
as saturation-excess runoff, so a smaller ``H0`` saturates (and spills) sooner.

Usage (from the repository root)::

    python docs/figures/plot_pdm.py

Writes ``docs/source/_static/pdm_saturation.png``.
"""

import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
OUT = os.path.join(ROOT, "docs", "source", "_static", "pdm_saturation.png")


def main():
    H = np.linspace(0, 150, 400)
    fig, ax = plt.subplots(figsize=(5.4, 3.4))
    for H0, c in [(20, "#d62728"), (50, "#ff7f0e"), (100, "#1f77b4")]:
        ax.plot(H, 1.0 - np.exp(-H / H0), color=c, lw=1.9,
                label=f"H0 = {H0} mm")
    ax.set_xlabel("shallow-reservoir storage H (mm)")
    ax.set_ylabel(r"saturated fraction $f_\mathrm{sat}$")
    ax.set_ylim(0, 1)
    ax.set_title(r"PDM saturation excess: $f_\mathrm{sat}=1-e^{-H/H_0}$",
                 fontsize=10)
    ax.legend(fontsize=8, title="characteristic depth")
    ax.grid(True, ls=":", lw=0.4, alpha=0.5)
    fig.savefig(OUT, dpi=130, bbox_inches="tight")
    print(f"wrote {os.path.relpath(OUT, ROOT)}")


if __name__ == "__main__":
    main()
