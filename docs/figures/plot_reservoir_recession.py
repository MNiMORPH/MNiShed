"""
Generate the power-law reservoir recession figure for the documentation.

Discharges a single reservoir from a common initial storage with no recharge for
several recession exponents ``b``, showing how ``b`` sets the recession shape:
``b = 1`` is the linear reservoir (a straight line in log-discharge), while
``b > 1`` gives a faster early and a slower late recession (a fatter tail).

Usage (from the repository root)::

    python docs/figures/plot_reservoir_recession.py

Writes ``docs/source/_static/reservoir_recession.png``.
"""

import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from mnished import Reservoir  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
OUT = os.path.join(ROOT, "docs", "source", "_static", "reservoir_recession.png")

TAU = 20.0       # recession coefficient (days)
H0 = 100.0       # initial storage (mm)
NDAYS = 150


def _recession(b):
    r = Reservoir(recession_coeff=TAU, H0=H0)
    r.recession_exponent = b
    r.recession_H_ref = H0
    q = np.empty(NDAYS)
    for t in range(NDAYS):
        r.discharge(1.0)
        q[t] = r.H_discharge
        r.H_excess = r.H_deficit = r.H_exfiltrated = 0.0
    return q


def main():
    t = np.arange(NDAYS)
    colors = {1.0: "#1f77b4", 1.5: "#2ca02c", 2.0: "#ff7f0e", 3.0: "#d62728"}
    fig, ax = plt.subplots(figsize=(5.8, 3.6))
    for b, c in colors.items():
        label = f"b = {b:g}" + (" (linear)" if b == 1.0 else "")
        ax.semilogy(t, _recession(b), color=c, lw=1.8, label=label)
    ax.set_xlabel("days since recession start")
    ax.set_ylabel("discharge (mm/day)")
    ax.set_ylim(1e-2, 10)
    ax.set_title(r"Power-law reservoir recession ($\tau$ = 20 d)", fontsize=10)
    ax.legend(fontsize=8, title="recession exponent")
    ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)
    fig.savefig(OUT, dpi=130, bbox_inches="tight")
    print(f"wrote {os.path.relpath(OUT, ROOT)}")


if __name__ == "__main__":
    main()
