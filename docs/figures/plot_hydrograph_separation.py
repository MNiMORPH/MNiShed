"""
Generate the hydrograph-separation figure for the documentation.

Fits ``HydrographSeparation`` to the Cannon example discharge and plots the
record decomposed into its fitted timescale components — a fast (event) reservoir
above a slow (baseflow) store — the separation used to seed reservoir timescales
and initial conditions.

Usage (from the repository root)::

    python docs/figures/plot_hydrograph_separation.py

Writes ``docs/source/_static/hydrograph_separation.png``.
"""

import os
import warnings

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

from mnished import HydrographSeparation  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
CANNON_CSV = os.path.join(ROOT, "examples", "cannon_forward", "CannonTestInput.csv")
OUT = os.path.join(ROOT, "docs", "source", "_static", "hydrograph_separation.png")


def main():
    df = pd.read_csv(CANNON_CSV)
    date = pd.to_datetime(df["Date"], format="mixed", errors="coerce")
    qcol = [c for c in df.columns if "Discharge" in c and "m^3" in c][0]
    Q = pd.to_numeric(df[qcol], errors="coerce").to_numpy()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        hs = HydrographSeparation(Q).fit()

    fast = np.asarray(hs._Q_components_fast)[0]
    slow = np.asarray(hs._Q_residual)
    tau_fast = float(np.min(hs.tau))

    fig, ax = plt.subplots(figsize=(6.6, 3.3))
    ax.fill_between(date, 0, slow, color="#7fb37f", label="slow / baseflow")
    ax.fill_between(date, slow, slow + np.clip(fast, 0, None), color="#9ecae1",
                    label=f"fast (event, τ ≈ {tau_fast:.0f} d)")
    ax.plot(date, Q, color="k", lw=0.6, alpha=0.8, label="total discharge")
    ax.set_ylabel("discharge (m³/s)")
    ax.set_title("Hydrograph separation into timescale components (Cannon)",
                 fontsize=10)
    ax.legend(fontsize=7.5, loc="upper right", framealpha=0.9)
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    fig.savefig(OUT, dpi=130, bbox_inches="tight")
    print(f"wrote {os.path.relpath(OUT, ROOT)}  (tau_fast {tau_fast:.1f} d)")


if __name__ == "__main__":
    main()
