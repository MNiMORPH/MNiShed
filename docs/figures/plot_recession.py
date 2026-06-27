"""
Generate the Brutsaert-Nieber recession-analysis figure for the documentation.

Fits the Brutsaert-Nieber recession model to the Cannon example discharge and
uses the built-in plot: the recession "cloud" of -dQ/dt against Q in log-log
space with the fitted power law, from which a reservoir recession exponent is
read.

Usage (from the repository root)::

    python docs/figures/plot_recession.py

Writes ``docs/source/_static/recession_bn.png``.
"""

import os
import warnings

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from mnished import BrutsaertNieber  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
CANNON_CSV = os.path.join(ROOT, "examples", "cannon_forward", "CannonTestInput.csv")
OUT = os.path.join(ROOT, "docs", "source", "_static", "recession_bn.png")


def main():
    df = pd.read_csv(CANNON_CSV)
    qcol = [c for c in df.columns if "Discharge" in c and "m^3" in c][0]
    Q = pd.to_numeric(df[qcol], errors="coerce").to_numpy()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bn = BrutsaertNieber(Q, min_recession_days=3).fit()

    fig, ax = plt.subplots(figsize=(5.4, 4.0))
    bn.plot(ax=ax)
    ax.set_title("Brutsaert-Nieber recession analysis (Cannon discharge)",
                 fontsize=10)
    fig.savefig(OUT, dpi=130, bbox_inches="tight")
    print(f"wrote {os.path.relpath(OUT, ROOT)}  (a={bn.a_:.3g}, n={bn.n_:.3g})")


if __name__ == "__main__":
    main()
