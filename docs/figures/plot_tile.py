"""
Generate the multipath (threshold-activated) tile-drain figure for the docs.

Recedes a reservoir from high storage with and without a multipath drain. The
multipath drain is *storage-state dependent*: it adds a fast outflow
``max(0, H - H_thr) / tau_mp`` only while storage is above the drain elevation
``H_thr``, giving a steep early recession that reverts to the slow matrix
recession once storage falls below the threshold. (The other tile representation,
the fractional-bypass ``f_tile`` / ``tau_tile``, instead diverts a constant
fraction of inter-reservoir flow and so needs the full cascade to act.)

Usage (from the repository root)::

    python docs/figures/plot_tile.py

Writes ``docs/source/_static/tile_multipath.png``.
"""

import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from mnished import Reservoir  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
OUT = os.path.join(ROOT, "docs", "source", "_static", "tile_multipath.png")

TAU = 40.0
H0 = 100.0
H_THR = 40.0
NDAYS = 150


def _recession(**kw):
    r = Reservoir(recession_coeff=TAU, H0=H0, **kw)
    q, h = np.empty(NDAYS), np.empty(NDAYS)
    for t in range(NDAYS):
        h[t] = r.Hwater
        r.discharge(1.0)
        q[t] = r.H_discharge
        r.H_excess = r.H_deficit = r.H_exfiltrated = 0.0
    return q, h


def main():
    t = np.arange(NDAYS)
    q_plain, _ = _recession()
    q_mp, h_mp = _recession(multipath_threshold=H_THR, multipath_timescale=8.0)
    cross = int(np.argmax(h_mp < H_THR))      # day storage falls below threshold

    fig, ax = plt.subplots(figsize=(5.8, 3.6))
    ax.semilogy(t, q_plain, color="#888", lw=1.8, label="matrix recession only")
    ax.semilogy(t, q_mp, color="#1f6fb2", lw=1.9,
                label="+ multipath (threshold)")
    ax.axvline(cross, color="#c66", ls="--", lw=1)
    ax.text(cross + 2, 4, "storage drops\nbelow drain depth", fontsize=7.5,
            color="#a44")
    ax.set_xlabel("days since recession start")
    ax.set_ylabel("discharge (mm/day)")
    ax.set_ylim(1e-2, 20)
    ax.set_title("Threshold-activated (multipath) tile drainage", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)
    fig.savefig(OUT, dpi=130, bbox_inches="tight")
    print(f"wrote {os.path.relpath(OUT, ROOT)}  (threshold crossed day {cross})")


if __name__ == "__main__":
    main()
