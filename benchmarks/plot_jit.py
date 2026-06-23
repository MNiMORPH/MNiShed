"""
Generate the JIT-scaling figure from a bench_jit results file.

Usage (from the repository root)::

    python benchmarks/plot_jit.py [results_file]

With no argument, uses the most recent file in ``benchmarks/results/``.
Writes ``docs/source/_static/jit_scaling.png``.
"""

import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "results")
STATIC_DIR = os.path.join(HERE, "..", "docs", "source", "_static")


def _latest_results():
    files = [f for f in glob.glob(os.path.join(RESULTS_DIR, "*.txt"))]
    if not files:
        sys.exit("No results files in benchmarks/results/ — run bench_jit.py first.")
    return max(files, key=os.path.getmtime)


def _parse(path):
    """Return (Ns, t_py_ms, t_jit_ms) arrays from the '# DATA' table."""
    rows = []
    in_data = False
    with open(path) as f:
        for line in f:
            if line.startswith("# DATA"):
                in_data = True
                continue
            if in_data and line.strip():
                n, tp, tj = line.split()[:3]
                rows.append((float(n), float(tp), float(tj)))
    arr = np.array(rows)
    return arr[:, 0], arr[:, 1], arr[:, 2]


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else _latest_results()
    Ns, t_py, t_jit = _parse(path)
    yrs = Ns / 365.0
    m_py, c_py = np.polyfit(Ns, t_py, 1)
    m_jit, c_jit = np.polyfit(Ns, t_jit, 1)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))

    # Panel 1: both runtimes are linear (mx + b); twin axes because they
    # differ by ~300x.
    ax1.plot(yrs, t_py, "o", color="C3", label="pure-Python")
    ax1.plot(yrs, m_py * Ns + c_py, "-", color="C3", lw=1)
    ax1.set_xlabel("record length [years]")
    ax1.set_ylabel("pure-Python runtime [ms]", color="C3")
    ax1.tick_params(axis="y", labelcolor="C3")
    ax1.set_title("Runtime is linear in record length")
    ax1.grid(alpha=0.3)
    axr = ax1.twinx()
    axr.plot(yrs, t_jit, "s", color="C0", label="Numba JIT")
    axr.plot(yrs, m_jit * Ns + c_jit, "-", color="C0", lw=1)
    axr.set_ylabel("Numba JIT runtime [ms]", color="C0")
    axr.tick_params(axis="y", labelcolor="C0")
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = axr.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc="upper left")

    # Panel 2: speedup rises to the per-step ratio m_py / m_jit.
    ax2.plot(yrs, t_py / t_jit, "o-", color="k")
    ax2.axhline(m_py / m_jit, ls="--", color="0.5",
                label=r"asymptote $m_{\mathrm{py}}/m_{\mathrm{jit}}$ = "
                      f"{m_py / m_jit:.0f}×")
    ax2.set_xlabel("record length [years]")
    ax2.set_ylabel("speedup (×)")
    ax2.set_title("Speedup rises to the per-step ratio")
    ax2.legend(loc="lower right")
    ax2.grid(alpha=0.3)

    fig.suptitle("MNiShed Numba JIT scaling — daily time loop, 3 reservoirs")
    fig.tight_layout()
    os.makedirs(STATIC_DIR, exist_ok=True)
    out = os.path.normpath(os.path.join(STATIC_DIR, "jit_scaling.png"))
    fig.savefig(out, dpi=130, bbox_inches="tight")
    print(f"parsed {os.path.basename(path)}")
    print(f"pure-Python {m_py * 1e3:.1f} us/day | JIT {m_jit * 1e3:.3f} us/day "
          f"| asymptotic {m_py / m_jit:.0f}x")
    print(f"saved {out}")


if __name__ == "__main__":
    main()
