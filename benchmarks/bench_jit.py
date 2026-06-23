"""
Numba-JIT performance benchmark for MNiShed's daily time loop.

Run from the repository root, **in an environment with Numba** (the ``jit``
extra) so the JIT path is actually exercised::

    pip install '.[jit]'
    python benchmarks/bench_jit.py
    python benchmarks/plot_jit.py

It times one forward ``Buckets.run()`` over a sweep of record lengths,
comparing the Numba JIT path against the *same code* forced to pure Python
(via ``mnished.mnished._numba_available``).  Each model is a synthetic
3-reservoir watershed with daily forcing; only the record length varies, so
the timings isolate the per-step cost of the time loop.

Results are printed to stdout and saved to ``benchmarks/results/`` as
``{short_commit}_{UTC_timestamp}.txt``.  The header records the commit,
working-tree state, and hardware/software so runs are reproducible and
comparable; ``plot_jit.py`` reads the most recent results file.
"""

import datetime
import os
import platform
import subprocess
import sys
import time
import warnings

import numpy as np
import pandas as pd

import mnished.mnished as M
from mnished import Buckets

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, "results")

# Record lengths to sweep (years × 365 days), 3 reservoirs throughout.
RECORD_YEARS = (2, 3, 5, 8, 12, 20, 30, 45, 65, 90, 130, 180)


def _git():
    def run(*args):
        try:
            return subprocess.check_output(
                ["git", *args], cwd=HERE,
                stderr=subprocess.DEVNULL).decode().strip()
        except Exception:
            return "?"
    full = run("rev-parse", "HEAD")
    short = full[:7] if full != "?" else "nogit"
    branch = run("rev-parse", "--abbrev-ref", "HEAD")
    dirty = run("status", "--porcelain") != ""
    return full, short, branch, dirty


def _hardware():
    cpu = platform.processor() or "unknown"
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("model name"):
                    cpu = line.split(":", 1)[1].strip()
                    break
    except Exception:
        pass
    mem = "?"
    try:
        with open("/proc/meminfo") as f:
            mem = f"{int(f.readline().split()[1]) / 1048576:.1f} GiB"
    except Exception:
        pass
    return cpu, mem


def make_config(n_days, n_res, tag):
    csv = os.path.join(RESULTS_DIR, f"_synth_{tag}.csv")
    yml = os.path.join(RESULTS_DIR, f"_synth_{tag}.yml")
    rng = np.random.default_rng(7)
    dates = pd.date_range("1700-01-01", periods=n_days, freq="D")
    doy = dates.dayofyear.values
    temp = 8 + 15 * np.sin(2 * np.pi * (doy - 100) / 365) + rng.normal(0, 3, n_days)
    precip = np.clip(rng.gamma(0.4, 6.0, n_days), 0, None)
    et = np.clip(2.5 + 2.0 * np.sin(2 * np.pi * (doy - 150) / 365), 0, None)
    disch = np.clip(5 + 3 * np.sin(2 * np.pi * (doy - 120) / 365)
                    + rng.normal(0, 1, n_days), 0.1, None)
    pd.DataFrame({"Date": dates, "Precipitation [mm/day]": precip,
                  "Discharge [m^3/s]": disch, "Mean Temperature [C]": temp,
                  "Evapotranspiration [mm/day]": et}).to_csv(csv, index=False)
    rc = [16.0, 200.0, 3650.0][:n_res]
    ef = [0.6] * (n_res - 1) + [1.0]
    hm = ", ".join([".inf"] * n_res)
    with open(yml, "w") as f:
        f.write(f"""
timeseries: {{datafile: {csv}}}
initial_conditions:
    water_reservoir_effective_depths__mm: {[10.0] * n_res}
    snowpack__mm_SWE: 0
catchment:
    drainage_basin_area__km2: 3800.0
    evapotranspiration_method: datafile
    water_year_start_month: 10
general: {{spin_up_cycles: 1, enforce_water_balance: none}}
reservoirs:
    recession_coefficients: {rc}
    exfiltration_fractions: {ef}
    maximum_effective_depths__mm: [{hm}]
snowmelt: {{PDD_melt_factor: 2.0, fdd_threshold: .inf}}
modules: {{snowpack: true, frozen_ground: false}}
""")
    return yml


def init_from(yml):
    b = Buckets()
    b.initialize(yml)
    return b


def time_run(b, reps):
    t0 = time.perf_counter()
    for _ in range(reps):
        b.run()
    return (time.perf_counter() - t0) / reps


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    full, short, branch, dirty = _git()
    cpu, mem = _hardware()
    import numba
    if not M._numba_available:
        print("WARNING: Numba is unavailable, so the 'JIT' timings will fall "
              "back to pure Python and the speedup will read as ~1x. Run in a "
              "`pip install '.[jit]'` environment.", file=sys.stderr)

    rows = []
    n_first = None
    t_first = None
    for yr in RECORD_YEARS:
        n = 365 * yr
        yml = make_config(n, 3, n)
        M._numba_available = True
        b = init_from(yml)
        if t_first is None:
            t_first = time_run(b, 1)   # first run absorbs the one-time compile
            n_first = n
        else:
            b.run()                    # warm
        t_jit = time_run(b, 8)
        M._numba_available = False
        b2 = init_from(yml)
        reps = 2 if n > 30000 else (3 if n > 10000 else 5)
        t_py = time_run(b2, reps)
        M._numba_available = True
        rows.append((n, t_py * 1e3, t_jit * 1e3))

    Ns = np.array([r[0] for r in rows], float)
    t_py = np.array([r[1] for r in rows])
    t_jit = np.array([r[2] for r in rows])
    m_py, c_py = np.polyfit(Ns, t_py, 1)
    m_jit, c_jit = np.polyfit(Ns, t_jit, 1)

    def r2(x, y, m, c):
        yp = m * x + c
        return 1 - np.sum((y - yp) ** 2) / np.sum((y - y.mean()) ** 2)

    ts = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = os.path.join(RESULTS_DIR, f"{short}_{ts}.txt")
    lines = [
        "# MNiShed Numba-JIT benchmark",
        f"timestamp : {ts}",
        f"commit    : {full}{' (dirty working tree)' if dirty else ''}",
        f"branch    : {branch}",
        f"cpu       : {cpu}",
        f"cpu cores : {os.cpu_count()}",
        f"memory    : {mem}",
        f"os        : {platform.platform()}",
        f"python    : {platform.python_version()}  "
        f"numpy {np.__version__}  numba {numba.__version__}",
        "",
        f"fit pure-Python : T[ms] = {m_py:.5g} * N + {c_py:.4g}  "
        f"(R2={r2(Ns, t_py, m_py, c_py):.5f})  per-step {m_py * 1e3:.2f} us/day",
        f"fit Numba JIT   : T[ms] = {m_jit:.5g} * N + {c_jit:.4g}  "
        f"(R2={r2(Ns, t_jit, m_jit, c_jit):.5f})  per-step {m_jit * 1e3:.3f} us/day",
        f"asymptotic per-step speedup m_py/m_jit = {m_py / m_jit:.0f}x",
        f"one-time JIT compile (first run, {n_first} days) ~ {t_first:.1f} s",
        "",
        "# DATA  n_days  t_pyPython_ms  t_jit_ms  speedup",
    ]
    for n, tp, tj in rows:
        lines.append(f"{int(n):8d}  {tp:13.3f}  {tj:9.4f}  {tp / tj:7.1f}")
    text = "\n".join(lines) + "\n"
    with open(out, "w") as f:
        f.write(text)
    print(text)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
