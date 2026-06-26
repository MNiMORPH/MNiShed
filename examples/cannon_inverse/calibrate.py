#!/usr/bin/env python3
"""Generic in-process calibration — any basin, from params.yml alone.

The model run is declared in config, not in code: each parameter's `target`
field (see params.yml) says where it maps in run_and_score, so this one runner
calibrates any basin with no per-basin Python. It builds the model once with
mnished.ScoringModel and reuses it for every evaluation.

    conda activate mnished-jit          # numba JIT + spotpy
    python calibrate.py sceua [reps]            # best-fit (SCE-UA)
    python calibrate.py dream [reps] [iid|ar1]  # posterior / UQ (DREAM)

Phase 1: flat single-cascade parameters (name / name[i] targets). Sub-catchment
and lake targets (nested) are a planned extension (MNiMORPH/MNiShed#20).
"""

import re
import sys
import time

import numpy as np
import yaml
import spotpy

from mnished import ParameterSet, ScoringModel, log_flow_residual_terms

with open('params.yml') as f:
    _CFG = yaml.safe_load(f)
P, DRV, MOD = _CFG['parameters'], _CFG['driver'], _CFG.get('modules', {})
PSET = ParameterSet.from_params_yml(P)
NAMES = PSET.names
START, END = DRV.get('decade_start'), DRV.get('decade_end')


def build_kwargs(theta):
    """Assemble run_and_score keywords from each parameter's `target`.

    ``name[i]`` -> a list keyword at position i (grouped across parameters);
    bare ``name`` -> a scalar keyword; a ``log__`` prefix applies 10**.
    """
    lists, kw = {}, {}
    for name, spec in P.items():
        target = spec.get('target')
        if not target:
            continue
        val = theta.get(name, spec['fixed'])
        if name.startswith('log__'):
            val = 10.0 ** val
        m = re.fullmatch(r'(\w+)\[(\d+)\]', target)
        if m:
            lists.setdefault(m.group(1), {})[int(m.group(2))] = val
        else:
            kw[target] = val
    for key, idx in lists.items():
        kw[key] = [idx[i] for i in range(max(idx) + 1)]
    return kw


def _ar1_loglike(r, phi):
    n = len(r)
    one_m = max(1.0 - phi * phi, 1e-12)
    e = np.empty(n)
    e[0] = r[0] * np.sqrt(one_m)
    e[1:] = r[1:] - phi * r[:-1]
    sse = float(np.dot(e, e))
    return -1e300 if (not np.isfinite(sse) or sse <= 0) else \
        0.5 * np.log(one_m) - 0.5 * n * np.log(sse)


class Calibration:
    """SPOTPY setup driven entirely by params.yml targets."""

    def __init__(self, likelihood=None):
        self.likelihood = likelihood          # None=optimise KGE; 'iid'/'ar1'=DREAM
        self.evals = 0
        self.model = ScoringModel(
            DRV['config_template'],
            enforce_water_balance=DRV.get('enforce_water_balance', 'water-year'))
        self.params = [spotpy.parameter.Uniform(p.name, p.lower, p.upper,
                                                optguess=p.value) for p in PSET]
        if likelihood == 'ar1':
            self.params.append(spotpy.parameter.Uniform('likelihood_phi1',
                                                        0.0, 0.99, optguess=0.9))
        if likelihood:                        # observed log-flows (fixed): compute once
            ref = self._run({n: P[n]['initial'] for n in NAMES})
            self.log_obs = log_flow_residual_terms(
                ref, start=START, end=END)['log_obs'].to_numpy()

    def _run(self, theta):
        return self.model.score(modules=MOD, routing_N=DRV['routing_N'],
                                spin_up_cycles=DRV['spin_up_cycles'],
                                metric=DRV['metric'], start=START, end=END,
                                **build_kwargs(theta))

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        self.evals += 1
        if self.likelihood == 'ar1':
            self._phi = float(vector[len(NAMES)])
            vector = vector[:len(NAMES)]
        theta = dict(zip(NAMES, vector))
        try:
            res = self._run(theta)
        except Exception:
            return list(self.log_obs - 10.0) if self.likelihood else [-10.0]
        if self.likelihood:
            return list(log_flow_residual_terms(
                res, start=START, end=END)['log_mod'].to_numpy())
        return [res.score if np.isfinite(res.score) else -10.0]

    def evaluation(self):
        return list(self.log_obs) if self.likelihood else [1.0]

    def objectivefunction(self, simulation, evaluation):
        if self.likelihood == 'ar1':
            return _ar1_loglike(np.asarray(simulation) - np.asarray(evaluation),
                                self._phi)
        if self.likelihood == 'iid':
            return spotpy.likelihoods.gaussianLikelihoodMeasErrorOut(
                evaluation, simulation)
        return 1.0 - simulation[0]            # minimise 1 - KGE


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'sceua'
    reps = int(sys.argv[2]) if len(sys.argv) > 2 else 4000
    t0 = time.time()
    if mode == 'sceua':
        setup = Calibration()
        sampler = spotpy.algorithms.sceua(setup, dbname='calib_sceua',
                                          dbformat='ram')
        sampler.sample(reps)
        dt = time.time() - t0
        print(f"\nsceua: {setup.evals} evals, wall={dt:.1f}s "
              f"({dt / max(setup.evals, 1) * 1000:.0f} ms/eval)")
        print(f"best {DRV['metric']} = "
              f"{1.0 - sampler.status.objectivefunction_min:.4f}")
    elif mode == 'dream':
        like = sys.argv[3] if len(sys.argv) > 3 else 'iid'
        setup = Calibration(likelihood=like)
        sampler = spotpy.algorithms.dream(setup, dbname=f'calib_dream_{like}',
                                          dbformat='csv', parallel='seq')
        sampler.sample(reps, nChains=8)
        dt = time.time() - t0
        print(f"\ndream-{like}: {setup.evals} evals, wall={dt:.1f}s "
              f"({dt / max(setup.evals, 1) * 1000:.0f} ms/eval)")
    else:
        raise SystemExit("mode must be 'sceua' or 'dream'")


if __name__ == '__main__':
    main()
