#!/usr/bin/env python3
"""Generic in-process calibration — driven by params.yml, no per-basin Python.

A thin SPOTPY adapter over ``mnished.Calibrator`` (the standard config-driven
model setup): Calibrator reads each parameter's ``target`` from params.yml,
builds the model once (ScoringModel), and scores a parameter set; this script
wires it to SPOTPY's SCE-UA (best-fit) or DREAM (UQ).

    conda activate mnished-jit                 # numba JIT + spotpy
    python calibrate.py sceua [reps]           # best-fit (SCE-UA)
    python calibrate.py dream [reps] [iid|ar1] # posterior / UQ (DREAM)

Targets may be flat (``name`` / ``name[i]``) or nested sub-catchment / lake
(``sub_catchments[I].key[j]``, with I a shared-index list like ``0,1``). See
``mnished.Calibrator`` and MNiMORPH/MNiShed#20.
"""

import sys
import time

import numpy as np
import spotpy

from mnished import Calibrator, log_flow_residual_terms

CAL = Calibrator.from_yaml('params.yml')
NAMES = CAL.names
WINDOWS = CAL.windows                   # one window, or several (driver `decades:`)


def _ar1_loglike(r, phi):
    n = len(r)
    one_m = max(1.0 - phi * phi, 1e-12)
    e = np.empty(n)
    e[0] = r[0] * np.sqrt(one_m)
    e[1:] = r[1:] - phi * r[:-1]
    sse = float(np.dot(e, e))
    return -1e300 if (not np.isfinite(sse) or sse <= 0) else \
        0.5 * np.log(one_m) - 0.5 * n * np.log(sse)


class Setup:
    """SPOTPY setup: model = Calibrator.score; objective = KGE or a likelihood."""

    def __init__(self, likelihood=None):
        self.likelihood = likelihood          # None=optimise KGE; 'iid'/'ar1'=DREAM
        self.evals = 0
        self.params = [spotpy.parameter.Uniform(p.name, p.lower, p.upper,
                                                optguess=p.value)
                       for p in CAL.parameter_set]
        if likelihood == 'ar1':
            self.params.append(spotpy.parameter.Uniform('likelihood_phi1',
                                                        0.0, 0.99, optguess=0.9))
        if likelihood:                        # observed log-flows (fixed): compute once
            refs = CAL.score_windows(
                {p.name: p.value for p in CAL.parameter_set})
            parts = [log_flow_residual_terms(r, start=w['start'],
                                             end=w['end'])['log_obs'].to_numpy()
                     for r, w in zip(refs, WINDOWS)]
            self.log_obs = np.concatenate(parts)
            self._window_lens = [len(p) for p in parts]   # for per-window AR(1)

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        self.evals += 1
        if self.likelihood == 'ar1':
            self._phi = float(vector[len(NAMES)])
            vector = vector[:len(NAMES)]
        try:
            results = CAL.score_windows(dict(zip(NAMES, vector)))
        except Exception:
            return list(self.log_obs - 10.0) if self.likelihood else [-10.0]
        if self.likelihood:                   # residuals concatenated across windows
            return list(np.concatenate([
                log_flow_residual_terms(r, start=w['start'], end=w['end'])['log_mod'].to_numpy()
                for r, w in zip(results, WINDOWS)]))
        scores = [r.score if np.isfinite(r.score) else -10.0 for r in results]
        return [float(np.mean(scores))]       # mean KGE over windows

    def evaluation(self):
        return list(self.log_obs) if self.likelihood else [1.0]

    def objectivefunction(self, simulation, evaluation):
        if self.likelihood == 'ar1':
            # Whiten each window's residuals on their own and sum the
            # log-likelihoods: windows are disjoint in time, so the AR(1) lag
            # must not cross a boundary (that would invent a correlation
            # between residuals decades apart).
            r = np.asarray(simulation) - np.asarray(evaluation)
            ll, i = 0.0, 0
            for n in self._window_lens:
                ll += _ar1_loglike(r[i:i + n], self._phi)
                i += n
            return ll
        if self.likelihood == 'iid':
            return spotpy.likelihoods.gaussianLikelihoodMeasErrorOut(
                evaluation, simulation)
        return 1.0 - simulation[0]            # minimise 1 - KGE


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'sceua'
    reps = int(sys.argv[2]) if len(sys.argv) > 2 else 4000
    t0 = time.time()
    if mode == 'sceua':
        setup = Setup()
        sampler = spotpy.algorithms.sceua(setup, dbname='calib_sceua',
                                          dbformat='ram')
        sampler.sample(reps)
        dt = time.time() - t0
        print(f"\nsceua: {setup.evals} evals, wall={dt:.1f}s "
              f"({dt / max(setup.evals, 1) * 1000:.0f} ms/eval)")
        print(f"best {CAL.driver['metric']} = "
              f"{1.0 - sampler.status.objectivefunction_min:.4f}")
    elif mode == 'dream':
        like = sys.argv[3] if len(sys.argv) > 3 else 'iid'
        setup = Setup(likelihood=like)
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
