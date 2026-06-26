#!/usr/bin/env python3
"""In-process Bayesian (DREAM) calibration via SPOTPY — no Dakota fork.

Runs the DREAM algorithm (Vrugt) with the model called **in-process** (warm
JIT, no per-eval Python startup). Reuses the same observed log-flow data
(``cannon_logq_obs.dat``) and modelled-log-flow mapping
(``driver_bayes.modeled_log_flows``) as the Dakota path, so the likelihood is
identical to that path's intent — only the evaluation interface changes.

Two likelihoods:
  * ``iid`` — Gaussian on the log-flow residuals, errors independent
    (SPOTPY ``gaussianLikelihoodMeasErrorOut``; the measurement error is
    marginalised out). Treats the 1096 daily residuals as independent.
  * ``ar1`` — Gaussian with a first-order autoregressive error model: the
    residuals are whitened with an AR(1) coefficient ``phi1`` (sampled by
    DREAM) before the Gaussian, which accounts for the strong day-to-day
    autocorrelation streamflow residuals carry. This is the autocorrelation
    term of the Schoups & Vrugt (2010) generalized likelihood.

    python make_bayes_data.py            # writes cannon_logq_obs.dat (once)
    python run_spotpy.py [reps] [iid|ar1]

Run in an environment with the Numba JIT active (here: mnished-jit). Serial
('seq') is fastest: each eval is ~15-20 ms but returns a 1096-value vector, so
multiprocessing loses more to inter-process pickling than it saves on compute.
"""

import sys
import time

import numpy as np
import yaml
import spotpy

from mnished import ParameterSet
from driver_bayes import modeled_log_flows   # validated param -> log_mod mapping

with open('params.yml') as f:
    _cfg = yaml.safe_load(f)
P = _cfg['parameters']
PSET = ParameterSet.from_params_yml(P)            # uniform priors over bounds
LOG_OBS = np.loadtxt('cannon_logq_obs.dat')       # fixed observed log-flows


def _ar1_marginal_loglike(r, phi):
    """AR(1) Gaussian log-likelihood with the noise scale marginalised out.

    Whitens the residuals ``r`` with AR(1) coefficient ``phi`` and integrates
    sigma under a Jeffreys prior, giving ``0.5*log(1-phi^2) - (N/2)*log(SSE)``
    up to a constant. ``phi -> 0`` recovers the iid case.
    """
    n = len(r)
    one_m = max(1.0 - phi * phi, 1e-12)
    e = np.empty(n)
    e[0] = r[0] * np.sqrt(one_m)
    e[1:] = r[1:] - phi * r[:-1]
    sse = float(np.dot(e, e))
    if not np.isfinite(sse) or sse <= 0:
        return -1e300
    return 0.5 * np.log(one_m) - 0.5 * n * np.log(sse)


class CannonDREAM:
    """SPOTPY setup: model = run_and_score; likelihood = iid or AR(1) Gaussian."""

    def __init__(self, likelihood='iid'):
        self.likelihood = likelihood
        self.names = PSET.names
        self.params = [spotpy.parameter.Uniform(p.name, p.lower, p.upper,
                                                optguess=p.value)
                       for p in PSET]
        if likelihood == 'ar1':
            # AR(1) coefficient as an extra sampled (nuisance) parameter
            self.params.append(spotpy.parameter.Uniform(
                'likelihood_phi1', 0.0, 0.99, optguess=0.9))

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        if self.likelihood == 'ar1':
            self._phi = float(vector[-1])         # read in objectivefunction
            vec = vector[:len(self.names)]
        else:
            vec = vector
        theta = dict(zip(self.names, vec))
        get = lambda n: theta.get(n, P[n]['fixed'])   # noqa: E731
        try:
            return list(modeled_log_flows(get))
        except Exception:
            return list(LOG_OBS - 10.0)           # far-off -> low likelihood

    def evaluation(self):
        return list(LOG_OBS)

    def objectivefunction(self, simulation, evaluation):
        if self.likelihood == 'ar1':
            r = np.asarray(simulation) - np.asarray(evaluation)
            return _ar1_marginal_loglike(r, self._phi)
        return spotpy.likelihoods.gaussianLikelihoodMeasErrorOut(
            evaluation, simulation)


def main():
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
    like = sys.argv[2] if len(sys.argv) > 2 else 'iid'
    t0 = time.time()
    sampler = spotpy.algorithms.dream(
        CannonDREAM(like), dbname=f'cannon_dream_{like}', dbformat='csv',
        parallel='seq')
    sampler.sample(reps, nChains=8)
    dt = time.time() - t0
    print(f"\nlikelihood={like} reps={reps}  wall={dt:.1f}s "
          f"({dt/reps*1000:.0f} ms/eval)")


if __name__ == '__main__':
    main()
