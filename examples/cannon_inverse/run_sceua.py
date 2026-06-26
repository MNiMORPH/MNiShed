#!/usr/bin/env python3
"""In-process best-fit calibration via SPOTPY SCE-UA — the optimizer analogue
of run_spotpy.py's DREAM. Same setup pattern (ParameterSet + run_and_score),
maximises KGE in-process (warm JIT) instead of fanning out through Dakota's
fork interface. For an A/B against the Dakota EGO + pattern-search workflow.

    python run_sceua.py [repetitions]
"""
import sys
import time

import numpy as np
import yaml
import spotpy

from mnished import ParameterSet, ScoringModel

with open('params.yml') as f:
    _cfg = yaml.safe_load(f)
P, DRV, MOD = _cfg['parameters'], _cfg['driver'], _cfg.get('modules', {})
PSET = ParameterSet.from_params_yml(P)
NAMES = PSET.names


def _g(n, th):
    return th.get(n, P[n]['fixed'])


class CannonOpt:
    def __init__(self):
        self.params = [spotpy.parameter.Uniform(p.name, p.lower, p.upper,
                                                optguess=p.value)
                       for p in PSET]
        self.evals = 0
        # Build the model once; every evaluation reuses it (no per-eval CSV
        # re-read or reconstruction). Bit-identical to run_and_score.
        self.model = ScoringModel('cannon_cfg_template.yml')

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        self.evals += 1
        th = dict(zip(NAMES, vector))
        try:
            s = self.model.score(
                recession_coeff=[10 ** _g('log__t_recession_shallow', th),
                                 10 ** _g('log__t_recession_soil', th),
                                 10 ** _g('log__t_recession_karst', th)],
                f_to_discharge=[_g('f_exfiltration_shallow', th),
                                _g('f_exfiltration_soil', th)],
                melt_factor=_g('PDD_melt_factor', th),
                fdd_threshold=10 ** _g('log__fdd_threshold', th),
                Hmax=[10 ** _g('log__Hmax_shallow', th)],
                modules=MOD, routing_K=10 ** _g('log__routing_K', th),
                routing_N=DRV['routing_N'],
                spin_up_cycles=DRV['spin_up_cycles'], metric='KGE',
                start=DRV['decade_start'], end=DRV['decade_end']).score
            return [s if np.isfinite(s) else -10.0]
        except Exception:
            return [-10.0]

    def evaluation(self):
        return [1.0]                       # ideal KGE

    def objectivefunction(self, simulation, evaluation):
        return 1.0 - simulation[0]         # minimise 1-KGE -> maximise KGE


def main():
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
    setup = CannonOpt()
    t0 = time.time()
    sampler = spotpy.algorithms.sceua(setup, dbname='cannon_sceua',
                                      dbformat='ram')
    sampler.sample(reps)
    dt = time.time() - t0
    print(f"\nSCE-UA in-process: {setup.evals} evals, wall={dt:.1f}s "
          f"({dt / max(setup.evals, 1) * 1000:.0f} ms/eval)")
    print(f"best KGE = {1.0 - sampler.status.objectivefunction_min:.4f}")


if __name__ == '__main__':
    main()
