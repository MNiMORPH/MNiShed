#!/usr/bin/env python3
"""Write the observed log-flow vector for Dakota bayes_calibration.

Runs the model once and writes the scored observed log-flows to
``cannon_logq_obs.dat`` (one value per line, Dakota ``freeform``).  The
scoring mask (which days are scored) and the observed column depend only on
observation availability — not on the calibrated parameters — so this file is
fixed for the whole Bayesian run.  ``driver_bayes.py`` returns the matching
modelled log-flows as ``calibration_terms``; Dakota differences the two.

The model run here uses each parameter's ``initial`` value purely to populate
the modelled column (which is discarded); only the observed column is written.

Usage:
    python make_bayes_data.py
"""

import yaml
import numpy as np

from mnished import run_and_score
from mnished.calibration import log_flow_residual_terms

with open('params.yml') as f:
    cfg = yaml.safe_load(f)
P, drv, MODULES = cfg['parameters'], cfg['driver'], cfg.get('modules', {})


def _init(name):
    return P[name]['initial']


result = run_and_score(
    'cannon_cfg_template.yml',
    recession_coeff       = [10 ** _init('log__t_recession_shallow'),
                              10 ** _init('log__t_recession_soil'),
                              10 ** _init('log__t_recession_karst')],
    f_to_discharge        = [_init('f_exfiltration_shallow'),
                              _init('f_exfiltration_soil')],
    melt_factor           =  _init('PDD_melt_factor'),
    fdd_threshold         =  10 ** _init('log__fdd_threshold'),
    Hmax                  = [10 ** _init('log__Hmax_shallow')],
    modules               =  MODULES,
    routing_K             =  10 ** _init('log__routing_K'),
    routing_N             =  drv['routing_N'],
    start                 =  drv['decade_start'],
    end                   =  drv['decade_end'],
    spin_up_cycles        =  drv['spin_up_cycles'],
    metric                =  'KGE',
)

terms = log_flow_residual_terms(result, start=drv['decade_start'],
                                end=drv['decade_end'])
np.savetxt('cannon_logq_obs.dat', terms['log_obs'].to_numpy())
eps = 0.01 * terms['obs'].mean()
print(f"Wrote cannon_logq_obs.dat  ({len(terms)} scored days; eps = {eps:.6g})")
print("Now run: python generate_dakota_in.py --bayes")
