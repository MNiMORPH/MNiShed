#!/usr/bin/env python3
"""Dakota DREAM (Bayesian) driver for MNiShed log-flow calibration.

Where ``driver.py`` returns the scalar ``1 - KGE`` for an optimiser, this
driver returns the **modelled log-flow vector** as Dakota ``calibration_terms``.
Dakota differences it against the observed log-flows in the
``calibration_data_file`` and forms the Gaussian likelihood, with
``calibrate_error_multipliers`` inferring the error scale.  Parameters and run
settings are read from ``params.yml`` (shared with the optimisation driver).

Generate the data file and the Bayesian Dakota input first:

    python make_bayes_data.py
    python generate_dakota_in.py --bayes
    dakota -i dakota_bayes.in -o dakota_bayes.out

Note: DREAM needs many thousands of model evaluations, so the Numba JIT should
be active (the pure-Python fallback is ~100x slower).
"""

import yaml
import numpy as np

from mnished import run_and_score
from mnished.calibration import log_flow_residual_terms

with open('params.yml') as f:
    _cfg = yaml.safe_load(f)
_driver    = _cfg['driver']
_param_cfg = _cfg['parameters']
SPIN_UP_CYCLES = _driver['spin_up_cycles']
ROUTING_N      = _driver['routing_N']
DECADE_START   = _driver['decade_start']
DECADE_END     = _driver['decade_end']
MODULES        = _cfg.get('modules', {})

# Mirror generate_dakota_in.py's module auto-fix so active flags match the
# Dakota variables block.
_MODULE_PARAMS = {
    'snowpack':      ['PDD_melt_factor'],
    'frozen_ground': ['log__fdd_threshold', 'snow_insulation_k'],
    'direct_runoff': ['f_direct_runoff'],
    'rain_on_snow':  [],
}
for _mod, _names in _MODULE_PARAMS.items():
    if not MODULES.get(_mod, True):
        for _name in _names:
            if _name in _param_cfg:
                _param_cfg[_name]['active'] = False


def modeled_log_flows(get):
    """Modelled log-flow vector over the scored window.

    ``get(name)`` returns the active Dakota value or the fixed fallback,
    exactly as in ``driver.py``.  Returned vector aligns 1:1 with the
    observed log-flows written by ``make_bayes_data.py`` (the scoring mask
    is parameter-independent), so it is the Dakota ``calibration_terms``.
    """
    result = run_and_score(
        'cannon_cfg_template.yml',
        recession_coeff       = [10 ** get('log__t_recession_shallow'),
                                  10 ** get('log__t_recession_soil'),
                                  10 ** get('log__t_recession_karst')],
        f_to_discharge        = [get('f_exfiltration_shallow'),
                                  get('f_exfiltration_soil')],
        melt_factor           =  get('PDD_melt_factor'),
        fdd_threshold         =  10 ** get('log__fdd_threshold'),
        Hmax                  = [10 ** get('log__Hmax_shallow')],
        direct_runoff_fraction=  get('f_direct_runoff'),
        modules               =  MODULES,
        routing_K             =  10 ** get('log__routing_K'),
        routing_N             =  ROUTING_N,
        start                 =  DECADE_START,
        end                   =  DECADE_END,
        spin_up_cycles        =  SPIN_UP_CYCLES,
        metric                =  'KGE',
    )
    terms = log_flow_residual_terms(result, start=DECADE_START, end=DECADE_END)
    return terms['log_mod'].to_numpy()


def main():
    import dakota.interfacing as di
    params, results = di.read_parameters_file()

    def get(name):
        p = _param_cfg[name]
        return params[name] if p['active'] else p['fixed']

    labels = list(results)
    try:
        log_mod = modeled_log_flows(get)
        ok = bool(np.all(np.isfinite(log_mod))) and len(log_mod) == len(labels)
    except Exception:
        ok = False

    if ok:
        for lab, val in zip(labels, log_mod):
            results[lab].function = float(val)
    else:
        # A far-off prediction → negligible likelihood, so the chain rejects
        # this sample rather than crashing the run.
        for lab in labels:
            results[lab].function = -10.0
    results.write()


if __name__ == '__main__':
    main()
