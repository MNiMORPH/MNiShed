#!/usr/bin/env python3
"""
Find the best-fit parameters from a completed Dakota calibration run,
re-run MNiShed with those parameters, and produce a diagnostic plot.

Figure layout
-------------
Left column  : precipitation (top, inverted) + observed/modelled discharge
Right column : flow duration curve (log scale) with observed BFI annotated

Usage (from cannon_inverse/):
    python plot_best.py                      # uses dakota.dat, saves best_fit.png
    python plot_best.py --dat dakota_test.dat --save test_fit.png
"""

import argparse
import sys
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from mnished import run_and_score
from mnished.calibration import _nse, _kge, _log_kge, _kge_logfdc

CFG_TEMPLATE  = 'cannon_cfg_template.yml'
OBJECTIVE_COL = 'neg_kge'
ROUTING_N     = 2      # Nash-cascade shape; must match driver.py ROUTING_N

_MODULE_PARAMS_MAP = {
    'snowpack':      ['PDD_melt_factor'],
    'frozen_ground': ['log__fdd_threshold', 'snow_insulation_k'],
    'direct_runoff': ['f_direct_runoff'],
    'rain_on_snow':  [],
}

def _load_params_yml(path):
    try:
        with open(path) as f:
            pcfg = yaml.safe_load(f)
        modules = pcfg.get('modules', {})
        params  = pcfg['parameters']
        # Mirror generate_dakota_in.py auto-fix so active flags match dakota.in.
        for mod, names in _MODULE_PARAMS_MAP.items():
            if not modules.get(mod, True):
                for name in names:
                    if name in params:
                        params[name]['active'] = False
        return (pcfg['driver']['metric'], modules, params)
    except FileNotFoundError:
        return 'KGE_logKGE_logFDC', {}, {}

# Deferred: populated after --params arg is parsed
METRIC  = None
MODULES = None
_PARAMS = None


def _is_active(name):
    return (_PARAMS or {}).get(name, {}).get('active', False)


def _get(row, name):
    """Return active param value from dat row, or fixed fallback from params.yml."""
    if _is_active(name):
        return float(row[name])
    return float((_PARAMS or {}).get(name, {}).get('fixed', 0.0))


def read_best_params(dat_file):
    try:
        df = pd.read_csv(dat_file, sep=r'\s+')
    except FileNotFoundError:
        sys.exit(f'Error: {dat_file} not found. Run Dakota first.')
    df = df.rename(columns={'%eval_id': 'eval_id'})
    for col in df.columns:
        if col != 'interface':
            df[col] = pd.to_numeric(df[col], errors='coerce')
    return df.loc[df[OBJECTIVE_COL].idxmin()]


def run_model(row):
    return run_and_score(
        CFG_TEMPLATE,
        recession_coeff       = [10 ** _get(row, 'log__t_recession_shallow'),
                                  10 ** _get(row, 'log__t_recession_soil'),
                                  10 ** _get(row, 'log__t_recession_karst')],
        f_to_discharge        = [_get(row, 'f_exfiltration_shallow'),
                                  _get(row, 'f_exfiltration_soil')],
        melt_factor           =  _get(row, 'PDD_melt_factor'),
        fdd_threshold         =  10 ** _get(row, 'log__fdd_threshold'),
        snow_insulation_k     =  _get(row, 'snow_insulation_k'),
        Hmax                  = [10 ** _get(row, 'log__Hmax_shallow')],
        direct_runoff_fraction=  _get(row, 'f_direct_runoff'),
        routing_K             =  10 ** _get(row, 'log__routing_K'),
        routing_N             =  ROUTING_N,
        modules               =  MODULES,
        metric                =  METRIC,
    )


def make_plot(result, params, save_path, metric=METRIC):
    b     = result.buckets
    score = result.score
    aic   = result.aic

    mask  = (b.hydrodata['Specific Discharge (modeled) [mm/day]'].notna()
             & b.hydrodata['Specific Discharge [mm/day]'].notna())
    m_all = np.asarray(b.hydrodata.loc[mask, 'Specific Discharge (modeled) [mm/day]'])
    o_all = np.asarray(b.hydrodata.loc[mask, 'Specific Discharge [mm/day]'])
    nse        = _nse(m_all, o_all)
    kge        = _kge(m_all, o_all)
    log_kge    = _log_kge(m_all, o_all)
    kge_logfdc = result.kge_logfdc

    dates = b.hydrodata['Date']

    # --- Figure layout: left column (time series) + right column (FDC) ---
    fig = plt.figure(figsize=(14, 7))
    gs  = fig.add_gridspec(2, 2, width_ratios=[3, 1], height_ratios=[1, 2.5],
                           hspace=0.05, wspace=0.25)
    ax_p   = fig.add_subplot(gs[0, 0])
    ax_q   = fig.add_subplot(gs[1, 0], sharex=ax_p)
    ax_fdc = fig.add_subplot(gs[:, 1])

    # --- Precipitation (inverted) ---
    ax_p.bar(dates, b.hydrodata['Precipitation [mm/day]'],
             width=1, color='steelblue', alpha=0.7)
    ax_p.set_ylabel('Precip.\n[mm/day]')
    ax_p.invert_yaxis()
    ax_p.yaxis.set_label_position('right')
    ax_p.yaxis.tick_right()
    plt.setp(ax_p.get_xticklabels(), visible=False)

    # --- Discharge time series ---
    ax_q.plot(dates, b.hydrodata['Specific Discharge [mm/day]'],
              color='royalblue', lw=1.5, label='Observed')
    ax_q.plot(dates, b.hydrodata['Specific Discharge (modeled) [mm/day]'],
              color='k', lw=1.5, label='Modelled')
    ax_q.set_ylabel('Specific discharge [mm/day]')
    ax_q.set_xlabel('Date')
    ax_q.set_ylim(bottom=0)
    ax_q.legend(loc='upper right', fontsize=9)
    ax_q.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    plt.setp(ax_q.get_xticklabels(), rotation=30, ha='right')

    # Annotation box
    t_shallow  = 10 ** _get(params, 'log__t_recession_shallow')
    t_soil     = 10 ** _get(params, 'log__t_recession_soil')
    t_karst    = 10 ** _get(params, 'log__t_recession_karst')
    routing_K  = 10 ** _get(params, 'log__routing_K')

    score_str = f'logKGE = {log_kge:.3f}   NSE = {nse:.3f}   KGE = {kge:.3f}   KGE$_{{logFDC}}$ = {kge_logfdc:.3f}   AIC = {aic:.1f}'
    param_lines = (
        f'BFI: obs = {result.bfi_obs:.3f},  mod = {result.bfi_mod:.3f}\n'
        f'$\\tau_{{sh}}$ = {t_shallow:.1f} d,  '
        f'$\\tau_{{soil}}$ = {t_soil:.0f} d,  '
        f'$\\tau_{{karst}}$ = {t_karst:.0f} d\n'
        f'$f_{{sh}}$ = {_get(params, "f_exfiltration_shallow"):.3f},  '
        f'$f_{{soil}}$ = {_get(params, "f_exfiltration_soil"):.3f}'
    )
    if _is_active('PDD_melt_factor'):
        param_lines += f',  PDD = {_get(params, "PDD_melt_factor"):.2f} mm °C$^{{-1}}$ d$^{{-1}}$'
    param_lines += f'\n$H_{{max}}$ = {10 ** _get(params, "log__Hmax_shallow"):.0f} mm'
    if _is_active('log__fdd_threshold'):
        param_lines += f',  FDD$_{{thresh}}$ = {10 ** _get(params, "log__fdd_threshold"):.0f} °C·d'
    if _is_active('snow_insulation_k'):
        param_lines += f',  $k_{{ins}}$ = {_get(params, "snow_insulation_k"):.4f} mm$^{{-1}}$ SWE'
    param_lines += f'\n$K_{{route}}$ = {routing_K:.2f} d  (N={ROUTING_N})'
    if _is_active('f_direct_runoff'):
        param_lines += f',  $\\gamma_{{direct}}$ = {_get(params, "f_direct_runoff"):.3f}'

    ann = score_str + '\n' + param_lines
    ax_q.text(0.02, 0.97, ann, transform=ax_q.transAxes,
              va='top', fontsize=8.5,
              bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))

    # --- Flow duration curve ---
    ax_fdc.semilogy(result.fdc_obs.index, result.fdc_obs.values,
                    color='royalblue', lw=1.5, label='Observed')
    ax_fdc.semilogy(result.fdc_mod.index, result.fdc_mod.values,
                    color='k', lw=1.5, label='Modelled')
    ax_fdc.set_xlabel('Exceedance probability [%]')
    ax_fdc.set_ylabel('Specific discharge [mm/day]')
    ax_fdc.set_xlim(0, 100)
    ax_fdc.legend(fontsize=9)
    ax_fdc.set_title('Flow duration curve', fontsize=10)
    ax_fdc.grid(True, which='both', alpha=0.3)

    fig.suptitle('MNiShed – Cannon River best-fit calibration', fontsize=13)

    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f'Figure saved to {save_path}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dat',    default='dakota.dat',   help='Dakota tabular data file')
    parser.add_argument('--save',   default='best_fit.png', help='Output figure path')
    parser.add_argument('--params', default='params.yml',   help='params.yml config file')
    parser.add_argument('--no-show', action='store_true',   help='Save only; skip plt.show()')
    args = parser.parse_args()

    METRIC, MODULES, _PARAMS = _load_params_yml(args.params)

    best = read_best_params(args.dat)

    print(f'\nBest evaluation: {int(best["eval_id"])}')
    print(f'  metric           = {METRIC}')
    print(f'  t_recession_shallow = {10 ** _get(best, "log__t_recession_shallow"):.1f} days')
    print(f'  t_recession_soil    = {10 ** _get(best, "log__t_recession_soil"):.0f} days')
    print(f'  t_recession_karst   = {10 ** _get(best, "log__t_recession_karst"):.0f} days')
    print(f'  f_exfilt_shallow = {_get(best, "f_exfiltration_shallow"):.4f}')
    print(f'  f_exfilt_soil    = {_get(best, "f_exfiltration_soil"):.4f}')
    if _is_active('PDD_melt_factor'):
        print(f'  PDD_melt_factor  = {_get(best, "PDD_melt_factor"):.4f} mm/°C/day')
    print(f'  Hmax_shallow     = {10 ** _get(best, "log__Hmax_shallow"):.1f} mm')
    if _is_active('log__fdd_threshold'):
        print(f'  fdd_threshold    = {10 ** _get(best, "log__fdd_threshold"):.1f} °C·day')
    if _is_active('snow_insulation_k'):
        print(f'  snow_insulation_k= {_get(best, "snow_insulation_k"):.4f} mm⁻¹ SWE')
    rK = 10 ** _get(best, 'log__routing_K')
    print(f'  routing_K        = {rK:.3f} days  (N={ROUTING_N},'
          f' mean travel time = {ROUTING_N * rK:.2f} days)')
    if _is_active('f_direct_runoff'):
        print(f'  f_direct_runoff  = {_get(best, "f_direct_runoff"):.4f}')

    result = run_model(best)
    b     = result.buckets
    mask  = (b.hydrodata['Specific Discharge (modeled) [mm/day]'].notna()
             & b.hydrodata['Specific Discharge [mm/day]'].notna())
    m_all = np.asarray(b.hydrodata.loc[mask, 'Specific Discharge (modeled) [mm/day]'])
    o_all = np.asarray(b.hydrodata.loc[mask, 'Specific Discharge [mm/day]'])
    print(f'  logKGE           = {_log_kge(m_all, o_all):.4f}')
    print(f'  NSE              = {_nse(m_all, o_all):.4f}')
    print(f'  KGE              = {_kge(m_all, o_all):.4f}')
    print(f'  KGE_logFDC       = {result.kge_logfdc:.4f}')
    print(f'  AIC              = {result.aic:.2f}')
    print(f'  BFI obs          = {result.bfi_obs:.4f}')
    print(f'  BFI mod          = {result.bfi_mod:.4f}')

    make_plot(result, best, save_path=args.save, metric=METRIC)
    if not args.no_show:
        plt.show()
