#!/usr/bin/env python3
"""
Forward run of hydroRaVENS for the Cannon River catchment.

Runs the model with the configuration in cannon_cfg.yml, prints
goodness-of-fit metrics, and saves a diagnostic hydrograph plot.

Usage (from this directory):
    python run_forward.py
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from mnished import run_and_score

result = run_and_score(
    'cannon_cfg.yml',
    spin_up_cycles = 3,
    routing_N      = 2,
    routing_K      = 3.0,   # Nash-cascade storage time constant [days]
    metric         = 'KGE_logKGE_logFDC',
)

b = result.buckets
print(f'Score (KGE_logKGE_logFDC) = {result.score:.3f}')
print(f'KGE_logFDC                = {result.kge_logfdc:.3f}')
print(f'AIC                       = {result.aic:.1f}')
print(f'BFI  obs / mod            = {result.bfi_obs:.3f} / {result.bfi_mod:.3f}')

# --- Diagnostic plot ---
fig, (ax_p, ax_q) = plt.subplots(2, 1, figsize=(12, 6), sharex=True,
                                  gridspec_kw={'height_ratios': [1, 2.5],
                                               'hspace': 0.05})

dates = b.hydrodata['Date']

ax_p.bar(dates, b.hydrodata['Precipitation [mm/day]'],
         width=1, color='steelblue', alpha=0.7)
ax_p.set_ylabel('Precip.\n[mm/day]')
ax_p.invert_yaxis()
ax_p.yaxis.set_label_position('right')
ax_p.yaxis.tick_right()

ax_q.plot(dates, b.hydrodata['Specific Discharge [mm/day]'],
          color='royalblue', lw=1.5, label='Observed')
ax_q.plot(dates, b.hydrodata['Specific Discharge (modeled) [mm/day]'],
          color='k', lw=1.5, label='Modelled')
ax_q.set_ylabel('Specific discharge [mm/day]')
ax_q.set_xlabel('Date')
ax_q.set_ylim(bottom=0)
ax_q.legend()
ax_q.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
plt.setp(ax_q.get_xticklabels(), rotation=30, ha='right')

fig.suptitle('hydroRaVENS – Cannon River forward run')
plt.savefig('forward_run.png', dpi=150, bbox_inches='tight')
print('Plot saved to forward_run.png')
plt.show()
