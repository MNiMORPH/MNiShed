#!/bin/bash
# Run a full calibration, plot the best fit, and archive results.
#
# Usage:
#   bash run.sh <short-description>
#   e.g.: bash run.sh kge_3res_nogamma
#
# Requires dakota and python (with mnished + dakota.interfacing) on PATH.
# If using a conda environment: conda activate <your-env> before running.

set -euo pipefail

DESC="${1:?Usage: bash run.sh <short-description>}"
TIMESTAMP=$(date +%Y-%m-%d_%H%M%S)
RUN_NAME="${TIMESTAMP}_${DESC}"

echo "=== Run: $RUN_NAME ==="

# Clean previous ephemeral outputs
rm -rf out dakota.dat dakota.out dakota.rst fort.13 LHS_*.out

# Optimise
dakota -i dakota.in -o dakota.out

# Save diagnostic figure
if python plot_best.py --save best_fit.png --no-show; then
    echo "Best-fit plot saved."
else
    echo "Warning: plot_best.py failed; archiving without plot." >&2
fi

# Archive
bash archive_run.sh "$RUN_NAME"

echo "=== Archived to runs/$RUN_NAME ==="
