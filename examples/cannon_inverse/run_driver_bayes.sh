#!/bin/bash
# Dakota calls this script for each Bayesian (DREAM) function evaluation.
# Requires python and dakota.interfacing to be on PATH / importable.
# If using a conda environment: conda activate <your-env> before running Dakota.
python driver_bayes.py "$@"
