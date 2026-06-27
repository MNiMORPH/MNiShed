[![DOI](https://zenodo.org/badge/199317220.svg)](https://doi.org/10.5281/zenodo.6787390)
[![Documentation Status](https://readthedocs.org/projects/mnished/badge/?version=latest)](https://mnished.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://img.shields.io/pypi/v/MNiShed)](https://pypi.org/project/MNiShed/)
[![Tests](https://github.com/MNiMORPH/MNiShed/actions/workflows/tests.yml/badge.svg)](https://github.com/MNiMORPH/MNiShed/actions/workflows/tests.yml)

# MNiShed
**Watershed hydrological model**

*MN* = Minnesota · *Mni* = water (Dakota) · *Mini* = small (lumped) · *Shed* = watershed

<!-- start-intro -->

MNiShed is a lumped, daily-timestep conceptual hydrological model. It routes precipitation through an optional snowpack stage and then through a cascade of one or more reservoirs (soil zone, groundwater, etc.), producing streamflow. A basin can be split into parallel sub-catchments — including open-water lakes — to represent contrasting hydrologic zones. Evapotranspiration is either read from a data file or computed with the Thornthwaite–Chang 2019 equation, optionally reshaped by a growing-degree-day vegetation-phenology coefficient, and scaled to close the long-run water balance. The model follows the [CSDMS Basic Model Interface (BMI)](https://csdms.colorado.edu/wiki/BMI), and ships with a fast in-process calibration stack and post-fit diagnostics.

<!-- end-intro -->

---

[Read the full documentation on ReadTheDocs](https://mnished.readthedocs.io/)

---

## Features

- **Process options** — snowpack (degree-day melt, rain-on-snow), frozen ground, a configurable reservoir cascade (power-law or linear recession), tile drainage (two formulations), PDM saturation-excess, and direct runoff.
- **Parallel sub-catchments & lakes** — split a basin into zones with independent physics; open-water lake elements with a threshold stage–discharge outlet, direct precipitation minus evaporation, and a bidirectional lake↔groundwater exchange.
- **Evapotranspiration** — measured or Thornthwaite–Chang, with an optional growing-degree-day vegetation-phenology coefficient and water-balance closure.
- **Calibration** — fast in-process calibration (build-once `ScoringModel`, config-driven `Calibrator` with declarative parameter targets), multi-window/decadal objectives, and SPOTPY (SCE-UA, DREAM) or Dakota drivers.
- **Analysis** — post-fit parameter-identifiability diagnostics (likelihood profiles, a curvature eigenspectrum, ridges), a seasonal mass-balance diagnostic, Brutsaert–Nieber recession analysis, and data-driven priors.
- **Performance & interoperability** — an optional Numba JIT time loop (verified identical to pure Python) and a CSDMS BMI wrapper.

---

## Installation

```bash
pip install mnished
```

To install from source for development:

```bash
git clone https://github.com/MNiMORPH/MNiShed.git
cd MNiShed
pip install -e '.[bmi,jit]'   # optional extras: bmi, jit, docs, lint
```

The optional extras are `bmi` (CSDMS BMI wrapper), `jit` (Numba JIT
acceleration), `docs` (build the documentation), and `lint` (Ruff).

## Quick start

**Python API**

```python
import mnished

b = mnished.Buckets()
b.initialize('config.yml')
b.run()
b.compute_NSE(verbose=True)
b.plot()
```

**Command-line interface**

```bash
mnished -y config.yml
```

See the [Quick Start guide](https://mnished.readthedocs.io/en/latest/quickstart.html) for the configuration file format and input data requirements.

## Citation

If you use MNiShed, please cite it using the metadata in [CITATION.cff](CITATION.cff) or via the Zenodo record:

> Wickert, A. D. (2026). MNiShed: Watershed hydrological model. https://doi.org/10.5281/zenodo.6787390

## Contact

Please report bugs and request features via [GitHub Issues](https://github.com/MNiMORPH/MNiShed/issues).
