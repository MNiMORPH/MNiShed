[![DOI](https://zenodo.org/badge/199317220.svg)](https://doi.org/10.5281/zenodo.6787390)
[![Documentation Status](https://readthedocs.org/projects/mnished/badge/?version=latest)](https://mnished.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://img.shields.io/pypi/v/MNiShed)](https://pypi.org/project/MNiShed/)
[![Tests](https://github.com/MNiMORPH/MNiShed/actions/workflows/tests.yml/badge.svg)](https://github.com/MNiMORPH/MNiShed/actions/workflows/tests.yml)

# :black_nib: MNiShed
**Rain and Variable Evapotranspiration, Nieve, and Streamflow**

<!-- start-intro -->

MNiShed is a lumped, daily-timestep conceptual hydrological model. It routes precipitation through an optional snowpack stage and then through a cascade of one or more linear reservoirs (soil zone, groundwater, etc.), producing streamflow. Evapotranspiration is either read from a data file or computed with the Thornthwaite–Chang 2019 equation, and is scaled per water year so that the long-run water balance closes. The model follows the [CSDMS Basic Model Interface (BMI)](https://csdms.colorado.edu/wiki/BMI).

<!-- end-intro -->

---

[Read the full documentation on ReadTheDocs](https://mnished.readthedocs.io/)

---

## Installation

```bash
pip install mnished
```

To install from source for development:

```bash
git clone https://github.com/MNiMORPH/MNiShed.git
cd MNiShed
pip install -e .
```

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

> Wickert, A. D. (2026). MNiShed: Rain and Variable Evapotranspiration, Nieve, and Streamflow. https://doi.org/10.5281/zenodo.6787390

## Contact

Please report bugs and request features via [GitHub Issues](https://github.com/MNiMORPH/MNiShed/issues).
