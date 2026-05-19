[![DOI](https://zenodo.org/badge/199317220.svg)](https://doi.org/10.5281/zenodo.6787390)
[![Documentation Status](https://readthedocs.org/projects/hydroravens/badge/?version=latest)](https://hydroravens.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://img.shields.io/pypi/v/hydroRaVENS)](https://pypi.org/project/hydroRaVENS/)
[![Tests](https://github.com/MNiMORPH/hydroRaVENS/actions/workflows/tests.yml/badge.svg)](https://github.com/MNiMORPH/hydroRaVENS/actions/workflows/tests.yml)

# :black_nib: HydroRaVENS
**Rain and Variable Evapotranspiration, Nieve, and Streamflow**

<!-- start-intro -->

HydroRaVENS is a lumped, daily-timestep conceptual hydrological model. It routes precipitation through an optional snowpack stage and then through a cascade of one or more linear reservoirs (soil zone, groundwater, etc.), producing streamflow. Evapotranspiration is either read from a data file or computed with the Thornthwaite–Chang 2019 equation, and is scaled per water year so that the long-run water balance closes. The model follows the [CSDMS Basic Model Interface (BMI)](https://csdms.colorado.edu/wiki/BMI).

<!-- end-intro -->

---

[Read the full documentation on ReadTheDocs](https://hydroravens.readthedocs.io/)

---

## Installation

```bash
pip install hydroravens
```

To install from source for development:

```bash
git clone https://github.com/MNiMORPH/hydroRaVENS.git
cd hydroRaVENS
pip install -e .
```

## Quick start

**Python API**

```python
import hydroravens

b = hydroravens.Buckets()
b.initialize('config.yml')
b.run()
b.compute_NSE(verbose=True)
b.plot()
```

**Command-line interface**

```bash
hydroravens -y config.yml
```

See the [Quick Start guide](https://hydroravens.readthedocs.io/en/latest/quickstart.html) for the configuration file format and input data requirements.

## Citation

If you use HydroRaVENS, please cite it using the metadata in [CITATION.cff](CITATION.cff) or via the Zenodo record:

> Wickert, A. D. (2026). HydroRaVENS: Rain and Variable Evapotranspiration, Nieve, and Streamflow. https://doi.org/10.5281/zenodo.6787390

## Contact

Please report bugs and request features via [GitHub Issues](https://github.com/MNiMORPH/hydroRaVENS/issues).
