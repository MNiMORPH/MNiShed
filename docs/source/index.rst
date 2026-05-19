HydroRaVENS Documentation
==========================

.. image:: https://zenodo.org/badge/199317220.svg
   :target: https://zenodo.org/badge/latestdoi/199317220
   :alt: Zenodo DOI

A simple, flexible reservoir-based hydrological model for water balance simulation and streamflow prediction.

**HydroRaVENS** (Rain and Variable Evapotranspiration, Nieve, and Streamflow) is a
conceptual daily-timestep model that routes precipitation through an optional snowpack
and a cascade of linear reservoirs. Ideal for ungauged basins, climate impact studies,
and educational applications.

**Key Features**
~~~~~~~~~~~~~~~~

* Optional snowpack module -- positive-degree-day melt with rain-on-snow sensible heat
* Optional frozen ground module -- Molnau & Bissell FGI blocks deep infiltration
* Cascading linear reservoirs -- stack multiple reservoirs from soil to groundwater
* Flexible ET -- read from data or compute with the Thornthwaite--Chang equation
* Exact annual water balance -- ET scaled so P - Q - ET = 0 over each water year
* Calibration module -- KGE, NSE, and log-KGE scoring; AIC; baseflow index; flow
  duration curve; Nash-cascade channel routing; decadal chaining
* Python API and command-line interface
* Lightweight -- minimal dependencies (NumPy, Pandas, Matplotlib, PyYAML)

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   quickstart
   tutorial

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   configuration
   model_description

.. toctree::
   :maxdepth: 2
   :caption: Reference

   api
   calibration
   recession_analysis
   priors
   references

.. toctree::
   :maxdepth: 2
   :caption: CSDMS Basic Model Interface

   bmi

Quick Example
~~~~~~~~~~~~~

**Python API:**

.. code-block:: python

    import hydroravens

    model = hydroravens.Buckets()
    model.initialize('config.yml')
    model.run()
    nse = model.compute_NSE(verbose=True)
    model.plot()

**Command-line:**

.. code-block:: bash

    hydroravens -y config.yml

Getting Help
~~~~~~~~~~~~

* **Report bugs:** `GitHub Issues <https://github.com/MNiMORPH/hydroRaVENS/issues>`_
* **Discuss:** `GitHub Discussions <https://github.com/MNiMORPH/hydroRaVENS/discussions>`_
* **Learn more:** `CSDMS Model Page <https://csdms.colorado.edu/wiki/Model:HydroRaVENS>`_

About
~~~~~

HydroRaVENS is developed and maintained by the `MNiMORPH group <https://github.com/MNiMORPH>`_.
It is published under the GNU General Public License v3.0.

**Citation:**

If you use HydroRaVENS in your research, please cite it using the information in
``CITATION.cff`` at the repository root, or via the Zenodo DOI badge above.
