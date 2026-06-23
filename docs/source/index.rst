MNiShed Documentation
==========================

.. image:: https://zenodo.org/badge/199317220.svg
   :target: https://zenodo.org/badge/latestdoi/199317220
   :alt: Zenodo DOI

A simple, flexible reservoir-based hydrological model for water balance simulation and streamflow prediction.

**MNiShed** (*MN* = Minnesota, *Mni* = water in Dakota, *Mini* = small/lumped,
*Shed* = watershed) is a miniature lumped conceptual model that routes precipitation
through an optional snowpack and a cascade of reservoirs to produce streamflow.
Ideal for long water-balance studies, climate impact assessments, and ungauged basins.

**Key Features**
~~~~~~~~~~~~~~~~

* Optional snowpack module -- positive-degree-day melt with rain-on-snow sensible heat
* Optional frozen ground module -- Molnau & Bissell FGI blocks deep infiltration
* Cascading reservoirs -- linear or nonlinear (power-law) recession, stacked from soil to groundwater
* Flexible ET -- read from data or compute with the Thornthwaite--Chang equation
* Exact annual water balance -- ET scaled so P - Q - ET = 0 over each water year
* Calibration module -- KGE, NSE, and log-KGE scoring; AIC; baseflow index; flow
  duration curve; Nash-cascade channel routing; decadal chaining
* Python API and command-line interface
* Lightweight -- minimal dependencies (NumPy, Pandas, SciPy, Matplotlib, PyYAML)

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
   :caption: Integration

   bmi

.. toctree::
   :maxdepth: 2
   :caption: Performance

   benchmarks

Quick Example
~~~~~~~~~~~~~

**Python API:**

.. code-block:: python

    import mnished

    model = mnished.Buckets()
    model.initialize('config.yml')
    model.run()
    nse = model.compute_NSE(verbose=True)
    model.plot()

**Command-line:**

.. code-block:: bash

    mnished -y config.yml

Getting Help
~~~~~~~~~~~~

* **Report bugs:** `GitHub Issues <https://github.com/MNiMORPH/MNiShed/issues>`_
* **Discuss:** `GitHub Discussions <https://github.com/MNiMORPH/MNiShed/discussions>`_
* **Learn more:** `CSDMS Model Page <https://csdms.colorado.edu/wiki/Model:MNiShed>`_

About
~~~~~

MNiShed is developed and maintained by the `MNiMORPH group <https://github.com/MNiMORPH>`_.
It is published under the GNU General Public License v3.0.

**Citation:**

If you use MNiShed in your research, please cite it using the information in
``CITATION.cff`` at the repository root, or via the Zenodo DOI badge above.
