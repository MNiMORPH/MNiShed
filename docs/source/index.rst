MNiShed Documentation
==========================

.. image:: https://zenodo.org/badge/199317220.svg
   :target: https://zenodo.org/badge/latestdoi/199317220
   :alt: Zenodo DOI

A simple, flexible reservoir-based hydrological model for water balance simulation and streamflow prediction.

**MNiShed** (*MN* = Minnesota, *Mni* = water in Dakota, *Mini* = small/lumped,
*Shed* = watershed) is a lumped conceptual model that routes precipitation
through an optional snowpack and one or more parallel sub-catchments — each a
cascade of reservoirs — to produce streamflow.
Ideal for long water-balance studies, climate impact assessments, and ungauged basins.

**Key Features**
~~~~~~~~~~~~~~~~

* Cascading reservoirs -- linear or nonlinear (power-law) recession, soil to
  groundwater, with optional leakance/threshold junctions and tile drainage
* Parallel sub-catchments -- partition a basin into spatially distinct
  hydraulic zones (e.g. till uplands vs. lake-clay lowlands), each its own
  reservoir cascade, area-weighted to the outlet
* Calibration -- KGE, NSE, and log-KGE scoring; AIC model comparison; baseflow
  index; flow duration curve; Nash-cascade routing; decadal chaining
* Fast -- Numba JIT-compiled daily time loop, roughly 100--400× faster than
  pure Python (see :doc:`benchmarks`)
* Data-driven setup -- estimate recession exponents (Brutsaert--Nieber) and
  parameter priors directly from observed streamflow
* Exact annual water balance -- ET scaled so P - Q - ET = 0 over each water year
* Flexible ET -- from data or the Thornthwaite--Chang equation; optional
  storage-dependent and reservoir-draw ET
* Optional snowpack module -- positive-degree-day melt with rain-on-snow sensible heat
* Optional frozen ground module -- Molnau & Bissell FGI blocks deep infiltration
* CSDMS Basic Model Interface -- couple MNiShed with other models
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
   example_crow_wing

.. toctree::
   :maxdepth: 2
   :caption: Reference

   api
   input_contract
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
