API Reference
=============

This page documents the complete public API exported from the
``hydroravens`` package.  See also:

- :doc:`recession_analysis` — :class:`~hydroravens.BrutsaertNieber` (Brutsaert–Nieber recession analysis)
- :doc:`priors` — :class:`~hydroravens.Priors` and :func:`~hydroravens.suggest_priors` (data-driven priors)
- :doc:`calibration` — :func:`~hydroravens.run_and_score` and :class:`~hydroravens.CalibResult` (calibration workflow)

.. contents:: On this page
   :local:
   :depth: 1

Buckets
-------

The primary class for running a watershed simulation.

.. autoclass:: hydroravens.Buckets
   :members:
   :undoc-members:
   :show-inheritance:

Reservoir
---------

A single power-law (or linear) reservoir.

.. autoclass:: hydroravens.Reservoir
   :members:
   :undoc-members:
   :show-inheritance:

Snowpack
--------

Snowpack accumulation and positive-degree-day melt.

.. autoclass:: hydroravens.Snowpack
   :members:
   :undoc-members:
   :show-inheritance:

run_and_score
-------------

Run the model with a given parameter set and return goodness-of-fit metrics.

.. autofunction:: hydroravens.run_and_score

CalibResult
-----------

Named tuple returned by :func:`~hydroravens.run_and_score`.

.. autoclass:: hydroravens.CalibResult
   :no-inherited-members:

HydrographSeparation
--------------------

Spectral and time-domain decomposition of a discharge record into
reservoir-timescale components.

.. autoclass:: hydroravens.HydrographSeparation
   :members: fit, get_initial_conditions, get_parameter_priors, summary
   :member-order: bysource
   :show-inheritance:
