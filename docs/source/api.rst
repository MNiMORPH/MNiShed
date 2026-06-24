API Reference
=============

This page documents the complete public API exported from the
``mnished`` package.  See also:

- :doc:`recession_analysis` — :class:`~mnished.BrutsaertNieber` (Brutsaert–Nieber recession analysis)
- :doc:`priors` — :class:`~mnished.Priors` and :func:`~mnished.suggest_priors` (data-driven priors)
- :doc:`calibration` — :func:`~mnished.run_and_score` and :class:`~mnished.CalibResult` (calibration workflow)
- :doc:`bmi` — :class:`~mnished.BmiMNiShed` (CSDMS Basic Model Interface wrapper)

.. contents:: On this page
   :local:
   :depth: 1

Buckets
-------

The primary class for running a watershed simulation.

.. autoclass:: mnished.Buckets
   :members:
   :undoc-members:
   :show-inheritance:

Reservoir
---------

A single power-law (or linear) reservoir.

.. autoclass:: mnished.Reservoir
   :members:
   :undoc-members:
   :show-inheritance:

Snowpack
--------

Snowpack accumulation and positive-degree-day melt.

.. autoclass:: mnished.Snowpack
   :members:
   :undoc-members:
   :show-inheritance:

SubCatchment
------------

A parallel hydraulic compartment of the basin with its own reservoir cascade
and per-compartment snowpack and frozen-ground state.

.. autoclass:: mnished.SubCatchment
   :members:
   :undoc-members:
   :show-inheritance:

run_and_score
-------------

Run the model with a given parameter set and return goodness-of-fit metrics.

.. autofunction:: mnished.run_and_score

CalibResult
-----------

Named tuple returned by :func:`~mnished.run_and_score`.

.. autoclass:: mnished.CalibResult
   :no-inherited-members:

HydrographSeparation
--------------------

Spectral and time-domain decomposition of a discharge record into
reservoir-timescale components.

.. autoclass:: mnished.HydrographSeparation
   :members: fit, get_initial_conditions, get_parameter_priors, summary
   :member-order: bysource
   :show-inheritance:
