Data-Driven Priors
==================

:func:`~mnished.suggest_priors` combines
:class:`~mnished.BrutsaertNieber` recession analysis with
:class:`~mnished.HydrographSeparation` to produce a coherent set of
parameter starting points before any model run or calibration.

See the :doc:`tutorial` for a full worked example.

.. autofunction:: mnished.suggest_priors

.. autoclass:: mnished.Priors
   :members: summary, to_yaml_snippet
   :member-order: bysource
