Recession Analysis
==================

MNiShed includes tools for analysing the nonlinear storage–discharge
relationship directly from observed streamflow records, before running or
calibrating the model.

.. contents:: On this page
   :local:
   :depth: 2

Theory
------

During a baseflow recession (no rainfall or snowmelt), water drains from
storage at a rate set by the storage–discharge relationship.  If discharge
is a power law of storage depth,

.. math::

    Q = c \, H^b,

then substituting :math:`\mathrm{d}H/\mathrm{d}t = -Q` gives

.. math::

    -\frac{\mathrm{d}Q}{\mathrm{d}t} = K \, Q^n,
    \qquad n = \frac{2b - 1}{b},

where :math:`b` is the MNiShed ``recession_exponent`` and :math:`n`
is the slope on a log–log plot of :math:`-\mathrm{d}Q/\mathrm{d}t`
versus :math:`Q` — the Brutsaert & Nieber (1977) recession plot.

Inverting:

.. math::

    b = \frac{1}{2 - n}

The two exponents have the following special values:

.. list-table::
   :widths: 20 20 50
   :header-rows: 1

   * - B–N slope *n*
     - MNiShed *b*
     - Physical interpretation
   * - 1
     - 1
     - Linear reservoir; exponential recession.
   * - 3/2
     - 2
     - Long-time Boussinesq solution for horizontal, unconfined aquifer
       (Brutsaert & Nieber 1977).
   * - → 2
     - → ∞
     - Extreme nonlinearity; storage empties nearly instantaneously.

The conversion :math:`b = 1/(2-n)` is valid only for :math:`n < 2`.
The Boussinesq short-time solution gives :math:`n = 3`, which is outside
this framework; in practice, early fast-response recessions (surface
runoff, tile drains) produce steep tails on the B–N cloud that should not
be used to set the baseflow recession exponent.

Brutsaert–Nieber Analysis
--------------------------

.. autoclass:: mnished.BrutsaertNieber
   :members: fit, to_reservoir_exponent, summary, plot
   :member-order: bysource

Using Results in MNiShed
-----------------------------

The B–N analysis produces a single catchment-integrated exponent *b*.
In a multi-reservoir MNiShed model, the reservoirs have different
physical roles and the exponent should be assigned accordingly:

.. list-table::
   :widths: 20 20 60
   :header-rows: 1

   * - Reservoir
     - Suggested *b*
     - Rationale
   * - Soil (fastest)
     - B–N estimate
     - The full recession cloud is dominated by the soil/fast response.
       Use as a calibration starting point; typical calibrated values
       are *b* ≈ 3–4 for agricultural catchments.
   * - Intermediate
     - B–N lower-envelope estimate (fix, do not calibrate)
     - The lower envelope of the B–N cloud represents the *slowest
       observable* recession pathway — typically the shallow Quaternary
       zone (outwash, fractured regolith; MRT days to weeks).  The truly
       deep reservoir drains too slowly and at too low a flux to register
       in the cloud, so B–N is capturing the intermediate, not the deep.
       The long-time Boussinesq value (:math:`b = 2.0`) is a reasonable
       lower bound; basin-specific fits typically yield :math:`b \approx
       2.0`–2.5. Fix rather than calibrate (see equifinality note in
       :doc:`model_description`).
   * - Deep (slowest)
     - 1.0 (fixed)
     - Mean residence time of decades to centuries; contributes too
       little flux to constrain *b* from streamflow.  The linear
       approximation (Darcy flow in a confined or semi-confined aquifer)
       is physically appropriate and avoids adding a free parameter.

The B–N slope *n* from the **lower envelope** of the cloud corresponds
to long-duration recessions dominated by the intermediate subsurface
zone and is the most physically meaningful region for setting the
intermediate-reservoir exponent.
The upper scatter reflects short, event-driven recessions (tile drains,
surface runoff) and should not be used to set the intermediate exponent.

For a complete prior-estimation workflow that combines the B–N analysis
with timescale estimation, see :func:`~mnished.suggest_priors` and
the :doc:`tutorial`.

Workflow
--------

A typical workflow fits the recession curve and uses the result as a
prior for model calibration:

.. code-block:: python

    import pandas as pd
    from mnished import BrutsaertNieber

    df = pd.read_csv('streamflow.csv', parse_dates=['Date'])
    Q  = df['Specific Discharge [mm/day]'].values

    bn = BrutsaertNieber(Q, min_recession_days=3).fit()
    bn.summary()
    bn.plot()

    b_prior = bn.to_reservoir_exponent()
    print(f"Suggested recession_exponent prior: {b_prior:.2f}")

The output might look like::

    Brutsaert–Nieber recession analysis
      Recession pairs used : 412
      Fitted slope  n      : 1.5831
      Fitted coeff  a      : 0.0147
      MNiShed   b      : 2.367
      Reference (long-time Boussinesq): n = 1.5, b = 2.0

The lower envelope of the B–N cloud corresponds to long-duration
recessions dominated by the intermediate subsurface zone (shallow
Quaternary units, fractured regolith; MRT days to weeks) and is the most
physically meaningful region for setting the intermediate-reservoir
exponent.  The truly deep reservoir (MRT decades–centuries) produces
discharge too small to dominate any part of the cloud and should be
fixed at :math:`b = 1` independently.  The upper scatter reflects
short, event-driven recessions (tile drains, surface runoff) that are
better captured by a calibrated soil-zone exponent.

Goodness of fit and the power-law assumption
--------------------------------------------

:meth:`~mnished.BrutsaertNieber.fit` automatically assesses whether
the recession cloud is well described by a power law and issues a
``UserWarning`` when it is not.

Two diagnostics are computed and stored as attributes:

* **R² of the power-law fit** (:attr:`~mnished.BrutsaertNieber.r2_`)
  — coefficient of determination of the linear log–log regression.
  The B–N cloud is inherently noisy (all types of recession contribute),
  so R² values of 0.4–0.7 are typical for real catchments.  R² < 0.4
  indicates a poor overall fit; inspect the plot.

* **Curvature test** (:attr:`~mnished.BrutsaertNieber.r2_quadratic_`)
  — R² of a degree-2 polynomial fit in log–log space.  If the quadratic
  fit is substantially better than the linear fit (ΔR² > 0.05), the
  recession cloud is systematically curved.  A straight line in log–log
  space is the signature of a power-law storage–discharge relationship;
  systematic curvature means this assumption does not hold.

**What to do when the power-law assumption fails:**

A curved recession cloud in log–log space signals that the catchment
sensitivity function g(Q) — as defined by Kirchner (2009) — is not
a power law.  This can arise from:

* Multiple reservoirs draining simultaneously, each with different
  exponents, producing a composite curve.  Consider whether a
  multi-reservoir model with different fixed exponents captures the
  cloud's shape better than a single-exponent fit.
* A genuine threshold or non-power-law drainage process (e.g. a
  perched water table that drains rapidly when full and slowly when
  nearly empty).
* Data quality issues (precipitation misattribution, rating curve
  errors at low flow).

In the first case, the lower envelope still provides a useful prior for
the intermediate reservoir exponent.  In the second case, a different
reservoir functional form may be needed — one that MNiShed does not
currently implement, following Kirchner's (2009) more general framework.
In the third case, inspect and correct the input data before re-fitting.

Interpretation guide
--------------------

* **n ≈ 1.5, b ≈ 2** — recession consistent with the Boussinesq
  long-time solution for a horizontal unconfined aquifer.  A reasonable
  fixed prior for the **intermediate** reservoir (shallow Quaternary
  units, fractured regolith).  The deep reservoir should be fixed at
  :math:`b = 1` regardless of the B–N result.
* **n ≈ 1.3–1.6, b ≈ 1.5–2.5** — typical mixed catchment range.
* **n > 1.7, b > 3** — rapid nonlinear drainage, often reflecting
  near-surface soil processes.  Calibrate rather than fix.
* **n ≥ 2** — unphysical in this framework; refit using only the lower
  envelope or with a longer ``min_recession_days``.

References
----------

Brutsaert, W. and Nieber, J. L. (1977). Regionalized drought flow
hydrographs from a mature glaciated plateau. *Water Resources Research*,
13(3), 637–643. https://doi.org/10.1029/WR013i003p00637

Kirchner, J. W. (2009). Catchments as simple dynamical systems:
Catchment characterization, rainfall-runoff modeling, and doing hydrology
backward. *Water Resources Research*, 45, W02429.
https://doi.org/10.1029/2008WR006912
