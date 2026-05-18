Model Description
==================

Theory & Mathematical Formulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

HydroRaVENS is a lumped, daily-timestep conceptual hydrological model. Water 
moves through three sequential processes each day: optional snowpack accumulation/melt, 
routing through cascading linear reservoirs, and evapotranspiration.

All fluxes are expressed as depths over the drainage basin (mm/day).
Mass is conserved to within numerical precision.

Daily Water Balance
~~~~~~~~~~~~~~~~~~~

On each day, the model computes:

.. math::

    P + M - E = \Delta S + Q

where:

* :math:`P` = precipitation (mm/day)
* :math:`M` = snowmelt (mm/day)
* :math:`E` = evapotranspiration (mm/day)
* :math:`\Delta S` = change in storage (mm/day)
* :math:`Q` = streamflow (mm/day)

Optional Process Modules
~~~~~~~~~~~~~~~~~~~~~~~~~

Several processes can be individually enabled or disabled through the
``modules`` block in the configuration file (see :doc:`configuration`).
All modules that require temperature data are silently inactive when
``Mean Temperature [C]`` is absent from the input CSV, regardless of
their flag setting.

.. list-table::
   :widths: 20 10 60
   :header-rows: 1

   * - Module
     - Default
     - Purpose
   * - ``snowpack``
     - on
     - Accumulation, PDD melt, sublimation, and rain-on-snow.
   * - ``frozen_ground``
     - on
     - Frozen ground index (FGI) that blocks deep infiltration.
       Requires ``snowpack: true``.
   * - ``rain_on_snow``
     - on
     - Sensible-heat contribution of warm rain on snowpack.
       Requires ``snowpack: true``.
   * - ``direct_runoff``
     - off
     - Hortonian-inspired bypass fraction. See below.
   * - ``dtr_fgi_decay``
     - on
     - DTR-based FGI decay. When T_min/T_max columns are present, the
       per-day decay coefficient varies with the diurnal temperature
       range. Disable to revert to constant ``fgi_decay_coeff``
       (original Molnau & Bissell behaviour). See Frozen Ground Module.

Snowpack Module (Optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~

If mean air temperature is provided, snowpack processes are enabled.

**Accumulation:**
  When :math:`T \leq 0°C`, net water input (precipitation minus ET) accumulates
  as snow (stored as SWE):

  .. math::

      \text{SWE}_{t+1} = \text{SWE}_t + (P_t - E_t)

**Melt:**
  When :math:`T > 0°C`, melt is computed using the positive-degree-day (PDD) approach:

  .. math::

      M_t = \min(\text{SWE}_t,\ \alpha \cdot T_t \cdot \Delta t + M_{\text{ROS},t})

  where :math:`\alpha` is the melt factor (mm SWE °C⁻¹ day⁻¹) and
  :math:`M_{\text{ROS}}` is the rain-on-snow sensible-heat contribution
  (see below). All melt is routed directly to the top reservoir.

**Rain-on-snow (ROS) sensible heat:**
  Rain arriving at temperature :math:`T > 0°C` carries thermal energy that
  can melt additional snow:

  .. math::

      M_{\text{ROS},t} = \frac{c_p}{L_f} \cdot T_t \cdot P_t

  where :math:`c_p / L_f \approx 0.01253\ °C^{-1}` is the ratio of the
  specific heat of water to the latent heat of fusion. This term is
  typically small relative to the PDD term except during warm rain-on-snow
  events (McCabe et al. 2007; Würzer et al. 2016).

**ET deficit:**
  When precipitation minus ET is negative, the deficit first sublimates
  snow:

  .. math::

      \text{Sublimation} = \min\!\left(\text{SWE}_t,\ \max(0,\ E_t - P_t)\right)

  Any deficit exceeding available SWE is passed to the top subsurface
  reservoir and, if still unmet, carried forward to the next time step.

Frozen Ground Module (Optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When a ``fdd_threshold`` is set, the model tracks a **frozen ground index**
(FGI; Molnau & Bissell 1983) that accumulates freezing degree-days,
decays passively each day, and additionally decays during warm periods:

.. math::

    \text{FGI}_t = \max\!\left(0,\ A \cdot \text{FGI}_{t-1} - T_{\text{eff},t} - D_t\right)

where :math:`A_t` is the daily decay coefficient (see below),
:math:`T_{\text{eff},t} = T_t \cdot e^{-k \cdot \text{SWE}_t}` is the
snow-insulation-adjusted temperature (°C; negative values increase FGI,
positive values reduce it), and :math:`D_t` is an additional thaw credit
described below. When :math:`\text{FGI}_t` exceeds ``fdd_threshold``,
the top reservoir's exfiltration fraction is set to 1.0 so that all
drainage becomes direct runoff, simulating frozen-soil blockage of deep
infiltration.

**DTR-based decay coefficient** :math:`A_t`:
  The decay coefficient :math:`A_t` represents sub-daily heat input to
  frozen soil that is not resolved by the daily-mean temperature forcing.
  Its dominant physical driver is the fraction of the day when air
  temperature is above 0°C even though the daily mean is negative —
  a process that is frequent in maritime climates (Pacific Northwest,
  Atlantic Europe) and rare in continental ones (central North America,
  Siberia).

  When daily minimum and maximum temperature columns
  (``Minimum Temperature [C]``, ``Maximum Temperature [C]``) are present
  in the input CSV, :math:`A_t` is computed from the diurnal temperature
  range (DTR):

  .. math::

      f_{\text{above}} = \frac{\max(0,\ T_{\max,t})}
                              {T_{\max,t} - T_{\min,t}}

      A_t = 1 - (1 - A_0)\,f_{\text{above}}

  where :math:`A_0` is ``fgi_decay_coeff`` (default 0.97, M&B 1983) and
  :math:`f_{\text{above}}` is the fraction of the day above 0°C under a
  linear diurnal-cycle assumption. On days entirely below freezing
  (:math:`T_{\max} \leq 0`), :math:`A_t = 1.0` — no passive decay;
  FGI accumulates unimpeded. When the diurnal cycle straddles 0°C,
  :math:`A_t` falls toward :math:`A_0`. This naturally produces
  near-unity :math:`A` for continental winters (where :math:`T_{\max}`
  rarely crosses 0°C during cold spells) and M&B-like decay for maritime
  climates with frequent overnight freeze–daytime thaw cycles.

  ``fgi_decay_coeff`` thus represents the *maximum* daily decay rate,
  reached when the entire cold day oscillates across 0°C. The steady-state
  :math:`\text{FGI}^* = |T| / (1 - A_t)` is now climate-dependent:
  large (persistent frost) for continental conditions, bounded for
  maritime ones.

  If T_min/T_max are absent, :math:`A_t` falls back to the constant
  ``fgi_decay_coeff`` (original M&B behaviour).

**Coupling snowmelt and frozen-ground thaw via the melt factor:**
  The PDD melt factor :math:`\alpha` has units of mm SWE per °C·day,
  making it a natural conversion factor between the thermal forcing
  (°C·day) and a water-equivalent depth (mm SWE):

  .. math::

      \text{mm SWE} = \alpha \cdot \text{°C·day}
      \qquad \Longleftrightarrow \qquad
      \text{°C·day} = \frac{\text{mm SWE}}{\alpha}

  When total melt energy in a timestep exceeds the available SWE, the
  leftover energy :math:`\Delta E` is expressed in mm SWE. Dividing by
  :math:`\alpha` converts it back to degree-days — the same currency the
  FGI uses — so the residual can be credited toward thawing frozen ground:

  .. math::

      D_t = \frac{\max(0,\ \alpha T_t \Delta t + M_{\text{ROS},t} - \text{SWE}_t)}{\alpha}

  This means that once the snowpack is fully depleted, any remaining
  thermal energy continues to thaw frozen soil rather than being
  discarded. The melt factor thus serves as the bridge between the two
  empirical degree-day representations.

  Note that :math:`\alpha` characterises the snow surface, not the soil.
  Applying the same conversion to the soil implicitly assumes that the
  atmosphere-to-surface thermal coupling is the same in both cases, which
  is a simplification. Within a degree-day framework, however, the
  approach is internally self-consistent.

**Snow insulation:**
  A deep snowpack buffers the soil surface from cold air, reducing the
  effective freezing degree-days that accumulate in the FGI. This is
  represented by an exponential insulation factor:

  .. math::

      T_{\text{eff},t} = T_t \cdot e^{-k \cdot \text{SWE}_t}

  where :math:`k` (mm⁻¹ SWE) is the snow insulation decay constant
  (``snow_insulation_k`` in the ``snowmelt`` config section; default
  0.0). The factor is applied to both freezing and thawing temperature
  forcing; meltwater heat delivery (``excess_dd``) is not scaled because
  meltwater reaches the soil surface directly. The parameterisation
  originates in Molnau & Bissell (1983), and was adopted by LISFLOOD
  (van der Knijff et al. 2010) and GSSHA (Downer & Ogden 2004).

  .. note::

      ``snow_insulation_k`` and ``fdd_threshold`` are correlated: both
      control how much frozen-ground effect the model sees. Calibrating
      them simultaneously from streamflow alone leads to equifinality —
      the optimizer trades a near-zero threshold against moderate
      insulation rather than finding physically meaningful values for
      either. Recommended practice:

      * Fix ``snow_insulation_k`` at a literature or field-derived value
        and calibrate only ``fdd_threshold``, or
      * Leave ``snow_insulation_k = 0`` (default) and treat
        ``fdd_threshold`` as the sole free FGI parameter.

      The insulation term is most useful when independent observations
      (soil temperature, frost depth) are available to constrain
      :math:`k`, or for deep-snowpack alpine catchments where the
      insulation effect is large relative to the threshold uncertainty.

Regional Groundwater Import (Optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When ``baseflow_Q > 0`` in the ``catchment`` configuration section, a
constant daily flux (mm/day) is added to modeled discharge after all
reservoir routing:

.. math::

    Q_{\text{out},t} = Q_{\text{routed},t} + Q_{\text{base}}

This represents regional groundwater inflow from outside the surface
catchment — for example, deep confined-aquifer discharge or inter-basin
transfer that is decoupled from the local shallow water balance.
:math:`Q_{\text{base}}` is not mass-balanced against precipitation or
ET; it adds water to the stream without a corresponding source in the
reservoir cascade.

Use with care: ``baseflow_Q`` is most physically justified when
independent hydrogeological evidence supports an external groundwater
source (artesian springs, regional flow systems). Calibrating it from
streamflow alone risks compensating for other structural deficiencies
in the model.

Direct Runoff Bypass (Optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When ``direct_runoff: true`` in the ``modules`` block, a fixed fraction
:math:`\gamma` of positive daily recharge is intercepted before it enters
the top reservoir and routed directly to the stream:

.. math::

    Q_{\text{direct},t} = \gamma \cdot \max(0,\ R_t)

    R_{\text{remaining},t} = (1 - \gamma) \cdot R_t

where :math:`R_t` is the recharge available after ET (mm/day) and
:math:`\gamma` is ``direct_runoff_fraction`` in the ``general``
configuration section.

This formulation is *Hortonian-inspired* — conceptually motivated by
infiltration-excess overland flow — but cannot be a rigorous
physical representation at the daily timestep. Without sub-daily
intensity data, it is impossible to determine whether the daily
precipitation total actually exceeds hydraulic conductivity at any
moment. The bypass may nonetheless provide a useful empirical
correction for catchments with significant impervious area, compacted
soils, or urban fractions, and for days dominated by a single intense
storm event (Yilmaz et al. 2008).

The module is off by default (``direct_runoff: false``). Unless the
calibration objective clearly demands it, leaving it off is recommended
to avoid equifinality: the bypass can mimic effects already represented
by the shallow-reservoir timescale, causing parameters to lose their
intended physical interpretation.

Linear Reservoir Cascade
~~~~~~~~~~~~~~~~~~~~~~~~

Water drains through a stack of reservoirs (top = shallowest, bottom = deepest).
Each reservoir first receives its recharge input, then drains by exponential decay:

.. math::

    Q_i(t) = \bigl(H_i(t) + Q_{\text{recharge},i}\bigr) \cdot (1 - e^{-\Delta t / \tau_i})

where:

* :math:`H_i` = water depth in reservoir :math:`i` at the start of the time step (mm)
* :math:`Q_{\text{recharge},i}` = water input to reservoir :math:`i` this time step (mm)
* :math:`\tau_i` = e-folding residence time (days)
* :math:`\Delta t` = time step (1 day)

**Discharge Partitioning:**
  Of the water drained from reservoir :math:`i`, a fraction :math:`f_i` exits as 
  river discharge; the remainder infiltrates to the next layer:
  
  .. math::
  
      Q_{\text{discharge},i} = f_i \cdot Q_i \\
      Q_{\text{infiltrate},i} = (1 - f_i) \cdot Q_i

**Constraint:**
  The bottom reservoir should fully discharge (:math:`f_{\text{bottom}} = 1.0`).
  A warning is issued if not, as this violates mass conservation.

**Storage Update:**
  Recharge is applied first, then exponential drainage:

  .. math::

      H_i(t+1) = \bigl(H_i(t) + Q_{\text{recharge},i}\bigr)\,e^{-\Delta t/\tau_i}

  Overflow above :math:`H_{\max}` exits immediately as direct runoff; any
  deficit is passed to the next-deeper reservoir. ET is not subtracted
  separately at the reservoir level — it is already incorporated into the
  recharge input to the top reservoir.

**Physical interpretation:**
  No reservoir is fixed to a particular process; meaning is set by
  parameter choice. Successive reservoirs naturally span progressively
  longer timescales — interflow (days), soil moisture (months),
  groundwater (years) — but that mapping is the user's choice, analogous
  to the multi-component runoff structure of HBV (Bergström 1976).

Nonlinear (Power-Law) Recession (Optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each reservoir can optionally use a power-law storage–discharge relationship
in place of the default linear (exponential) decay. The discharge rate scales
as a power of storage:

.. math::

    Q_i = \frac{H_i}{\tau_i} \left(\frac{H_i}{H_{\text{ref}}}\right)^{b_i - 1}

where :math:`b_i \geq 1` is the recession exponent for reservoir :math:`i`
and :math:`H_{\text{ref}}` is a reference storage depth at which :math:`\tau_i`
retains its standard linear meaning. Setting :math:`b_i = 1` recovers the
linear reservoir exactly.

For :math:`b > 1`, the relationship is superlinear: high-storage states drain
faster than the linear equivalent, and as storage depletes the drainage rate
slows more rapidly. This behaviour is observed in natural catchments and has
theoretical grounding in subsurface flow geometry (Brutsaert & Nieber 1977;
Kirchner 2009). Brutsaert & Nieber (1977) derived :math:`b \approx 2.2` from
analysis of recession hydrographs; individual catchments and reservoirs may
deviate substantially from this value.

.. note::

    In calibration against streamflow, the recession exponent of the deepest
    calibrated reservoir (typically the karst or bedrock zone) tends to be
    poorly constrained and may drift toward high values that mimic
    threshold-like behavior (see below). Fixing deep-reservoir exponents at
    the theoretical Brutsaert-Nieber value (:math:`b = 2.2`) is recommended
    unless independent evidence supports a different value.

Mean Residence Time for Nonlinear Reservoirs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Because :math:`\tau` is defined as the drainage timescale specifically at
the 1 mm reference storage :math:`H_{\text{ref}}`, it is **not directly
comparable** across reservoirs with different recession exponents. Two
reservoirs with identical :math:`\tau` but different :math:`b` drain at
very different rates at any realistic storage depth above 1 mm, because
the power-law nonlinearity amplifies (for :math:`b > 1`) or suppresses the
drainage rate relative to the linear case.

A physically meaningful alternative is the **mean residence time** (MRT)
at a representative steady-state discharge :math:`Q_{\text{ref}}`. At
steady state the outflow equals the input, so
:math:`Q_{\text{ref}} = H_{ss}^b / \tau`, giving:

.. math::

    H_{ss} = \bigl(Q_{\text{ref}}\cdot\tau\bigr)^{1/b}

    \text{MRT} = \frac{H_{ss}}{Q_{\text{ref}}} =
        \frac{\tau^{1/b}}{Q_{\text{ref}}^{\,1-1/b}}

For :math:`b = 1` (linear) this reduces exactly to :math:`\tau`,
consistent with the standard result for exponential decay. For
:math:`b > 1`, MRT is smaller than :math:`\tau` whenever
:math:`Q_{\text{ref}} > 1` mm/day, reflecting the accelerated drainage at
operating storage depths well above :math:`H_{\text{ref}}`.

MRT is implemented as
:meth:`~hydroravens.hydroravens.Reservoir.mean_residence_time`::

    mrt = reservoir.mean_residence_time(Q_ref=0.5)  # [mm/day] → [days]

**Choosing** :math:`Q_{\text{ref}}`: use the long-term mean discharge
attributed to the reservoir. For a cascade this is approximately the
mean annual flux passing through that layer. The result is most useful
as a cross-experiment or cross-basin diagnostic; its absolute value
depends on the chosen reference flux.

.. note::

    Calibration studies (Wickert 2026) show that raw :math:`\tau` values
    for strongly nonlinear reservoirs (:math:`b \approx 4`–5, typical of
    tile-drain-influenced soil zones) can differ by more than an order of
    magnitude between experiments while MRT remains nearly constant at
    ~5–6 days, because the optimizer adjusts :math:`\tau` to compensate
    for changes in :math:`b` and the operating storage level. MRT is
    therefore the recommended diagnostic for comparing calibrated
    timescales across experiments or basins, and should be reported
    alongside :math:`\tau` and :math:`b` whenever the recession exponent
    differs among reservoirs.

Saturation-Excess Runoff (PDM, Optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The probability-distributed model (PDM; Moore 1985) represents spatial
variability in soil storage capacity within a reservoir. Storage capacities
are assumed to follow an exponential distribution with mean :math:`H_0`
(``pdm_H0``). When the reservoir storage :math:`H` is positive, the
fraction of the area already saturated is:

.. math::

    f_{\text{sat}} = 1 - e^{-H / H_0}

That fraction of each day's positive recharge runs off immediately as
saturation-excess overland flow; the remainder enters the reservoir normally.
Setting a very large :math:`H_0` renders the PDM effectively inactive.

.. note::

    **PDM and power-law recession are equifinal.** Both produce nonlinear
    storage–discharge behaviour, but through different mechanisms: the PDM
    uses a sharp storage threshold (saturation-excess), while power-law
    recession is a smooth, continuous relationship. When both are active and
    calibrated simultaneously against streamflow alone, the optimizer can trade
    one against the other without penalty, obscuring the physical meaning of
    each parameter.

    Calibration experience on the Cannon River (Wickert 2026) shows that when
    the recession exponent :math:`b` is free, the PDM characteristic depth
    :math:`H_0` is driven to large values (effectively inactive), and vice
    versa: when PDM is removed, :math:`b` rises to compensate. The two
    mechanisms are redundant for the purpose of matching a streamflow record.

    **Recommended practice:** use one or the other, not both.

    * Power-law recession (:math:`b > 1`) is preferred when the goal is a
      smooth, theoretically grounded nonlinearity. It has a direct connection
      to subsurface flow theory (Brutsaert & Nieber 1977; Kirchner 2009) and
      avoids the threshold artifact.
    * PDM is preferred when independent evidence (saturation mapping, soil
      moisture observations) supports a genuine threshold process.

**On shallow reservoirs, tile drains, and macropore flow:**
  The same equifinality argument extends to structural choices about reservoir
  count and fast-pathway modules. A shallow reservoir with a short timescale,
  a tile-drain bypass, or a macropore fraction all represent fast, near-surface
  runoff generation. If the soil-zone recession exponent is calibrated and found
  to be robustly high (:math:`b \approx 3`–:math:`4`), the soil reservoir is
  already capturing this fast-response behavior through its nonlinear
  storage–discharge relationship. Adding a separate shallow reservoir, tile, or
  macropore pathway then produces similar discharge dynamics with additional
  parameters, potentially leading to false conclusions because the parameter
  names imply distinct physical mechanisms that the data cannot distinguish.

  The same equifinality extends to the frozen-ground module. Frozen-ground
  blockage of infiltration during winter and early spring produces elevated
  fast discharge — the same behavioral signature as a nonlinear soil recession
  with high :math:`b`. When both are active and calibrated simultaneously,
  the optimizer may drive the frozen-ground parameters (``fdd_threshold``,
  ``snow_insulation_k``) toward values that effectively disable the module,
  with :math:`b_{\text{soil}}` compensating. Streamflow alone cannot
  distinguish frozen-soil bypassing from high-storage nonlinear drainage.

  Unless there is independent process evidence (tile-drain monitoring, tracer
  data, direct macropore observations, or soil temperature / frost-depth records)
  that requires a mechanistically separate fast or frozen pathway, prefer
  calibrating :math:`b_{\text{soil}}` and accepting the result as an effective
  representation of the aggregate fast-response behavior of the soil zone.

Evapotranspiration
~~~~~~~~~~~~~~~~~~

Two methods are supported:

**Method 1: From Data**
  ET is read directly from the input CSV and scaled to close the annual water balance:
  
  .. math::
  
      E_{\text{scaled}} = E_{\text{observed}} \cdot \frac{P_{\text{annual}} - Q_{\text{observed,annual}}}{ET_{\text{observed,annual}}}

**Method 2: Thornthwaite-Chang (2019)**
  Daily reference ET (:math:`ET_0`) is estimated from temperature and photoperiod:
  
  .. math::
  
      ET_0 = \text{f}(T_{\max}, T_{\min}, \text{photoperiod})
  
  Then scaled by water year to match observed P − Q balance.

In both cases, the annual scaling factor is stored and applied to ensure that 
:math:`P - Q - E = 0` over each water year.

Model Skill Evaluation
~~~~~~~~~~~~~~~~~~~~~~

Several goodness-of-fit metrics are available via the ``metric`` argument
of :func:`~hydroravens.calibration.run_and_score`.

**Kling-Gupta Efficiency (KGE)** (Gupta et al. 2009):

.. math::

    \text{KGE} = 1 - \sqrt{(r-1)^2 + (\alpha-1)^2 + (\beta-1)^2}

where :math:`r` is the Pearson correlation, :math:`\alpha = \sigma_m/\sigma_o`
is the variability ratio, and :math:`\beta = \mu_m/\mu_o` is the bias ratio.
KGE = 1 is perfect; KGE > 0.5 is generally considered satisfactory.

**Log-space KGE (logKGE):** KGE computed on :math:`\log(Q + \epsilon)`,
emphasising low-flow performance. :math:`\epsilon` is set to 1% of the
mean observed discharge.

**Nash-Sutcliffe Efficiency (NSE)**:

.. math::

    \text{NSE} = 1 - \frac{\sum_t (Q_{\text{mod},t} - Q_{\text{obs},t})^2}
                          {\sum_t (Q_{\text{obs},t} - \bar{Q}_{\text{obs}})^2}

NSE = 1 is perfect; NSE = 0 means the mean is as good as the model.

**Available composite metrics:**

.. list-table::
   :widths: 35 65
   :header-rows: 1

   * - ``metric`` key
     - Formula (equal weights)
   * - ``KGE``
     - :math:`\text{KGE}(Q, Q_{\text{obs}})`
   * - ``NSE``
     - :math:`\text{NSE}(Q, Q_{\text{obs}})`
   * - ``logKGE``
     - :math:`\text{KGE}(\log Q, \log Q_{\text{obs}})`
   * - ``KGE_logKGE``
     - :math:`\tfrac{1}{2}(\text{KGE} + \text{logKGE})`
   * - ``KGE_logKGE_logFDC``
     - :math:`\tfrac{1}{3}(\text{KGE} + \text{logKGE} + \text{KGE}_{\log\text{FDC}})`
   * - ``KGE_logKGE_logFDC_BFI``
     - :math:`\tfrac{1}{4}(\text{KGE} + \text{logKGE} + \text{KGE}_{\log\text{FDC}} + \text{BFI\_score})`
   * - ``logKGE_logFDC_BFI``
     - :math:`\tfrac{1}{3}(\text{logKGE} + \text{KGE}_{\log\text{FDC}} + \text{BFI\_score})`

:math:`\text{KGE}_{\log\text{FDC}}` is KGE computed on the log-transformed
flow-duration curve (sorted discharge). :math:`\text{BFI\_score} = 1 - |\text{BFI}_{\text{mod}}/\text{BFI}_{\text{obs}} - 1|`
penalises baseflow bias using the Eckhardt (2005) digital filter
(:math:`\alpha = 0.98`, :math:`\text{BFI}_{\max} = 0.80` for perennial streams).
All composite scores are bounded above by 1.0.

Model Assumptions & Limitations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Strengths:**

* ✅ Mass-conserving (exact balance)
* ✅ Minimal data requirements (daily P and Q)
* ✅ Fast computation (runs decades in seconds)
* ✅ Physically interpretable parameters
* ✅ Suitable for ungauged basins (transfer parameters)

**Limitations:**

* ⚠️ Lumped model (no spatial variability)
* ⚠️ Lumped nonlinearity — power-law recession and PDM both approximate
  threshold behavior but cannot distinguish the underlying mechanism
  (macropores, tile drains, saturation-excess) from streamflow alone
* ⚠️ Daily timestep by design — degree-day snowmelt, Thornthwaite ET, and
  linear reservoir drainage are daily-scale parameterisations that lose
  physical meaning at finer resolution
* ⚠️ Simplified groundwater–surface-water interaction
* ⚠️ Snowpack simplified (PDD model; ignores energy balance)
* ⚠️ No representation of lakes or artificial storage
* ⚠️ ET method choice matters; Thornthwaite-Chang is approximate

**Best suited for:**

* Long-term water balance studies (years–decades)
* Climate impact assessments
* Ungauged or poorly-gauged basins
* Parameter transfer to similar basins
* Educational modeling
