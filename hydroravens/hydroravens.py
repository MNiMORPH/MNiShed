#! /usr/bin/python3

########################################
# Then, methought, the air grew denser #
#         - Edgar Allan Poe            #
#              THE RAVEN               #
########################################

# Started by A. Wickert
# 25 July 2019
# Updated slightly by J. Jones
# 08 Oct 2019
# Significant Update by A. Wickert
# October 2022
# CLI added by A. Wickert
# November 2023

import argparse
import math
import sys
import warnings

import matplotlib.dates as mdates
import numpy as np
import pandas as pd
import yaml
from matplotlib import pyplot as plt

# c_p / L_f: water's specific heat divided by the latent heat of fusion.
# c_p = 4186 J kg⁻¹ °C⁻¹, L_f = 334 000 J kg⁻¹  →  ≈ 0.01253 °C⁻¹
# Gives mm SWE melted per mm of rain per °C of rain temperature.
_CP_LF = 4186.0 / 334000.0


class Reservoir(object):
    """
    Generic reservoir. Accepts new water (recharge), and sends it to other
    reservoirs and/or out of the system (discharge) at a rate that is
    proportional to the amount of water held in the reservoir.
    """

    def __init__(self, t_recession, f_to_discharge=1., Hmax=np.inf, pdm_H0=None,
                 H0=0., f_tile=0.0, tau_tile=None):
        """
        Initialize a reservoir.

        Parameters
        ----------
        t_recession : float
            Recession time scale [days]. For a linear reservoir (recession
            exponent b = 1) this equals the true e-folding time. For b > 1
            the actual mean residence time depends on storage level; use
            :meth:`mean_residence_time` to obtain a physically comparable
            timescale at a given reference discharge.
        f_to_discharge : float, optional
            Fraction of water lost each time step that exits as river
            discharge. The remainder (1 - f_to_discharge) infiltrates to
            the next-deeper reservoir. Default 1.0 (all to discharge).
        Hmax : float, optional
            Maximum water depth the reservoir can hold. Default np.inf.
        pdm_H0 : float or None, optional
            Characteristic storage depth [mm] for the probability-distributed
            model (PDM) of saturation-excess overland flow.  Storage capacity
            is assumed exponentially distributed across the catchment with
            mean pdm_H0; the saturated fraction when reservoir depth is H is
            f_sat = 1 - exp(-H / pdm_H0).  That fraction of each positive
            recharge pulse is shed immediately as overland flow.  Mutually
            exclusive with a finite Hmax.  Default None (PDM off).
        H0 : float, optional
            Initial water depth at the start of the simulation. Default 0.
        f_tile : float, optional
            Fraction of subsurface drainage (H_infiltrated) that is diverted
            to a fast tile-drain sub-reservoir instead of passing to the
            next-deeper reservoir.  The tile sub-reservoir drains directly to
            stream with e-folding time tau_tile.  Represents agricultural
            tile drainage or any fast lateral subsurface path.  Default 0
            (no tile drainage).  Requires tau_tile when > 0.
        tau_tile : float or None, optional
            E-folding residence time [days] of the tile-drain sub-reservoir.
            Typical values: 3–21 days for agricultural tile systems.  Required
            when f_tile > 0; ignored when f_tile == 0.  Default None.

        Raises
        ------
        ValueError
            If t_recession <= 0, f_to_discharge < 0 or > 1, Hmax < 0,
            pdm_H0 <= 0, f_tile < 0 or > 1, or f_tile > 0 with no tau_tile.
        """
        self.Hwater = H0
        self.Hmax = Hmax
        self.pdm_H0 = pdm_H0
        self.t_recession = t_recession
        self.f_to_discharge = f_to_discharge

        # Initialized here so all instance attributes exist before
        # recharge() and discharge() are first called
        self.H_excess = 0.
        self.H_deficit = 0.
        self.H_exfiltrated = 0.
        self.H_infiltrated = 0.
        self.H_discharge = 0.

        # Check values and note whether they are reasonable
        if t_recession <= 0:
            raise ValueError("t_recession must be > 0.")
        if f_to_discharge < 0:
            raise ValueError("Negative f_to_discharge not possible.")
        elif f_to_discharge > 1:
            raise ValueError("f_to_discharge: Cannot discharge >100% of water.")
        elif f_to_discharge == 0:
            warnings.warn("All water infiltrates when f_to_discharge is 0:"+
                          " you may have created a\n"+
                          "redundant pass-through water-storage layer")
        if Hmax < 0:
            raise ValueError("Hmax must be >= 0 (and >0 makes more sense)")
        if pdm_H0 is not None and pdm_H0 <= 0:
            raise ValueError("pdm_H0 must be > 0")
        if f_tile < 0 or f_tile > 1:
            raise ValueError("f_tile must be in [0, 1]")
        if f_tile > 0 and tau_tile is None:
            raise ValueError("tau_tile must be provided when f_tile > 0")

        self.f_tile = f_tile
        if f_tile > 0.0 and tau_tile is not None:
            self.tile_res = Reservoir(tau_tile, f_to_discharge=1.0)
        else:
            self.tile_res = None

        # Power-law recession: Q = (H/τ) * (H/recession_H_ref)^(recession_exponent-1).
        # recession_exponent=1 (default) recovers the linear reservoir exactly.
        # recession_H_ref is the storage at which τ has its usual linear meaning [mm].
        self.recession_exponent = 1.0
        self.recession_H_ref    = 1.0

    def recharge(self, H):
        """
        Add or remove water from the reservoir.

        Recharge H can be positive (net water input, e.g. P > ET) or
        negative (net deficit, e.g. ET > P). Sets H_excess if the
        reservoir overflows Hmax, or H_deficit if more water is removed
        than the reservoir holds.

        Parameters
        ----------
        H : float
            Depth of water added (positive) or removed (negative).

        Raises
        ------
        ValueError
            If Hwater is already negative before recharge is applied.
        """
        # Extra water above a maximum cap
        self.H_excess = 0.
        # Water that this layer cannot hold and cannot be passed to a deeper layer
        self.H_deficit = 0.

        # ERROR if water is less than 0 -- may be able to remove
        # this check later
        if self.Hwater < 0:
            raise ValueError("Hwater in reservoir < 0; non-physical")

        # PDM: exponential distribution of storage capacities.
        # Saturated fraction of catchment sheds positive recharge immediately
        # as saturation-excess overland flow before the remainder enters storage.
        if self.pdm_H0 is not None and H > 0:
            f_sat = 1.0 - np.exp(-self.Hwater / self.pdm_H0)
            self.H_excess = f_sat * H
            H = (1.0 - f_sat) * H

        # What if more water is lost during "recharge" than exists in reservoir?
        # Create a deficit and bring Hwater to 0
        if self.Hwater + H < 0:
            self.H_deficit += self.Hwater + H
            self.Hwater = 0.
        # What if more water is added than maximum reservoir capacity?
        # Mark excess (straight to runoff) and bring Hwater to Hmax
        elif self.Hwater + H > self.Hmax:
            self.H_excess += self.Hwater + H - self.Hmax
            self.Hwater = self.Hmax
        # Otherwise, we're in a range in which 0 <= H <= Hmax
        # Yay! Things are easier!
        else:
            self.Hwater += H

    def discharge(self, dt):
        """
        Discharge water from the reservoir over one time step.

        Computes water lost by exponential decay, partitions it between
        river discharge (H_exfiltrated) and infiltration to the next-deeper
        reservoir (H_infiltrated), and adds overflow from recharge()
        (H_excess) to H_discharge.

        Parameters
        ----------
        dt : float
            Time step duration (same units as t_recession; typically days).
        """
        b   = self.recession_exponent
        H0  = self.Hwater
        if b == 1.0 or H0 <= 0.0:
            dH = H0 * (1 - np.exp(-dt / self.t_recession))
        else:
            # Exact integration of dH/dt = -(H/t_recession)·(H/H_ref)^(b-1)
            #   = -H^b / (t_recession · H_ref^(b-1))
            # Substituting u = H^(1-b):  du/dt = (b-1)/tau_eff
            # => H(t+dt) = [H0^(1-b) + (b-1)·dt/tau_eff]^(1/(1-b))
            tau_eff = self.t_recession * self.recession_H_ref ** (b - 1.0)
            H_new   = (H0 ** (1.0 - b) + (b - 1.0) * dt / tau_eff) ** (1.0 / (1.0 - b))
            dH      = H0 - max(0.0, H_new)
        self.H_exfiltrated = dH * self.f_to_discharge
        self.H_discharge = self.H_excess + self.H_exfiltrated
        self.H_infiltrated = dH * (1 - self.f_to_discharge)
        self.Hwater -= dH

        # Tile drainage: divert f_tile of H_infiltrated to the fast sub-reservoir.
        # The remainder continues to the next-deeper reservoir as normal.
        if self.tile_res is not None:
            tile_in = self.f_tile * self.H_infiltrated
            self.H_infiltrated -= tile_in
            self.tile_res.recharge(tile_in)
            self.tile_res.discharge(dt)
            self.H_discharge += self.tile_res.H_discharge

    def mean_residence_time(self, Q_ref):
        """
        Mean residence time [days] at a reference steady-state discharge.

        For a nonlinear reservoir governed by
        :math:`Q = (H/\\tau)\\cdot(H/H_{\\mathrm{ref}})^{b-1}`, the
        steady-state storage at discharge *Q_ref* is
        :math:`H_{ss} = (Q_{\\mathrm{ref}} \\cdot \\tau)^{1/b}`, giving

        .. math::

            \\mathrm{MRT} = \\frac{H_{ss}}{Q_{\\mathrm{ref}}} =
                \\frac{\\tau^{1/b}}{Q_{\\mathrm{ref}}^{\\,1 - 1/b}}

        For a linear reservoir (*b* = 1) this reduces exactly to
        :math:`\\tau`.  For *b* > 1, MRT is smaller than :math:`\\tau`
        whenever :math:`Q_{\\mathrm{ref}} > 1` mm/day, reflecting the
        faster drainage at realistic operating storage depths.

        Unlike :attr:`t_recession`, which is the recession time scale only at the
        1 mm reference storage, MRT is a physically comparable timescale
        across reservoirs with different recession exponents.

        Parameters
        ----------
        Q_ref : float
            Representative steady-state discharge from this reservoir
            [mm/day]. Use the long-term mean flux attributed to this
            layer (e.g. mean annual discharge partitioned by reservoir).

        Returns
        -------
        float
            Mean residence time [days].

        Raises
        ------
        ValueError
            If Q_ref <= 0.
        """
        if Q_ref <= 0:
            raise ValueError("Q_ref must be > 0.")
        b = self.recession_exponent
        if b == 1.0:
            return self.t_recession
        return self.t_recession ** (1.0 / b) / Q_ref ** (1.0 - 1.0 / b)


class Snowpack(object):
    """
    Snowpack reservoir driven by temperature.

    Accumulates precipitation as snow when mean temperature is at or below
    0 °C. Melts at a positive-degree-day rate when temperature is above 0 °C.
    All melt is routed to the top subsurface reservoir as infiltration; direct
    discharge to the river is not modeled.

    Should precede the subsurface reservoir list in a watershed model.
    The melt factor is a positive-degree-day factor [mm/°C/day].
    """

    def __init__(self, melt_factor=None):
        """
        Initialize an empty snowpack.

        Parameters
        ----------
        melt_factor : float, optional
            Positive-degree-day melt factor (mm SWE °C⁻¹ day⁻¹).
            Can be set or updated later via set_melt_factor().
        """
        self.Hwater = 0.  # SWE
        self.melt_factor = melt_factor
        self.T = 0.
        self.H_infiltrated = 0.
        self.H_deficit = 0.

    def set_melt_factor(self, melt_factor):
        """
        Set or update the positive-degree-day melt factor.

        Parameters
        ----------
        melt_factor : float
            Melt rate per positive degree-day (mm SWE °C⁻¹ day⁻¹).
        """
        self.melt_factor = melt_factor

    def set_temperature(self, T):
        """
        Set the mean air temperature for the current time step.

        Parameters
        ----------
        T : float
            Mean air temperature (°C).
        """
        self.T = T

    def recharge(self, H):
        """
        Apply net water input or deficit to the snowpack.

        If T <= 0, positive H accumulates as snow (SWE). If T > 0,
        positive H bypasses the snowpack and is passed directly to the
        top subsurface reservoir via H_infiltrated. Negative H (ET > P)
        is removed from the snowpack as sublimation; any remainder that
        exceeds available SWE becomes H_deficit.

        Parameters
        ----------
        H : float
            Net water depth for this time step (mm). Positive = input
            (P - ET > 0); negative = deficit (ET - P > 0).
        """

        self.H_deficit = 0.  # Water deficit with neg ET; just this time step
        # If positive recharge
        if H >= 0:
            if self.T <= 0:
                self.Hwater += H
                self.H_infiltrated = 0.
            else:
                # Incoming precip component; melt sums with this
                # This is then directly passed to the first layer of the
                # set of hydrological reservoirs
                self.H_infiltrated = H
        # If negative recharge: remove water from snowpack via sublimation.
        # Any deficit beyond available SWE is passed down as H_deficit.
        else:
            # Sublimation (effectively) if snow present;
            # Otherwise pass water deficit
            if self.Hwater > -H:
                self.Hwater += H
            else:
                self.H_deficit += H + self.Hwater
                self.Hwater = 0
            self.H_infiltrated = 0.

    def melt(self, dt, P=0.0):
        """
        Compute positive-degree-day and rain-on-snow melt; update state.

        Both terms are routed to H_infiltrated (→ top soil reservoir). If
        total available energy exceeds the SWE present, the leftover is
        returned as equivalent degree-days so the caller can credit it
        toward frozen-soil thawing (FGI reduction) rather than losing it.

        Rain-on-snow sensible heat: water arriving at T_mean > 0 °C
        carries (c_p / L_f) · T · P mm SWE of latent-heat capacity.
        Spring snowpacks are near-isothermal at 0 °C, so cold-content
        corrections are negligible and the latent-heat term dominates.

        References
        ----------
        McCabe et al. (2007) doi:10.1175/BAMS-88-3-319
        Würzer et al. (2016) doi:10.1175/JHM-D-15-0181.1

        Parameters
        ----------
        dt : float
            Timestep [days].
        P : float, optional
            Raw liquid precipitation [mm/day]. Used to compute rain-on-snow
            sensible-heat melt. Default 0 (PDD only).

        Returns
        -------
        excess_dd : float
            Leftover melt energy after the snowpack is exhausted, expressed
            as degree-day equivalent [°C·day] = leftover mm SWE / melt_factor.
            Zero when SWE is not fully depleted.

            The melt factor (mm SWE °C⁻¹ day⁻¹) serves as the bridge
            between the PDD snowmelt representation and the frozen ground
            index (FGI): dividing excess melt depth (mm SWE) by melt_factor
            recovers the equivalent thermal forcing in °C·day, which is the
            currency the FGI uses. See Buckets._update_fgi().
        """
        if self.T <= 0:
            return 0.0

        pdd_avail   = self.melt_factor * self.T * dt    # [mm SWE]
        ros_avail   = _CP_LF * self.T * P               # [mm SWE] rain-on-snow
        total_avail = pdd_avail + ros_avail

        if total_avail <= self.Hwater:
            actual_melt = total_avail
            excess_dd   = 0.0
        else:
            actual_melt = self.Hwater
            # Leftover energy → °C·day equivalent for soil-thaw credit
            excess_dd = (total_avail - actual_melt) / self.melt_factor

        self.H_infiltrated += actual_melt
        self.Hwater        -= actual_melt
        return excess_dd


class Buckets(object):
    """
    Incorporates a list of reservoirs into a linear hierarchy that sends water
    either downwards or out to the surface. Reservoirs are ordered from top
    (nearest Earth's surface) to bottom (deepest groundwater); this order
    controls the direction of infiltration between layers.

    HydroRaVENS is designed as a daily-timestep model. This is a deliberate
    design choice: the physical parameterisations — degree-day snowmelt,
    Thornthwaite ET, and linear reservoir drainage — are climatological
    approximations that are well-founded at daily resolution but lose physical
    meaning at finer scales. The daily timestep is enforced in initialize().
    """

    def __init__(self, T_monthly_normals=None):
        """
        Initialize the watershed model.

        If using the ThorntwaiteChang2019 ET method, pass
        T_monthly_normals here so that the thermal index I and exponent
        a are computed once from climatological normals and remain fixed
        throughout the simulation.

        Parameters
        ----------
        T_monthly_normals : array-like of length 12, optional
            Long-term mean monthly temperatures (°C) used to compute the
            Thornthwaite thermal index I and exponent a per Chang et al.
            (2019), https://doi.org/10.1002/ird.2309. Required when
            evapotranspiration_method is 'ThorntwaiteChang2019'.
        """
        # Thornthwaite thermal index and exponent, per Chang et al. (2019)
        # https://doi.org/10.1002/ird.2309
        # I is climatologically imposed by the local normal temperature regime
        # and must remain fixed during simulation (not recomputed each timestep).
        if T_monthly_normals is not None:
            self.Chang_I = self._compute_Chang_I(T_monthly_normals)
            self.Chang_a = self._compute_Chang_a(self.Chang_I)

        # Frozen ground index (Molnau & Bissell 1983).  Disabled by default
        # (threshold = inf); overridden by snowmelt.fdd_threshold in YAML or
        # by run_and_score(fdd_threshold=).
        self.fdd_threshold = np.inf  # [°C·day]
        self._fgi          = 0.0    # current frozen ground index [°C·day]

    def _compute_Chang_I(self, T_monthly_normals):
        """
        Compute the Thornthwaite thermal index I from long-term monthly normal
        temperatures, per Chang et al. (2019), Eq. 1.
        https://doi.org/10.1002/ird.2309

        Parameters
        ----------
        T_monthly_normals : array-like, length 12
            Long-term mean monthly temperatures (°C). Negative values are
            treated as 0 per the Thornthwaite convention.

        Returns
        -------
        I : float
            Thermal index (dimensionless).
        """
        Tn = np.maximum(T_monthly_normals, 0)
        return np.sum((0.2 * Tn) ** 1.514)

    def _compute_Chang_a(self, I):
        """
        Compute the Thornthwaite exponent a from thermal index I, per
        Chang et al. (2019), Eq. 1.
        https://doi.org/10.1002/ird.2309

        Parameters
        ----------
        I : float
            Thermal index, as returned by _compute_Chang_I.

        Returns
        -------
        a : float
            Thornthwaite exponent (dimensionless).
        """
        return (6.75e-7 * I**3
                - 7.71e-5 * I**2
                + 1.7912e-2 * I
                + 0.49239)

    def export_Hlist(self):
        """
        Return the current water depths in all subsurface reservoirs.

        Useful for checkpointing reservoir state between a spin-up run
        and the main simulation, or for restarting a paused run.

        Returns
        -------
        list of float
            Water depth in each reservoir, ordered from shallowest
            (index 0) to deepest.
        """
        return [reservoir.Hwater for reservoir in self.reservoirs]

    def initialize(self, config_file=None, enforce_water_balance=None):
        """
        Set up the model from a YAML configuration file.

        Reads the configuration file, loads the input time series, builds
        the reservoir stack, instantiates snowpack if temperature data are
        present, optionally computes the water-year ET multiplier, and runs
        any requested spin-up cycles. Part of the CSDMS Basic Model Interface.

        Parameters
        ----------
        config_file : str, optional
            Path to the YAML configuration file. If None, all required
            values must be set on the object directly before calling
            update().
        enforce_water_balance : str or None, optional
            How to scale ET to close the water balance. When None (default),
            the value is read from ``general: enforce_water_balance`` in the
            YAML config, which itself defaults to ``'water-year'`` if absent.
            Accepted values:

            * ``'water-year'`` — Scale ET by a per-water-year multiplier so
              that P - Q - ET = 0 over each water year.
            * ``'global'`` — Scale ET by a single multiplier computed from
              sum(P - Q_obs) / sum(ET_raw) over the full record. No
              per-year overfitting; does not add hidden degrees of freedom
              to AIC comparisons.
            * ``'none'`` — Use raw ET without correction. Appropriate only
              when supplying trusted measured ET (e.g. eddy covariance).
              Using ``'none'`` with ThorntwaiteChang2019 will raise a
              warning because Thornthwaite ET carries large systematic
              biases.
        """
        if config_file is None:
            warnings.warn("No configuration file provided; all values needed "+
                          "for a model run therefore must be set independently.")

        # Parse YAML configuration file
        # And assign variables except for optimization bounds and plotting
        if config_file is not None:
            try:
                with open(config_file, "r") as yamlfile:
                    self.cfg = yaml.load(yamlfile, Loader=yaml.FullLoader)
            except FileNotFoundError:
                print("\nConfig file not found:", config_file, "\n")
                sys.exit(2)
            except yaml.YAMLError as e:
                print("\nCould not parse config file:", config_file, "\n", e)
                sys.exit(2)

        # Read input time series from the CSV path specified in the config
        self.hydrodata = pd.read_csv(
            self.cfg['timeseries']['datafile'],
            parse_dates=['Date'])

        # Set variables on reservoirs
        # First, check if all reservoirs have the same length
        for _key in self.cfg['reservoirs'].keys():
            if len(self.cfg['reservoirs'][_key]) == \
                    len(self.cfg['initial_conditions']['water_reservoir_effective_depths__mm']):
                pass
            else:
                raise ValueError(_key + ' within "reservoirs" contains a\n'+
                                 'different number of entries, implying'+
                                 'a different number of subsurface water\n'+
                                 'reservoirs, than '+
                                 'water_reservoir_effective_depths__mm'+
                                 ' within "initial_conditions".')

        # If all are the same length, then we will assign a number of reservoirs
        self.n_reservoirs = len(
            self.cfg['initial_conditions']['water_reservoir_effective_depths__mm'])
        # Using this, we will build a list of reservoir objects
        # and initialize them based on the provided inputs
        _pdm_H0   = self.cfg['reservoirs'].get('pdm_H0__mm',
                                               [None]  * self.n_reservoirs)
        _f_tile   = self.cfg['reservoirs'].get('tile_fractions',
                                               [0.0]   * self.n_reservoirs)
        _tau_tile = self.cfg['reservoirs'].get('tile_residence_times__days',
                                               [None]  * self.n_reservoirs)
        _recession_exp = self.cfg['reservoirs'].get('recession_exponents',
                                                   [1.0] * self.n_reservoirs)
        self.reservoirs = [
            Reservoir(
                t_recession    = self.cfg['reservoirs']['recession_timescales__days'][i],
                f_to_discharge = self.cfg['reservoirs']['exfiltration_fractions'][i],
                Hmax           = self.cfg['reservoirs']['maximum_effective_depths__mm'][i],
                pdm_H0         = _pdm_H0[i],
                H0             = self.cfg['initial_conditions']['water_reservoir_effective_depths__mm'][i],
                f_tile         = _f_tile[i],
                tau_tile       = _tau_tile[i],
            )
            for i in range(self.n_reservoirs)]
        for i, b_exp in enumerate(_recession_exp):
            self.reservoirs[i].recession_exponent = float(b_exp)

        # Check if bottom reservoir discharges all to river: conserve mass.
        # But allow through with a warning in case the user wants a
        # deep and non-discharging reservoir (although this could be set up
        # explicitly too).
        if self.reservoirs[-1].f_to_discharge < 1:
            warnings.warn("f_to_discharge of bottom water-storage layer < 1.\n"+
                          "You are not conserving mass.")

        # Set scalar variables based on yaml
        self.melt_factor         = self.cfg['snowmelt']['PDD_melt_factor']
        self.snow_insulation_k   = self.cfg['snowmelt'].get('snow_insulation_k',   0.0)
        self.fgi_decay_coeff     = self.cfg['snowmelt'].get('fgi_decay_coeff',     0.97)
        self.fdd_threshold       = self.cfg['snowmelt'].get('fdd_threshold',       np.inf)
        self.et_method = self.cfg['catchment']['evapotranspiration_method']
        if self.et_method == 'ThorntwaiteChang2019' and not hasattr(self, 'Chang_I'):
            n_years = len(self.hydrodata) / 365.25
            if n_years < 20:
                warnings.warn(
                    f"ThorntwaiteChang2019: monthly temperature normals were not provided "
                    f"and will be computed from the {n_years:.1f}-year input record. "
                    f"For short records this may not represent long-term climatology. "
                    f"Pass T_monthly_normals to Buckets() for reliable results.",
                    UserWarning, stacklevel=2,
                )
            T_normals = (self.hydrodata['Mean Temperature [C]']
                             .groupby(pd.DatetimeIndex(self.hydrodata['Date']).month)
                             .mean()
                             .reindex(range(1, 13))
                             .values)
            self.Chang_I = self._compute_Chang_I(T_normals)
            self.Chang_a = self._compute_Chang_a(self.Chang_I)
        self.water_year_start_month = self.cfg['catchment']['water_year_start_month']
        self.drainage_basin_area__km2 = self.cfg['catchment']['drainage_basin_area__km2']
        self.baseflow_Q = self.cfg['catchment'].get('baseflow_Q', 0.0)

        # Module enable/disable flags — read from config, default to on
        # (except direct_runoff, which defaults to off).
        _modules = self.cfg.get('modules', {})
        self.use_snowpack        = _modules.get('snowpack',        True)
        self.use_frozen_ground   = _modules.get('frozen_ground',   True)
        self.use_rain_on_snow    = _modules.get('rain_on_snow',    True)
        self.use_direct_runoff   = _modules.get('direct_runoff',   False)
        self.use_dtr_fgi_decay   = _modules.get('dtr_fgi_decay',   True)
        self.use_et_water_stress  = _modules.get('et_water_stress',  False)
        self.use_et_reservoir_draw = _modules.get('et_reservoir_draw', False)
        if self.use_et_water_stress and self.use_et_reservoir_draw:
            warnings.warn(
                "et_water_stress and et_reservoir_draw are mutually exclusive. "
                "et_reservoir_draw will be used; et_water_stress ignored.",
                UserWarning, stacklevel=2,
            )
            self.use_et_water_stress = False
        self.et_scale = 1.0   # global ET multiplier; calibrated when et_water_stress or
        #                       et_reservoir_draw (with enforce_water_balance='none') is on
        # Fraction of ET_pot drawn from reservoir 0 (shallow); 1-et_alpha from reservoir 1.
        # Read from general: et_alpha in config YAML; override via run_and_score(et_alpha=).
        self.et_alpha = self.cfg['general'].get('et_alpha', 1.0)
        # Wilting-point threshold for soil reservoir ET draw [mm].
        # wp_soil_sigma > 0 enables the spatially variable (Gaussian CDF) form.
        self.wp_soil = 0.0
        self.wp_soil_sigma = 0.0

        # Check if there is a mean temperature column for snowpack.
        # If not, note that no snowpack processes will be included
        self.has_snowpack = (self.use_snowpack and
                             'Mean Temperature [C]' in self.hydrodata.columns)

        # Detect optional daily T_min / T_max for DTR-based FGI decay.
        self._has_trange = (self.use_dtr_fgi_decay and
                            'Minimum Temperature [C]' in self.hydrodata.columns and
                            'Maximum Temperature [C]' in self.hydrodata.columns)
        if self.has_snowpack:
            # Instantiate snowpack
            self.snowpack = Snowpack(self.melt_factor)  # allow changes to melt factor later
        elif 'Mean Temperature [C]' in self.hydrodata.columns and not self.use_snowpack:
            pass  # snowpack deliberately disabled via modules config
        else:
            warnings.warn('"Mean Temperature [C]" has not been set. '
                          'No snowpack processes will be simulated.')

        # How many times to loop the full time series for the spin-up
        # Maybe I should permit a more sophisticated spin-up at some point!
        self.n_spin_up_cycles = self.cfg['general']['spin_up_cycles']

        # Resolve enforce_water_balance: keyword argument takes precedence over YAML,
        # which defaults to 'water-year' if the key is absent.
        if enforce_water_balance is None:
            enforce_water_balance = self.cfg['general'].get('enforce_water_balance',
                                                            'water-year')
        # Normalize legacy boolean values from old YAML configs or API calls.
        if enforce_water_balance is True:
            enforce_water_balance = 'water-year'
        elif enforce_water_balance is False:
            enforce_water_balance = 'none'
        self.enforce_water_balance = enforce_water_balance

        # Fraction of positive daily recharge that bypasses the reservoir
        # cascade and exits directly as runoff.  Conceptually inspired by
        # Hortonian (infiltration-excess) overland flow, but at a daily
        # timestep rainfall intensity is unavailable, so the fraction cannot
        # be a rigorous physical representation -- except in extreme events
        # where intense rainfall dominates the daily total.  In practice it
        # acts as a calibrated fast-bypass fraction, off by default.
        self.direct_runoff_fraction = self.cfg['general'].get(
            'direct_runoff_fraction', 0.0)

        # Initial conditions if resuming from prior run
        if self.has_snowpack:
            self.snowpack.Hwater = self.cfg['initial_conditions']['snowpack__mm_SWE']
        # Reservoir H0 values are set in the list comprehension above.

        # Enforce the daily timestep. HydroRaVENS is a daily model by design:
        # degree-day snowmelt, Thornthwaite ET, and linear reservoir drainage
        # are all daily-scale parameterisations.
        if (self.hydrodata['Date'].diff()[1:] == pd.Timedelta('1 day')).all():
            self.dt = 1.
        else:
            raise ValueError(
                "HydroRaVENS requires a continuous daily time series "
                "(exactly 1-day intervals throughout). Sub-daily or "
                "irregular timesteps are not supported; the physical "
                "parameterisations (degree-day snowmelt, Thornthwaite ET, "
                "linear reservoir drainage) are daily-scale approximations."
            )

        # Compute specific discharge from data
        self.hydrodata['Specific Discharge [mm/day]'] = (
            self.hydrodata['Discharge [m^3/s]']
            / (self.drainage_basin_area__km2*1E3) * 86400)

        # Create columns for model output
        self.hydrodata['Specific Discharge (modeled) [mm/day]'] = pd.NA
        self.hydrodata['Snowpack (modeled) [mm SWE]'] = pd.NA
        self.hydrodata['Subsurface storage (modeled total) [mm]'] = pd.NA

        # Start out at first timestep
        # Could modify this to pick up a run in the middle
        # Or start at the beginning of a water year
        # for example
        self._timestep_i = self.hydrodata.index[0]

        # Carry-over of any water deficit from the previous timestep that the
        # deepest reservoir could not satisfy (ET > P + all storage).  This is
        # the unpaid debt passed forward one step; distinct from
        # Reservoir.H_deficit and Snowpack.H_deficit, which are per-timestep only.
        self.H_deficit_carry = 0.

        # Compute the water years based on the input month for
        # water-year rollover
        self.compute_water_year()

        # Compute ET, optionally scaling to close the water balance.
        self.global_et_multiplier = 1.0   # default; overwritten below if global mode
        if self.enforce_water_balance == 'global':
            self.compute_global_ET_multiplier()
        elif self.enforce_water_balance == 'water-year':
            self.compute_ET_multiplier()
        elif self.et_method == 'ThorntwaiteChang2019':
            warnings.warn(
                "enforce_water_balance='none' with ThorntwaiteChang2019: Thornthwaite ET "
                "will not be rescaled to close the water balance. "
                "Thornthwaite ET carries large systematic biases; omitting "
                "the correction is likely to produce significant mass-balance "
                "errors. Consider enforce_water_balance='water-year' or 'global', "
                "or supply measured ET via evapotranspiration_method: datafile."
            )
        self.compute_ET()

        # Model spin-up, if requested
        for _ in range(self.n_spin_up_cycles):
            self.run()  # Spin-up; run() resets _timestep_i each call

    def compute_water_year(self):
        """
        Assign a water-year label to each row in self.hydrodata.

        Adds a 'Water Year' column. A water year begins on
        water_year_start_month and is labelled by the calendar year in
        which it ends. For example, with a start month of October (USGS
        convention), 1 Oct 2020 – 30 Sep 2021 is water year 2021.

        When water_year_start_month is 1 (January), the water year equals
        the calendar year and no offset is applied.
        """
        self.hydrodata['Water Year'] = pd.DatetimeIndex(self.hydrodata['Date']).year
        if self.water_year_start_month > 1:
            self.hydrodata['Water Year'] += (
                pd.DatetimeIndex(self.hydrodata['Date']).month
                >= self.water_year_start_month
            )

    def compute_global_ET_multiplier(self):
        """
        Compute a single ET scale factor to close the long-term water balance.

        Computes scale = sum(P - Q_obs) / sum(ET_raw) over all days where
        P, ET, and Q_obs are all finite, and stores it as self.global_et_multiplier.
        Unlike compute_ET_multiplier(), which fits one multiplier per water
        year, this uses a single ratio for the full record and does not
        overfit interannual variability. Appropriate when enforce_water_balance
        is set to 'global'.
        """
        if self.et_method == 'datafile':
            _raw_ET = np.asarray(self.hydrodata['Evapotranspiration [mm/day]'])
        else:
            _raw_ET = np.asarray(self.evapotranspiration_Chang2019())

        # Mask to days where all three water-balance terms are finite so that
        # P_sum, ET_sum, and Q_sum are computed over the same day-set.
        mask = (self.hydrodata['Specific Discharge [mm/day]'].notna()
                & self.hydrodata['Precipitation [mm/day]'].notna()
                & np.isfinite(_raw_ET))
        P_sum  = float(self.hydrodata.loc[mask, 'Precipitation [mm/day]'].sum())
        Q_sum  = float(self.hydrodata.loc[mask, 'Specific Discharge [mm/day]'].sum())
        ET_sum = float(_raw_ET[mask].sum())

        if ET_sum <= 0:
            warnings.warn("Global ET sum is zero or negative; multiplier set to 1.0.",
                          UserWarning, stacklevel=2)
            self.global_et_multiplier = 1.0
        else:
            self.global_et_multiplier = (P_sum - Q_sum) / ET_sum
            if self.global_et_multiplier <= 0:
                warnings.warn(
                    f"Global ET multiplier = {self.global_et_multiplier:.3f} (≤ 0); "
                    "P − Q ≤ 0 over the full record. Setting to 1.0.",
                    UserWarning, stacklevel=2)
                self.global_et_multiplier = 1.0

    def compute_ET_multiplier(self):
        """
        Compute per-water-year ET scaling factors to enforce water balance.

        For each water year, computes the ratio of required ET (P - Q) to
        measured or computed ET, and stores this as 'ET multiplier' in
        self.hydrodata_WY_means. This multiplier is later applied in
        compute_ET() to scale ET so that P - Q - ET = 0 over each water year.
        """
        # Originally used "sum", but then used "mean" so the headers would
        # still be sensible
        self.hydrodata_WY_means = self.hydrodata.groupby(
            self.hydrodata['Water Year']).mean(numeric_only=True)
        # Not needed, but no real harm in calculating
        self.hydrodata_WY_means['Runoff ratio'] = (
            self.hydrodata_WY_means['Specific Discharge [mm/day]']
            / self.hydrodata_WY_means['Precipitation [mm/day]'])
        _ET_required = -(self.hydrodata_WY_means['Specific Discharge [mm/day]'] -
                         self.hydrodata_WY_means['Precipitation [mm/day]'])
        self.hydrodata_WY_means['ET multiplier'] = (
            _ET_required / self.hydrodata_WY_means['Evapotranspiration [mm/day]'])

        _bad_wy = self.hydrodata_WY_means.index[
            self.hydrodata_WY_means['ET multiplier'] <= 0]
        if len(_bad_wy) > 0:
            warnings.warn(
                f"ET multiplier <= 0 in water year(s) {list(_bad_wy)}. "
                "Annual discharge exceeds precipitation for those years; "
                "scaled ET will be zero or negative (water-generating). "
                "Check gauge data or consider removing those years."
            )

    def compute_ET(self):
        """
        Build the ET time series used in the model.

        Obtains raw daily ET from the input data file or the Thornthwaite–Chang
        2019 equation (see evapotranspiration_Chang2019()). Five modes:

        1. et_water_stress=True: multiply raw ET by self.et_scale (a single
           global calibration parameter, analogous to PDD_melt_factor for snow).
           The per-water-year water-balance multiplier is bypassed because the
           dynamic stress factor applied in _et_stress_factor() at each time step
           would make the pre-computed per-year correction inconsistent.
        2. enforce_water_balance='water-year' (default, stress off): raw ET is
           multiplied by a per-water-year multiplier so P - Q - ET = 0 each year.
        3. enforce_water_balance='global' (stress off): raw ET multiplied by a
           single full-record multiplier from sum(P - Q_obs) / sum(ET_raw).
        4. enforce_water_balance='none' (stress off): raw ET used directly.
        5. et_reservoir_draw=True + enforce_water_balance='none': same as mode 1
           (et_scale applied), but ET is drawn from reservoirs after the cascade
           rather than pre-subtracted from precipitation. With 'global' WB, falls
           through to mode 3.

        The result is stored as 'ET for model [mm/day]' in self.hydrodata.
        """
        if self.et_method == 'datafile':
            _raw_ET = self.hydrodata['Evapotranspiration [mm/day]']
        elif self.et_method == 'ThorntwaiteChang2019':
            _raw_ET = self.evapotranspiration_Chang2019()
        else:
            raise ValueError('evapotranspiration_method must be "datafile" or '+
                             '"ThorntwaiteChang2019".')

        if self.use_et_water_stress or (self.use_et_reservoir_draw and
                                         self.enforce_water_balance == 'none'):
            # et_scale is the sole water-balance closure mechanism in both
            # et_water_stress mode and et_reservoir_draw without WB enforcement.
            self.hydrodata['ET for model [mm/day]'] = (
                np.asarray(_raw_ET) * self.et_scale)
        elif self.enforce_water_balance == 'global':
            self.hydrodata['ET for model [mm/day]'] = (
                np.asarray(_raw_ET) * self.global_et_multiplier)
        elif self.enforce_water_balance == 'water-year':
            # Merge per-water-year multiplier into hydrodata, then apply.
            # Drop any previous 'ET multiplier' column first so that calling
            # compute_ET() more than once (e.g. after a module flag override)
            # does not produce duplicate _x/_y suffix columns from pandas merge.
            # Use .to_numpy() to multiply by position rather than pandas index
            # so that any index reset from the merge cannot silently misalign rows.
            if 'ET multiplier' in self.hydrodata.columns:
                self.hydrodata = self.hydrodata.drop(columns=['ET multiplier'])
            self.hydrodata = self.hydrodata.merge(
                self.hydrodata_WY_means['ET multiplier'],
                on='Water Year')
            nan_wy = (self.hydrodata_WY_means['ET multiplier']
                          .index[self.hydrodata_WY_means['ET multiplier'].isna()]
                          .tolist())
            if nan_wy:
                warnings.warn(
                    f"ET multiplier is NaN for water year(s) {nan_wy} "
                    f"(no discharge observations). Raw ET used for those years "
                    f"(enforce_water_balance ineffective).",
                    UserWarning, stacklevel=2,
                )
                self.hydrodata['ET multiplier'] = (
                    self.hydrodata['ET multiplier'].fillna(1.0))
            self.hydrodata['ET for model [mm/day]'] = (
                np.asarray(_raw_ET) * self.hydrodata['ET multiplier'].to_numpy())
        else:
            self.hydrodata['ET for model [mm/day]'] = np.asarray(_raw_ET)

    def _et_stress_factor(self):
        """
        Water-availability multiplier applied to potential ET each time step.

        Returns 1 - exp(-H_shallow / H0), where H_shallow is the current water
        depth in the shallowest reservoir and H0 is its PDM characteristic
        storage depth (pdm_H0).  The multiplier is zero when the reservoir is
        empty and approaches 1 as it fills, so actual ET = potential ET * factor.

        Returns 1.0 (no stress) when et_water_stress is disabled or when the
        shallow reservoir has no pdm_H0 set.
        """
        if not self.use_et_water_stress:
            return 1.0
        H0 = self.reservoirs[0].pdm_H0
        if H0 is None:
            return 1.0
        return 1.0 - np.exp(-max(self.reservoirs[0].Hwater, 0.0) / H0)

    def _draw_et_from_reservoirs(self, ET_pot):
        """
        Remove actual ET from reservoirs after cascade drainage.

        ET_pot is partitioned between reservoir 0 (shallow) and reservoir 1
        (soil) by et_alpha.  Each draw is capped at available storage so
        Hwater never goes negative; unmet demand is implicitly lost as
        water-stress reduction of actual ET.

        Parameters
        ----------
        ET_pot : float
            Potential ET for this time step [mm/day], already scaled by
            et_scale or the global water-balance multiplier.
        """
        fractions = [self.et_alpha, 1.0 - self.et_alpha]
        for i, frac in enumerate(fractions):
            if i >= len(self.reservoirs):
                break
            demand = frac * ET_pot
            H = self.reservoirs[i].Hwater
            if i == 1 and self.wp_soil > 0.0:
                # Soil reservoir: scale demand by fraction of catchment above WP.
                if self.wp_soil_sigma > 0.0:
                    # Gaussian CDF: spatially variable WP — σ→0 recovers hard threshold.
                    f_avail = 0.5 * (1.0 + math.erf((H - self.wp_soil)
                                                     / (self.wp_soil_sigma * math.sqrt(2.0))))
                    demand = demand * f_avail
                else:
                    # Hard threshold: no ET extraction below wp_soil.
                    if H <= self.wp_soil:
                        continue
                    H = H - self.wp_soil   # only water above WP is available
            actual = min(demand, max(0.0, H))
            self.reservoirs[i].Hwater -= actual

    def _compute_snowpack(self, time_step):
        """
        Update the snowpack for one timestep; return excess melt energy.

        Sets temperature, applies net water input, then calls melt() with
        the raw precipitation so rain-on-snow sensible heat is included.
        Updates self.H_deficit_carry from the snowpack before returning.

        Returns
        -------
        excess_dd : float
            Leftover melt energy [°C·day] after SWE is fully depleted.
            Pass to _update_fgi() to credit toward frozen-soil thawing.
        """
        T = self.hydrodata['Mean Temperature [C]'][time_step]
        P = self.hydrodata['Precipitation [mm/day]'][time_step]
        self.snowpack.set_temperature(T)
        if self.use_et_reservoir_draw:
            _sp_recharge = P + self.H_deficit_carry
        else:
            _sp_recharge = (P
                            - self.hydrodata['ET for model [mm/day]'][time_step]
                            * self._et_stress_factor()
                            + self.H_deficit_carry)
        self.snowpack.recharge(_sp_recharge)
        excess_dd = self.snowpack.melt(self.dt,
                                       P=(P if self.use_rain_on_snow else 0.0))
        self.H_deficit_carry = self.snowpack.H_deficit
        return excess_dd

    def _update_fgi(self, time_step, excess_dd=0.0):
        """
        Update the frozen ground index; flag top reservoir as frozen if needed.

        FGI(t) = max(0, A_t · FGI(t-1) - T_eff - excess_dd)
          A_t   = daily FGI decay coefficient (climate-dependent; see below)
          T_eff = T_mean · exp(-snow_insulation_k · SWE)
          T_mean < 0  → FGI rises  (freezing degree-days accumulate)
          T_mean > 0  → FGI falls  (warm air thaws)
          excess_dd   → additional thaw credited from leftover snowmelt
                        energy [°C·day] = leftover mm SWE / melt_factor

        DTR-based decay (when Minimum/Maximum Temperature columns present):
          On days entirely below freezing (T_max ≤ 0), A_t = 1.0 — no
          sub-daily thawing, FGI accumulates unimpeded. When the diurnal
          cycle straddles 0°C (T_min < 0 < T_max), a fraction of the day
          is above freezing and drives partial thawing:

              f_above = T_max / (T_max - T_min)   [linear diurnal cycle]
              A_t     = 1 - (1 - fgi_decay_coeff) · f_above

          When f_above → 1, A_t → fgi_decay_coeff (M&B 0.97 maximum decay).
          When f_above = 0, A_t = 1.0 (no passive decay).
          This naturally gives near-unity A for continental climates (where
          T_max rarely crosses 0 during cold spells) and M&B-like decay for
          maritime climates with frequent freeze-thaw oscillations.

        Fallback (T_min / T_max absent): A_t = fgi_decay_coeff (constant,
          original M&B behaviour).

        When FGI exceeds fdd_threshold, the top reservoir's f_to_discharge
        is set to 1.0 so all drainage becomes direct runoff, simulating
        frozen-soil blockage of deep infiltration.

        References
        ----------
        Molnau & Bissell (1983) https://westernsnowconference.org/bibliography/1983Molnau.pdf
            (Western Snow Conference proceedings — original source for the
            FGI formulation and the exponential snow insulation factor)
        Shanley & Chalmers (1999) doi:10.1002/(SICI)1099-1085(199909)13:12/13
            <1843::AID-HYP879>3.0.CO;2-G
        Dunne & Black (1971) doi:10.1029/WR007i005p01160
        Snow insulation parameterisation (exponential form, original: M&B 1983):
            LISFLOOD: van der Knijff et al. (2010) doi:10.1080/02626660902852568
            GSSHA: Downer & Ogden (2004) doi:10.1061/(ASCE)1084-0699(2004)9:3(254)

        Parameters
        ----------
        time_step : int
            Current row index in self.hydrodata.
        excess_dd : float, optional
            Degree-day equivalent of leftover melt energy from
            _compute_snowpack() [°C·day]. Reduces FGI alongside air
            temperature. Default 0 (temperature-only FGI, per Molnau &
            Bissell). Computed as leftover_mm_SWE / melt_factor, where
            melt_factor (mm SWE °C⁻¹ day⁻¹) converts the residual melt
            depth back to the degree-day units that the FGI operates in.
            See Snowpack.melt(). Not scaled by snow insulation because
            meltwater delivers heat directly to the soil surface.

        Returns
        -------
        f0 : float
            Calibrated f_to_discharge of the top reservoir, saved before any
            frozen-ground override. Restore it after the discharge loop.
        """
        f0 = self.reservoirs[0].f_to_discharge
        if not self.use_frozen_ground or np.isinf(self.fdd_threshold):
            return f0

        if 'Mean Temperature [C]' not in self.hydrodata.columns:
            raise ValueError(
                "fdd_threshold is set but 'Mean Temperature [C]' is missing "
                "from the input data. FGI requires temperature data."
            )
        T     = self.hydrodata['Mean Temperature [C]'][time_step]
        T_eff = T * np.exp(-self.snow_insulation_k * self.snowpack.Hwater)

        if self._has_trange:
            T_max = self.hydrodata['Maximum Temperature [C]'][time_step]
            T_min = self.hydrodata['Minimum Temperature [C]'][time_step]
            DTR   = T_max - T_min
            if DTR > 0 and T_max > 0:
                f_above = min(1.0, T_max / DTR)
                A_t = 1.0 - (1.0 - self.fgi_decay_coeff) * f_above
            else:
                A_t = 1.0  # entirely below freezing; no sub-daily thawing
        else:
            A_t = self.fgi_decay_coeff

        self._fgi = max(0.0, A_t * self._fgi - T_eff - excess_dd)
        if self._fgi > self.fdd_threshold:
            self.reservoirs[0].f_to_discharge = 1.0
        return f0

    def update(self, time_step=None):
        """
        Advance the model by one time step.

        Routes precipitation minus ET through the snowpack (if present) and
        then through each subsurface reservoir in order from shallowest to
        deepest. Stores modeled specific discharge, snowpack SWE, and total
        subsurface storage in self.hydrodata for the current time step.
        Part of the CSDMS Basic Model Interface.

        Parameters
        ----------
        time_step : int, optional
            Index into self.hydrodata for the time step to update. If None,
            uses and then increments the internal counter self._timestep_i.
        """

        if time_step is None:
            time_step = self._timestep_i
            # Advance internal variable if external time step is not selected.
            # This should be a different variable and therefore not
            # modify the value of "time_step" by reference.
            self._timestep_i += 1

        # Skip timesteps with missing forcing: leave reservoir states unchanged
        # and record NaN output so the scoring mask excludes these days.
        _P  = self.hydrodata['Precipitation [mm/day]'][time_step]
        _ET = self.hydrodata['ET for model [mm/day]'][time_step]
        _skip = not np.isfinite(_P) or not np.isfinite(_ET)
        if not _skip and self.has_snowpack:
            _skip = not np.isfinite(
                self.hydrodata['Mean Temperature [C]'][time_step])
        if _skip:
            self.hydrodata.at[time_step,
                              'Specific Discharge (modeled) [mm/day]'] = np.nan
            if self.has_snowpack:
                self.hydrodata.at[time_step, 'Snowpack (modeled) [mm SWE]'] = (
                    self.snowpack.Hwater)
            self.hydrodata.at[time_step,
                              'Subsurface storage (modeled total) [mm]'] = (
                np.sum([res.Hwater for res in self.reservoirs])
                + np.sum([res.tile_res.Hwater for res in self.reservoirs
                          if res.tile_res is not None]))
            return

        excess_dd = self._compute_snowpack(time_step) if self.has_snowpack else 0.0
        f0 = self._update_fgi(time_step, excess_dd)

        qi = 0.0
        for i in range(len(self.reservoirs)):
            if i == 0:
                if self.has_snowpack:
                    _recharge = (self.snowpack.H_infiltrated
                                 + self.H_deficit_carry)
                else:
                    if self.use_et_reservoir_draw:
                        _recharge = (
                            self.hydrodata['Precipitation [mm/day]'][time_step] +
                            self.H_deficit_carry)
                    else:
                        _recharge = (
                            self.hydrodata['Precipitation [mm/day]'][time_step] -
                            self.hydrodata['ET for model [mm/day]'][time_step]
                            * self._et_stress_factor() +
                            self.H_deficit_carry)
                # Hortonian-inspired bypass: fraction exits without entering reservoirs.
                _q_direct = (max(0.0, _recharge) * self.direct_runoff_fraction
                             if self.use_direct_runoff else 0.0)
                qi += _q_direct
                self.reservoirs[i].recharge(_recharge - _q_direct)
            else:
                # Let water infiltrate to lower layers effectively
                # instantaneously; this isn't quite realistic, but
                # should be a simpler approach for parameter calibration
                # (Plus, this is just the water that did exit that above
                # container, which is already free to discharge, so this
                # seems more self-consistent.)
                # The amount of infiltrated water from above could be
                # negative; this represents ET in excess of what the
                # unsaturated zone ("soil zone"; top reservoir) holds.
                # Deeper loss of water could be due to plants tapping into
                # groundwater, direct lake evaporation, etc. -- or related
                # to this model not being physical or distributed, so just
                # needing to balance mass.
                self.reservoirs[i].recharge(
                    self.reservoirs[i-1].H_infiltrated
                    + self.reservoirs[i-1].H_deficit)
            self.reservoirs[i].discharge(self.dt)
            qi += self.reservoirs[i].H_discharge

        self.reservoirs[0].f_to_discharge = f0
        self.H_deficit_carry = self.reservoirs[-1].H_deficit

        if self.use_et_reservoir_draw:
            self._draw_et_from_reservoirs(
                self.hydrodata['ET for model [mm/day]'][time_step])

        self.hydrodata.at[time_step, 'Specific Discharge (modeled) [mm/day]'] = qi
        if self.has_snowpack:
            self.hydrodata.at[time_step, 'Snowpack (modeled) [mm SWE]'] = self.snowpack.Hwater
        self.hydrodata.at[time_step, 'Subsurface storage (modeled total) [mm]'] = (
            np.sum([res.Hwater for res in self.reservoirs])
            + np.sum([res.tile_res.Hwater for res in self.reservoirs
                      if res.tile_res is not None]))

    def evapotranspiration_Chang2019(self, Tmax=None, Tmin=None, photoperiod=None,
                                    k=0.69):
        """
        Modified daily Thornthwaite ET₀ equation.

        Chang et al. (2019), Eq. 1–4. https://doi.org/10.1002/ird.2309

        Parameters
        ----------
        Tmax : array-like
            Daily maximum temperature (°C).
        Tmin : array-like
            Daily minimum temperature (°C).
        photoperiod : array-like
            Photoperiod N (hours), computed from latitude and Julian day
            per Allen et al. (1998), Eqs. 2–4 of Chang et al. (2019).
        k : float
            Calibration coefficient in the T_ef formula. Default 0.69,
            recommended by Pereira & Pruitt (2004) for daily ET₀
            (https://doi.org/10.1016/j.agrformet.2004.01.005).
            Use 0.72 for monthly ET₀ per Camargo et al. (1999).

        Returns
        -------
        ET0 : array-like
            Daily reference evapotranspiration (mm day⁻¹).
        """
        if Tmax is None:
            Tmax = self.hydrodata['Maximum Temperature [C]']
        if Tmin is None:
            Tmin = self.hydrodata['Minimum Temperature [C]']
        if photoperiod is None:
            photoperiod = self.hydrodata['Photoperiod [hr]']

        Tef = 0.5 * k * (3 * Tmax - Tmin)
        C = photoperiod / 360.

        quadratic  = C * (-415.85 + 32.24 * Tef - 0.43 * Tef**2)
        power_law  = 16. * C * (10. * Tef / self.Chang_I) ** self.Chang_a

        ET0 = np.where(np.isnan(Tef), np.nan,
               np.where(Tef > 26,   quadratic,
               np.where(Tef > 0,    power_law,
                                    0.)))
        return ET0

    def run(self):
        """
        Advance the model through all time steps in self.hydrodata.

        Resets the internal time counter to the first row before iterating,
        so run() is safe to call more than once (e.g. spin-up then main run).
        Captures storage at the start of the run for check_mass_balance().
        Part of the CSDMS Basic Model Interface.
        """
        self._timestep_i = self.hydrodata.index[0]
        self._run_initial_storage = (
            sum(res.Hwater for res in self.reservoirs)
            + (self.snowpack.Hwater if self.has_snowpack else 0.0)
        )
        for _ in self.hydrodata.index:
            self.update()

    def finalize(self):
        """
        Report model skill and display output plots.

        Calls compute_NSE(verbose=True) to print the Nash–Sutcliffe Efficiency
        to stdout, then calls plot() to display a time-series comparison of
        observed and modeled specific discharge. Part of the CSDMS Basic Model
        Interface.
        """
        # Goodness of fit
        # Add options to print and/or save values later
        self.compute_NSE(verbose=True)
        # Plot
        # Add flag for plotting (or not) later
        self.plot()

    def plot(self):
        """
        Display a time-series comparison of precipitation and specific discharge.

        Produces a dual-axis figure: precipitation as a bar chart on the left
        axis and both observed and modeled specific discharge as line plots on
        the right axis.
        """
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%y-%m-%d'))
        plt.xlabel('Date', fontsize=14)
        plt.xticks(rotation=45, horizontalalignment='right')
        plt.ylabel('Precipitation [mm/day]', fontsize=14, color='C0')
        plt.bar(self.hydrodata['Date'].values,
                height=self.hydrodata['Precipitation [mm/day]'].values/self.dt,
                width=1., align='center', label='Precipitation [mm/day]',
                linewidth=0, color='C0', alpha=0.5)  # C0 is the default bar-plot color
        plt.twinx()
        plt.plot(self.hydrodata['Date'].values,
                 self.hydrodata['Specific Discharge [mm/day]'].values,
                 'royalblue', label='Data', linewidth=2, alpha=0.8)
        plt.plot(self.hydrodata['Date'].values,
                 self.hydrodata['Specific Discharge (modeled) [mm/day]'].values,
                 'k', label='Model', linewidth=2, alpha=0.8)
        plt.ylim(0, plt.ylim()[-1])
        plt.legend(title='Specific Discharge', fontsize=11,
                   title_fontsize=11, labelcolor='linecolor')
        plt.ylabel('Specific Discharge [mm/day]', fontsize=14, color='0.3')
        plt.tight_layout()
        plt.show()

    def check_mass_balance(self, time_step=None):
        """
        Compute the mass-balance discrepancy at a given time step.

        Compares cumulative inputs (P - ET) from the start of the record
        through time_step with cumulative outputs (discharge) plus current
        storage (snowpack + subsurface reservoirs) and any carried-over
        deficit. Returns the excess mass still in the model; a value near
        zero indicates good mass conservation.

        Parameters
        ----------
        time_step : int, optional
            Row index in self.hydrodata at which to evaluate the balance.
            Defaults to the last row.

        Returns
        -------
        excess_mass_in_model : float
            Excess water remaining in the model budget (mm). Should be ~0
            for a mass-conserving run.
        """
        if time_step is None:
            time_step = self.hydrodata.index[-1]
        # Additions equals discharge out; set up this way, and can check.
        total_additions = \
            self.hydrodata['Precipitation [mm/day]'][:time_step+1].sum() \
            - self.hydrodata['ET for model [mm/day]'][:time_step+1].sum()
        # Storage reservoirs; snowpack is 0 when not simulated
        snow_storage = (self.hydrodata['Snowpack (modeled) [mm SWE]'][time_step]
                        if self.has_snowpack else 0.)
        subsurface_storage = self.hydrodata['Subsurface storage (modeled total) [mm]'][time_step]
        # Mass removal
        outlet_discharge = self.hydrodata[
            'Specific Discharge (modeled) [mm/day]'][:time_step+1].sum()
        # Unpaid water deficit carried forward from the last timestep
        deficit = self.H_deficit_carry

        # Initial storage at the start of the last run() call (not at initialize()
        # time, since spin-up changes storage before the scored run begins).
        initial_storage = getattr(self, '_run_initial_storage', 0.0)

        # Discrepancy: inputs = outputs + ΔS, so excess ≈ 0 when mass is conserved.
        excess_mass_in_model = (outlet_discharge + subsurface_storage
                                + snow_storage - total_additions + deficit
                                - initial_storage)

        return excess_mass_in_model

    def compute_NSE(self, return_nse=True, verbose=False):
        """
        Compute the Nash–Sutcliffe Efficiency of the discharge simulation.

        Compares modeled and observed specific discharge for all time steps
        where both values are non-missing. Stores the result as self.NSE.

        Parameters
        ----------
        return_nse : bool, optional
            If True (default), return the NSE value.
        verbose : bool, optional
            If True, print the NSE value to stdout.

        Returns
        -------
        NSE : float or None
            Nash–Sutcliffe Efficiency coefficient. Returns None if
            return_nse is False. A value of 1 indicates perfect agreement;
            values below 0 indicate the model performs worse than the
            observed-mean predictor.
        """

        q_data  = self.hydrodata['Specific Discharge [mm/day]']
        q_model = self.hydrodata['Specific Discharge (modeled) [mm/day]']

        # Calculate NSE
        _realvalue = ~q_model.isna() & ~q_data.isna()
        NSE_num = np.sum((q_model[_realvalue] - q_data[_realvalue])**2)
        NSE_denom = np.sum((q_data[_realvalue] -
                            np.mean(q_data[_realvalue]))**2)
        if np.sum(~_realvalue):
            print("Excluded", np.sum(~_realvalue), "no-data points from NSE calculation")

        self.NSE = 1 - NSE_num / NSE_denom

        if verbose:
            print("NSE:", self.NSE)

        if return_nse:
            return self.NSE


def main():
    """Command-line entry point: parse ``-y``/``--configfile`` and run the model."""
    parser = argparse.ArgumentParser(
        description='Pass the configuration file path to run hydroRaVENS.')
    parser.add_argument('-y', '--configfile', type=str,
                        help='YAML file from which all inputs are read.')

    # Parse args if anything is passed.
    # If nothing is passed, then print help and exit.
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    b = Buckets()
    b.initialize(args.configfile)
    b.run()
    b.finalize()


if __name__ == "__main__":
    main()
