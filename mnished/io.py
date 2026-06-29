"""
mnished.io
~~~~~~~~~~

The MNiShed **input contract**: the authoritative, code-grounded specification of
the forcing CSV and config YAML that :meth:`mnished.Buckets.initialize` consumes,
plus a fast pre-flight validator.

Two roles:

* **Documents the data spec.** :data:`FORCING_COLUMNS` and
  :data:`CONFIG_SECTIONS` are a catalog (name, units, role) of every input MNiShed
  reads — a human- and producer-facing reference (e.g. for ``mnished-builder``).
* **Validates inputs up front, and *is* the enforcement.** Which columns/sections
  are required is config-dependent logic that lives in the validator functions
  (:func:`required_forcing_columns`, :func:`recommended_forcing_columns`,
  :func:`validate_config`, :func:`validate_forcing`) — the single source of truth
  for the conditions; a test keeps the catalog and the functions in sync.
  :func:`validate_inputs` checks a config + forcing pair and returns a
  :class:`ValidationReport` listing *all* problems at once.

This is a **contract**-level check (are the right columns and config sections
present, given the chosen options?). It deliberately does **not** duplicate the
deep model-consistency validation :meth:`~mnished.Buckets.initialize` performs
(sub-catchment area fractions summing to one, lake outflow coefficients, junction
types, mass conservation); run a trial ``initialize`` for that. The two are
complementary: this one is fast, needs no model run, and reports the whole list.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

# Canonical forcing column headers.
DATE = "Date"
PRECIP = "Precipitation [mm/day]"
DISCHARGE = "Discharge [m^3/s]"
TMEAN = "Mean Temperature [C]"
TMIN = "Minimum Temperature [C]"
TMAX = "Maximum Temperature [C]"
PHOTOPERIOD = "Photoperiod [hr]"
ET = "Evapotranspiration [mm/day]"

ET_METHODS = ("datafile", "ThornthwaiteChang2019")


@dataclass(frozen=True)
class ForcingColumn:
    """One column of the forcing CSV (see :data:`FORCING_COLUMNS`)."""

    name: str
    units: str
    description: str


#: Documentation catalog of every forcing column MNiShed may read, with its units
#: and role. This is a *catalog*, not the enforcement: which columns are required
#: under which options is config-dependent logic and lives in
#: :func:`required_forcing_columns` / :func:`recommended_forcing_columns` /
#: :func:`validate_forcing` (the single source of truth; a test keeps this catalog
#: and those functions in sync).
FORCING_COLUMNS = (
    ForcingColumn(DATE, "ISO date",
                  "Daily timestamps; must be a continuous 1-day series. Always "
                  "required."),
    ForcingColumn(PRECIP, "mm/day", "Daily precipitation. Always required."),
    ForcingColumn(DISCHARGE, "m^3/s",
                  "Observed streamflow at the gauge (converted to mm/day using "
                  "catchment.drainage_basin_area__km2). Always required."),
    ForcingColumn(TMEAN, "deg C",
                  "Daily mean temperature; drives snowmelt and the frozen-ground "
                  "index (needed when an fdd_threshold is active or snowpack is "
                  "on). Synthesized from Min+Max if absent, so supply this column "
                  "or both Min and Max. Thornthwaite ET uses Min/Max directly, "
                  "not this column."),
    ForcingColumn(TMIN, "deg C",
                  "Daily minimum temperature; required (with Max) by "
                  "ThornthwaiteChang2019 (the Chang effective temperature uses the "
                  "diurnal range), and supplies the snowmelt mean and the "
                  "dtr_fgi_decay range."),
    ForcingColumn(TMAX, "deg C",
                  "Daily maximum temperature; see Minimum Temperature."),
    ForcingColumn(PHOTOPERIOD, "hours",
                  "Day length, from latitude and date; required by the "
                  "Thornthwaite-Chang reference ET (which also covers phenology "
                  "photoperiod senescence)."),
    ForcingColumn(ET, "mm/day",
                  "Measured/reference evapotranspiration, used directly when "
                  "evapotranspiration_method is 'datafile' (required then)."),
)

#: The config-YAML spec: top-level sections, whether each is required, and its
#: required keys. Deep per-element validation is left to Buckets.initialize().
CONFIG_SECTIONS = {
    "timeseries":   {"required": True,  "keys": ["datafile"]},
    "catchment":    {"required": True,  "keys": ["drainage_basin_area__km2",
                                                 "evapotranspiration_method",
                                                 "water_year_start_month"]},
    "general":      {"required": True,  "keys": ["spin_up_cycles"]},
    "reservoirs":   {"required": "one-of:reservoirs|sub_catchments",
                     "keys": ["recession_coefficients", "exfiltration_fractions",
                              "maximum_effective_depths__mm"]},
    "sub_catchments": {"required": "one-of:reservoirs|sub_catchments", "keys": []},
    "initial_conditions": {"required": "with-reservoirs",
                           "keys": ["water_reservoir_effective_depths__mm"]},
    "snowmelt":     {"required": False, "keys": []},
    "modules":      {"required": False, "keys": []},
    "phenology":    {"required": False, "keys": []},
}


@dataclass
class ValidationReport:
    """
    The outcome of an input-contract check.

    Attributes
    ----------
    errors : list of str
        Problems that will make :meth:`~mnished.Buckets.initialize` fail or
        produce wrong results (e.g. a missing required column).
    warnings : list of str
        Inputs MNiShed will accept but with degraded behaviour (e.g. a missing
        temperature column that silently disables snowpack).
    """

    errors: list = field(default_factory=list)
    warnings: list = field(default_factory=list)

    @property
    def ok(self):
        """``True`` when there are no errors (warnings are allowed)."""
        return not self.errors

    def extend(self, other):
        """Merge another report's errors and warnings into this one."""
        self.errors.extend(other.errors)
        self.warnings.extend(other.warnings)
        return self

    def raise_if_errors(self):
        """Raise :class:`ValueError` if there are any errors, else return self."""
        if self.errors:
            raise ValueError(str(self))
        return self

    def __str__(self):
        lines = []
        for e in self.errors:
            lines.append(f"  ERROR:   {e}")
        for w in self.warnings:
            lines.append(f"  warning: {w}")
        if not lines:
            return "MNiShed inputs OK (no contract problems found)."
        head = (f"MNiShed input contract: {len(self.errors)} error(s), "
                f"{len(self.warnings)} warning(s):")
        return "\n".join([head, *lines])


def _load_yaml(config):
    """Accept a config dict or a path to a YAML file; return (dict, base_dir)."""
    if isinstance(config, dict):
        return config, None
    path = Path(config)
    if not path.is_file():
        raise FileNotFoundError(f"config file not found: {path}")
    with open(path) as f:
        cfg = yaml.safe_load(f)
    if not isinstance(cfg, dict):
        raise ValueError(f"config {path} did not parse to a mapping")
    return cfg, path.parent


def _et_method(config):
    return (config.get("catchment", {}) or {}).get("evapotranspiration_method")


def _modules(config):
    return config.get("modules", {}) or {}


def _is_number(value):
    """True for a real numeric value (incl. inf), excluding bool/str/None."""
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _fdd_threshold_active(config):
    """True when snowmelt.fdd_threshold is a finite number (a real threshold).

    ``inf`` (the default 'never frozen') and non-numeric values are not active.
    """
    fdd = (config.get("snowmelt", {}) or {}).get("fdd_threshold")
    return _is_number(fdd) and np.isfinite(fdd)


def _has_temperature(columns):
    """The model needs a daily mean temperature, which it reads from
    'Mean Temperature [C]' or synthesizes from Min+Max — so either suffices."""
    return TMEAN in columns or (TMIN in columns and TMAX in columns)


def required_forcing_columns(config):
    """
    The forcing columns that are *hard*-required for this config.

    Returns
    -------
    dict
        ``{column_name: reason}`` for every column whose absence is an error
        given the config's ET method and options.
    """
    req = {
        DATE: "always required",
        PRECIP: "always required",
        DISCHARGE: "always required",
    }
    et = _et_method(config)
    if et == "ThornthwaiteChang2019":
        # The Chang effective temperature Tef = 0.5k(3*Tmax - Tmin) uses the
        # diurnal range, plus the photoperiod — so Min, Max, and Photoperiod are
        # the hard requirements (the snowmelt/phenology mean temperature is then
        # derived from min+max). Photoperiod also covers phenology
        # photoperiod-senescence, which only acts under Thornthwaite.
        req[PHOTOPERIOD] = "ThornthwaiteChang2019 ET (uses day length)"
        req[TMIN] = "ThornthwaiteChang2019 ET (effective temp uses the diurnal range)"
        req[TMAX] = "ThornthwaiteChang2019 ET (effective temp uses the diurnal range)"
    elif et == "datafile":
        req[ET] = "evapotranspiration_method is datafile"
    return req


def recommended_forcing_columns(config, already_required=None):
    """
    Forcing columns that are *soft*-required (missing → degraded behaviour).

    Returns ``{column_name: reason}`` for columns MNiShed will run without, but
    with a fallback. (The temperature requirement — Mean, or Min+Max — is checked
    separately by :func:`validate_forcing` since it is a disjunction.)
    """
    rec = {}
    if _modules(config).get("dtr_fgi_decay", True):
        for col in (TMIN, TMAX):
            rec[col] = ("modules.dtr_fgi_decay is on; falls back to a constant "
                        "frozen-ground-index decay without it")
    return rec


def validate_config(config):
    """
    Validate the config YAML against the input contract.

    Checks that the required top-level sections and keys are present, that
    ``evapotranspiration_method`` is recognised, that exactly one of
    ``reservoirs`` / ``sub_catchments`` defines the cascade, and that the
    reservoir list lengths are consistent. Deep model validation (area
    fractions, lake parameters, junctions) is left to
    :meth:`~mnished.Buckets.initialize`.

    Parameters
    ----------
    config : dict or str or pathlib.Path
        A parsed config mapping, or a path to a config YAML.

    Returns
    -------
    ValidationReport
    """
    cfg, _ = _load_yaml(config)
    report = ValidationReport()

    for section, spec in CONFIG_SECTIONS.items():
        if spec.get("required") is not True:    # one-of/with-reservoirs handled below
            continue
        block = cfg.get(section)
        if not isinstance(block, dict):
            report.errors.append(f"missing required config section '{section}:'")
            continue
        for key in spec["keys"]:
            if key not in block or block[key] is None:
                report.errors.append(f"config '{section}' is missing '{key}'")

    et = _et_method(cfg)
    if et is not None and et not in ET_METHODS:
        report.errors.append(
            f"catchment.evapotranspiration_method must be one of "
            f"{list(ET_METHODS)}; got {et!r}")

    fdd = (cfg.get("snowmelt", {}) or {}).get("fdd_threshold")
    if fdd is not None and not _is_number(fdd):
        report.errors.append(
            f"snowmelt.fdd_threshold must be a number [°C·day] (or omitted); "
            f"got {fdd!r}")

    has_res = isinstance(cfg.get("reservoirs"), dict)
    subs = cfg.get("sub_catchments")
    has_subs = isinstance(subs, list) and len(subs) > 0
    if has_res == has_subs:
        report.errors.append(
            "config must define the cascade with exactly one of 'reservoirs:' "
            "(single cascade) or a non-empty 'sub_catchments:' list "
            f"(got reservoirs={has_res}, sub_catchments={has_subs})")
    if has_res:
        report.extend(_validate_reservoir_block(cfg["reservoirs"], "reservoirs"))
        ic = cfg.get("initial_conditions", {}) or {}
        depths = ic.get("water_reservoir_effective_depths__mm")
        if depths is None:
            report.errors.append(
                "config with a top-level 'reservoirs:' needs "
                "'initial_conditions.water_reservoir_effective_depths__mm'")
        elif isinstance(depths, list):
            n = len(cfg["reservoirs"].get("recession_coefficients", []) or [])
            if n and len(depths) != n:
                report.errors.append(
                    f"initial_conditions.water_reservoir_effective_depths__mm has "
                    f"{len(depths)} entries but the cascade has {n} reservoirs")
    return report


def _validate_reservoir_block(block, where):
    report = ValidationReport()
    required = ["recession_coefficients", "exfiltration_fractions",
                "maximum_effective_depths__mm"]
    lengths = {}
    for key in required:
        val = block.get(key)
        if val is None:
            report.errors.append(f"'{where}' is missing required '{key}'")
        elif not isinstance(val, list):
            report.errors.append(f"'{where}.{key}' must be a list")
        else:
            lengths[key] = len(val)
    if len(set(lengths.values())) > 1:
        report.errors.append(
            f"'{where}' reservoir lists must all have the same length; got "
            + ", ".join(f"{k}={v}" for k, v in lengths.items()))
    return report


def validate_forcing(forcing, config=None):
    """
    Validate a forcing table against the input contract.

    Parameters
    ----------
    forcing : pandas.DataFrame or str or pathlib.Path
        The forcing data, or a path to the CSV.
    config : dict or str or pathlib.Path, optional
        The config the forcing will be used with. Without it, only the
        unconditional columns (Date, Precipitation, Discharge) and the
        daily-continuity check are enforced; with it, the conditionally required
        columns (ET method, frozen ground, phenology) are checked too.

    Returns
    -------
    ValidationReport
    """
    report = ValidationReport()
    if isinstance(forcing, pd.DataFrame):
        df = forcing
    else:
        path = Path(forcing)
        if not path.is_file():
            report.errors.append(f"forcing file not found: {path}")
            return report
        try:
            df = pd.read_csv(path)
        except Exception as exc:                       # malformed CSV
            report.errors.append(f"could not read forcing CSV {path}: {exc}")
            return report

    cfg = {}
    if config is not None:
        cfg, _ = _load_yaml(config)

    required = required_forcing_columns(cfg)
    for col, reason in required.items():
        if col not in df.columns:
            report.errors.append(f"forcing is missing column '{col}' ({reason})")
    for col, reason in recommended_forcing_columns(cfg, required).items():
        if col not in df.columns:
            report.warnings.append(f"forcing is missing column '{col}' ({reason})")

    # Snowmelt and frozen ground read a daily 'Mean Temperature [C]', which the
    # model synthesizes from Min+Max when absent — so either the mean column or
    # both min and max satisfies them. (ThornthwaiteChang2019 needs Min+Max in
    # their own right, already required above, so this only adds anything for the
    # datafile-ET case.) An error when an active fdd_threshold needs it; a warning
    # when only snowpack wants it (it disables silently without temperature).
    if (config is not None
            and _et_method(cfg) != "ThornthwaiteChang2019"
            and not _has_temperature(df.columns)):
        msg = ("forcing has no temperature input: needs 'Mean Temperature [C]', "
               "or both 'Minimum Temperature [C]' and 'Maximum Temperature [C]' "
               "to derive it")
        if _fdd_threshold_active(cfg):
            report.errors.append(
                f"{msg} — required by snowmelt.fdd_threshold (frozen ground)")
        elif _modules(cfg).get("snowpack", True):
            report.warnings.append(f"{msg}; snowpack is disabled without it")

    # Date continuity (mirrors Buckets.initialize's daily-series requirement).
    if DATE in df.columns:
        try:
            dates = pd.to_datetime(df[DATE])
        except (ValueError, TypeError):
            report.errors.append(f"'{DATE}' column is not parseable as dates")
        else:
            deltas = dates.diff().dropna()
            if len(deltas) and not (deltas == pd.Timedelta(days=1)).all():
                report.errors.append(
                    "forcing must be a continuous daily series (exactly 1-day "
                    "intervals); found gaps or non-daily spacing")

    # All-NaN required numeric columns are as bad as missing ones.
    for col in required:
        if col in df.columns and col != DATE and df[col].isna().all():
            report.errors.append(f"forcing column '{col}' is entirely empty (all NaN)")
    return report


def validate_inputs(config):
    """
    Validate a config *and* its forcing CSV together (the headline check).

    Loads the config, resolves ``timeseries.datafile`` relative to the config
    file, loads the forcing, and runs :func:`validate_config` and
    :func:`validate_forcing` (config-aware), returning one combined report.

    Parameters
    ----------
    config : dict or str or pathlib.Path
        A config mapping or a path to a config YAML. A path is needed to resolve
        a relative ``datafile``.

    Returns
    -------
    ValidationReport
    """
    cfg, base_dir = _load_yaml(config)
    report = validate_config(cfg)

    datafile = (cfg.get("timeseries", {}) or {}).get("datafile")
    if datafile:
        path = Path(datafile)
        if not path.is_absolute() and base_dir is not None:
            path = base_dir / path
        report.extend(validate_forcing(path, cfg))
    return report
