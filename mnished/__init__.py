from ._version import __version__
from .calibration import (
    CalibResult,
    Calibrator,
    ScoringModel,
    log_flow_residual_terms,
    run_and_score,
    target_kwargs,
)
from .diagnostics import SeasonalMassBalance
from .hydrograph_separation import HydrographSeparation
from .identifiability import (
    IdentifiabilityReport,
    Parameter,
    ParameterSet,
    eigenspectrum,
    profile,
    profile_all,
    ridge,
)
from .mnished import Buckets, Reservoir, Snowpack, SubCatchment
from .priors import Priors, suggest_priors
from .recession import BrutsaertNieber

__all__ = [
    "Reservoir",
    "Snowpack",
    "SubCatchment",
    "Buckets",
    "CalibResult",
    "run_and_score",
    "ScoringModel",
    "Calibrator",
    "target_kwargs",
    "log_flow_residual_terms",
    "SeasonalMassBalance",
    "HydrographSeparation",
    "BrutsaertNieber",
    "suggest_priors",
    "Priors",
    "Parameter",
    "ParameterSet",
    "profile",
    "profile_all",
    "eigenspectrum",
    "ridge",
    "IdentifiabilityReport",
    "__version__",
]

# BmiMNiShed requires the optional `bmipy` dependency; expose it only when
# available (pip install 'mnished[bmi]').
try:
    from .bmi import BmiMNiShed  # noqa: F401  (optional re-export; see __all__)
    __all__.append("BmiMNiShed")
except ImportError:
    pass
