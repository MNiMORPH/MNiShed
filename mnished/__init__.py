from ._version import __version__
from .mnished import Reservoir, Snowpack, Buckets
from .calibration import CalibResult, run_and_score
from .hydrograph_separation import HydrographSeparation
from .recession import BrutsaertNieber
from .priors import suggest_priors, Priors

__all__ = [
    "Reservoir",
    "Snowpack",
    "Buckets",
    "CalibResult",
    "run_and_score",
    "HydrographSeparation",
    "BrutsaertNieber",
    "suggest_priors",
    "Priors",
    "__version__",
]

# BmiMNiShed requires the optional `bmipy` dependency; expose it only when
# available (pip install 'mnished[bmi]').
try:
    from .bmi import BmiMNiShed  # noqa: F401  (re-export; recorded in __all__)
    __all__.append("BmiMNiShed")
except ImportError:
    pass
