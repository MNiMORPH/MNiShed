from ._version import __version__
from .calibration import CalibResult, run_and_score
from .hydrograph_separation import HydrographSeparation
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
    "HydrographSeparation",
    "BrutsaertNieber",
    "suggest_priors",
    "Priors",
    "__version__",
]

# BmiMNiShed requires the optional `bmipy` dependency; expose it only when
# available (pip install 'mnished[bmi]').
try:
    from .bmi import BmiMNiShed  # noqa: F401  (optional re-export; see __all__)
    __all__.append("BmiMNiShed")
except ImportError:
    pass
