"""lumpyrem -- lumped-parameter recharge modelling.

A Python port of the LUMPREM/LUMPREM2 Fortran models, adding array-based
forcing and output, a modern API, and compiled-speed execution.

Quick start::

    import pandas as pd
    from lumpyrem import Model, UpperStore, Forcing

    forcing = Forcing.from_csv("climate.csv")
    model = Model(upper=UpperStore(maxvol=0.3, ks=0.05, m=0.4, l=0.5))
    results = model.run(forcing, times="MS")
    results.df["total_rech"]

Layers, each usable on its own:

``lumpyrem.core``
    The kernel: a literal translation of ``rechmod2.f``, held to bit-fidelity
    against the Fortran.  Floats in, arrays out.
``lumpyrem.parameters``, ``lumpyrem.forcing``, ``lumpyrem.model``
    The designed API -- typed parameters, real dates, explicit resampling.
``lumpyrem.results``
    Labelled output, a balance check, a plot and a CSV writer.

Phases 0-2 of ``docs/conversion-plan.md`` are complete.  Array-based and
compiled execution arrive in Phase 3.
"""

from .forcing import Fill, Forcing, ForcingError
from .model import Model, ModelError
from .parameters import (
    InitialState,
    LowerStore,
    ParameterError,
    Solver,
    UpperStore,
    VolumeToElevation,
)
from .results import ConvergenceWarning, Results

__version__ = "0.0.2"

__all__ = [
    "ConvergenceWarning",
    "Fill",
    "Forcing",
    "ForcingError",
    "InitialState",
    "LowerStore",
    "Model",
    "ModelError",
    "ParameterError",
    "Results",
    "Solver",
    "UpperStore",
    "VolumeToElevation",
    "__version__",
]
