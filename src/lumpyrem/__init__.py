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
``lumpyrem.engine``, ``lumpyrem.compiled``
    Many cells at once, on the Numba-compiled kernel when Numba is installed
    and on ``core`` when it is not.  ``ModelGrid`` is the way in.
``lumpyrem.io``
    xarray Datasets and netCDF, for gridded forcing and results.  Needs
    xarray, and is imported only when used.

Many independent cells::

    from lumpyrem import ModelGrid

    grid = ModelGrid.from_arrays(upper=dict(maxvol=maxvol, ks=ks, m=0.4, l=0.5))
    results = grid.run(forcing, times="MS")
    results.to_xarray()["total_rech"]          # (time, cell)
"""

from .forcing import Fill, Forcing, ForcingError, GridForcing
from .model import Model, ModelError, ModelGrid
from .parameters import (
    InitialState,
    LowerStore,
    ParameterError,
    Solver,
    UpperStore,
    VolumeToElevation,
)
from .results import ConvergenceWarning, GridResults, Results

__version__ = "0.0.2"

__all__ = [
    "ConvergenceWarning",
    "Fill",
    "Forcing",
    "ForcingError",
    "GridForcing",
    "GridResults",
    "InitialState",
    "LowerStore",
    "Model",
    "ModelError",
    "ModelGrid",
    "ParameterError",
    "Results",
    "Solver",
    "UpperStore",
    "VolumeToElevation",
    "__version__",
]
