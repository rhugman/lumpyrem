"""xarray and netCDF: gridded results out, gridded forcing in.

Everything here needs xarray (``pip install lumpyrem[xarray]``), imported only
when a function is called, so the rest of the package does not depend on it.
Reading or writing a file also needs a netCDF backend -- netCDF4, h5netcdf or
scipy -- which xarray chooses.

The layout is one data variable per quantity on a ``time`` axis and, for a
grid, a ``cell`` axis: ``ds["total_rech"]`` is a ``(time, cell)`` array.
Results keep their own conventions (``lumpyrem.results``): ``time`` labels the
*end* of each reporting interval, flux variables are totals over the interval
ending there, and the first row is the initial state.

Plain CSV and DataFrame input and output stay on the objects themselves, as
:meth:`Forcing.from_csv <lumpyrem.forcing.Forcing.from_csv>` and
:meth:`Results.to_csv <lumpyrem.results.Results.to_csv>`.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import pandas as pd

from . import __version__
from .forcing import VARIABLES, Fill, Forcing, ForcingError, GridForcing
from .results import COLUMNS, GridResults, Results

__all__ = ["forcing_from_dataset", "read_forcing", "to_dataset", "to_netcdf"]

#: Columns that are states at their timestamp; every other column is a total
#: over the interval ending there.
STATE_COLUMNS: tuple[str, ...] = (
    "vol_upper", "vol_lower", "vol_drain", "vol_macro", "elevation", "depth_to_water",
)


def _xarray():
    try:
        import xarray
    except ImportError:
        raise ImportError(
            "lumpyrem.io needs xarray; install it with `pip install lumpyrem[xarray]`"
        ) from None
    return xarray


# ---------------------------------------------------------------------------
# Results out
# ---------------------------------------------------------------------------

def to_dataset(results: GridResults | Results):
    """Results as an xarray Dataset, one variable per column.

    A :class:`~lumpyrem.results.GridResults` gives variables on
    ``(time, cell)`` with per-cell convergence counts alongside; a single
    :class:`~lumpyrem.results.Results` gives variables on ``time`` alone.
    """
    xr = _xarray()
    grid = isinstance(results, GridResults)
    if not grid and not isinstance(results, Results):
        raise TypeError(f"expected Results or GridResults, got {type(results).__name__}")
    dims = ("time", "cell") if grid else ("time",)

    data_vars = {name: (dims, results.column(name), _column_attrs(name))
                 for name in COLUMNS}
    coords = {"time": results.times, "day": ("time", results.days)}
    attrs = {
        "source": f"lumpyrem {__version__}",
        "time_convention": "time labels the end of each reporting interval; the "
                           "first row is the initial state",
    }
    if grid:
        coords["cell"] = list(results.cells)
        data_vars["nonconverged_upper"] = ("cell", results.nonconverged[:, 0])
        data_vars["nonconverged_lower"] = ("cell", results.nonconverged[:, 1])
        data_vars["has_elevation"] = ("cell", results.has_elevation.astype(bool))
    else:
        attrs["nonconverged_upper"] = results.nonconverged_upper
        attrs["nonconverged_lower"] = results.nonconverged_lower
        attrs["has_elevation"] = int(results.has_elevation)
    return xr.Dataset(data_vars=data_vars, coords=coords, attrs=attrs)


def to_netcdf(results: GridResults | Results, path, **kwargs):
    """Write :func:`to_dataset` to a netCDF file.  Extra arguments go to xarray."""
    return to_dataset(results).to_netcdf(path, **kwargs)


def _column_attrs(name: str) -> dict:
    if name in STATE_COLUMNS:
        return {"cell_methods": "time: point"}
    return {"cell_methods": "time: sum (comment: over the interval ending at time)"}


# ---------------------------------------------------------------------------
# Forcing in
# ---------------------------------------------------------------------------

def forcing_from_dataset(ds, *, time: str = "time", cell: str = "cell",
                         fill: Mapping[str, str | Fill] | None = None
                         ) -> GridForcing | Forcing:
    """Forcing from an xarray Dataset whose variables are named as in
    :data:`~lumpyrem.forcing.VARIABLES`.

    A variable on ``(time, cell)`` differs between cells, and one on
    ``cell`` alone is a per-cell constant; one on ``time`` alone is shared, and
    a scalar is a constant everywhere.  Any other variable is
    rejected rather than ignored, as :meth:`Forcing.from_dataframe` does.

    Returns:
        A :class:`~lumpyrem.forcing.GridForcing` labelled with the ``cell``
        coordinate when any variable has a ``cell`` axis, otherwise a
        :class:`~lumpyrem.forcing.Forcing`.
    """
    unknown = [str(k) for k in ds.data_vars if k not in VARIABLES]
    if unknown:
        raise ForcingError(
            f"unknown forcing variable(s) {unknown}; known names are {list(VARIABLES)}"
        )
    if time not in ds.coords:
        raise ForcingError(f"the dataset has no {time!r} coordinate")

    arrays, per_cell = {}, False
    ntime = ds.sizes[time]
    for name, var in ds.data_vars.items():
        extra = [d for d in var.dims if d not in (time, cell)]
        if extra:
            raise ForcingError(f"{name}: dimensions {var.dims}; expected some of "
                               f"({time!r}, {cell!r})")
        if var.dims == (cell,):
            # A per-cell constant: the same value on every day, cell by cell.
            per_cell = True
            arrays[name] = np.broadcast_to(var.to_numpy(), (ntime, var.size)).copy()
        elif cell in var.dims:
            per_cell = True
            arrays[name] = var.transpose(time, cell).to_numpy()
        elif var.dims:
            arrays[name] = var.to_numpy()
        else:
            arrays[name] = var.item()

    dates = pd.DatetimeIndex(ds[time].to_numpy())
    if per_cell:
        labels = tuple(ds[cell].to_numpy().tolist()) if cell in ds.coords else None
        return GridForcing.from_arrays(dates, fill=fill, cells=labels, **arrays)
    series = {k: v for k, v in arrays.items() if np.ndim(v) == 1}
    constants = {k: v for k, v in arrays.items() if np.ndim(v) == 0}
    return Forcing.from_dataframe(pd.DataFrame(series, index=dates), fill=fill,
                                  **constants)


def read_forcing(path, *, time: str = "time", cell: str = "cell",
                 fill: Mapping[str, str | Fill] | None = None, **open_kwargs):
    """Read forcing from a netCDF file; see :func:`forcing_from_dataset`."""
    xr = _xarray()
    with xr.open_dataset(path, **open_kwargs) as ds:
        return forcing_from_dataset(ds.load(), time=time, cell=cell, fill=fill)
