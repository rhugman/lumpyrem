"""Model output: a labelled table, a balance check, a plot and a writer.

The flux columns are *interval totals*, not rates.  Row ``k`` covers the days
after the previous output time up to and including its own, so a monthly output
schedule gives monthly totals and a daily one gives daily totals.  The volume
columns are instantaneous, at the end of the interval.  Row 0 is the initial
state: its volumes are the initial conditions and every flux is zero.

Because the index labels interval *ends*, ``results.df.loc[a:b]`` reads as
"everything that happened after ``a``, up to ``b``", which is also how a
MODFLOW stress period reads.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .core import COLUMNS as KERNEL_COLUMNS

__all__ = ["Results", "GridResults", "ConvergenceWarning", "COLUMNS"]


class ConvergenceWarning(UserWarning):
    """A sub-step exhausted its Picard iterations.

    The Fortran prints this to stdout and carries on with whatever iterate it
    was holding, which is how a non-converged answer gets read as a converged
    one.  Results are still returned -- they are just not solutions.
    """


def _rename(name: str) -> str:
    return "depth_to_water" if name == "depth-to-water" else name


#: Output columns.  The kernel's names, with the one hyphen removed.
COLUMNS: tuple[str, ...] = tuple(_rename(c) for c in KERNEL_COLUMNS)

_BALANCE = COLUMNS.index("balance")
_RAINFALL = COLUMNS.index("rainfall")
_IRRIGATION = COLUMNS.index("irrigation")
_DEFAULT_PLOT = ("rainfall", "total_rech", "evap_upper", "runoff")


@dataclass(frozen=True)
class Results:
    """A completed run.

    Attributes:
        times: output instants, length ``nout + 1``; ``times[0]`` is the start
            of the run and carries the initial state.
        days: the same instants as 1-based simulation day counts, ``days[0]``
            being 0.
        values: ``(nout + 1, 26)`` float64, in :data:`COLUMNS` order.
        nonconverged_upper: sub-steps in which the upper store hit ``mxiter``.
        nonconverged_lower: the same for the lower store.
        has_elevation: False when the model carried no volume-to-elevation
            conversion, in which case the last two columns are NaN.
    """

    times: pd.DatetimeIndex
    days: np.ndarray
    values: np.ndarray
    nonconverged_upper: int = 0
    nonconverged_lower: int = 0
    has_elevation: bool = True
    columns: tuple[str, ...] = COLUMNS

    # ------------------------------------------------------------------
    @property
    def converged(self) -> bool:
        """True when every sub-step of the run met the solver tolerance."""
        return self.nonconverged_steps == 0

    @property
    def nonconverged_steps(self) -> int:
        return self.nonconverged_upper + self.nonconverged_lower

    def __len__(self) -> int:
        return len(self.days)

    def __getitem__(self, name: str) -> pd.Series:
        return pd.Series(self.column(name), index=self.times, name=name)

    def column(self, name: str) -> np.ndarray:
        """One column as a float64 array."""
        try:
            return self.values[:, self.columns.index(name)]
        except ValueError:
            raise KeyError(
                f"no column {name!r}; available columns are {list(self.columns)}"
            ) from None

    def to_dataframe(self) -> pd.DataFrame:
        """The whole table, indexed by output time, with ``day`` alongside."""
        frame = pd.DataFrame(self.values, index=self.times, columns=list(self.columns))
        frame.insert(0, "day", self.days)
        frame.index.name = "time"
        return frame

    @property
    def df(self) -> pd.DataFrame:
        """Shorthand for :meth:`to_dataframe`."""
        return self.to_dataframe()

    # ------------------------------------------------------------------
    def balance_error(self, *, relative: bool = False) -> float:
        """Worst water-balance residual over the run.

        The kernel computes the residual itself, per interval, in the
        ``balance`` column: inflow minus outflow minus the change in every
        store and buffer.  A well-posed run closes to rounding.

        Two known defects in the reference can open it legitimately -- the
        store-empties fix-up being bypassed by irrigation, and buffer water
        claimed on a first pass never reaching the lower store.  Both are
        reproduced here on purpose, so a non-zero result is not automatically
        a bug in the caller's setup; see ``docs/conversion-plan.md``, traps 13
        and 14.

        Args:
            relative: divide by total inflow (rainfall plus irrigation) instead
                of returning the residual in depth units.

        Returns:
            The largest absolute residual, or 0.0 for an empty run.
        """
        return _balance_error(self.values, relative)

    # ------------------------------------------------------------------
    def to_csv(self, path, **kwargs):
        """Write the table as CSV.  Extra arguments go to pandas."""
        return self.to_dataframe().to_csv(path, **kwargs)

    def to_xarray(self):
        """An xarray Dataset on ``time``, one variable per column.  Needs xarray."""
        from .io import to_dataset
        return to_dataset(self)

    def plot(self, columns=None, *, ax=None, **kwargs):
        """Plot columns against time and return the Axes.

        Args:
            columns: column names to draw; a useful flux set by default.
            ax: existing Axes to draw on.
            **kwargs: passed to :meth:`pandas.DataFrame.plot`.
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError(
                "Results.plot needs matplotlib; install it with "
                "`pip install lumpyrem[plot]` or `pip install matplotlib`"
            ) from None
        names = list(_DEFAULT_PLOT if columns is None else columns)
        for name in names:
            self.column(name)      # raises KeyError before anything is drawn
        if ax is None:
            _, ax = plt.subplots()
        frame = self.to_dataframe()[names]
        frame.plot(ax=ax, **kwargs)
        ax.set_xlabel("")
        ax.set_ylabel("depth per interval")
        return ax

    def __repr__(self) -> str:
        if len(self.times) == 0:
            return "<Results: empty>"
        state = "converged" if self.converged else (
            f"NOT converged ({self.nonconverged_steps} sub-steps)")
        return (f"<Results: {len(self.days) - 1} output times, "
                f"{self.times[0].date()} to {self.times[-1].date()}, {state}>")


def _balance_error(values: np.ndarray, relative: bool) -> float:
    """Worst ``balance`` residual of one cell's ``(nrows, 26)`` table."""
    resid = np.abs(values[:, _BALANCE])
    if resid.size == 0:
        return 0.0
    worst = float(resid.max())
    if not relative:
        return worst
    inflow = float(np.abs(values[:, _RAINFALL]).sum() + np.abs(values[:, _IRRIGATION]).sum())
    return worst / inflow if inflow > 0.0 else 0.0


# ---------------------------------------------------------------------------
# Many cells
# ---------------------------------------------------------------------------

@dataclass(frozen=True, eq=False)
class GridResults:
    """A completed :class:`~lumpyrem.ModelGrid` run: every cell on one schedule.

    Columns and conventions are those of :class:`Results`, with a ``cell``
    axis added: ``values[t, c, k]`` is column ``k`` of cell ``c`` at output
    time ``t``.  ``results["total_rech"]`` is then a ``(time, cell)`` table,
    and :meth:`to_xarray` gives the same as a Dataset -- one variable per
    column on ``(time, cell)`` -- which drops a recharge array straight into
    a gridded groundwater model.

    Attributes:
        times: output instants, ``times[0]`` the start of the run.
        days: the same instants as simulation day counts, ``days[0] == 0``.
        values: ``(ntime, ncell, 26)`` float64, in :data:`COLUMNS` order.
        cells: the cell labels, in order.
        nonconverged: ``(ncell, 2)`` sub-step counts, upper then lower store,
            that hit ``mxiter``.
        has_elevation: ``(ncell,)``; False where a cell carried no
            volume-to-elevation conversion and its last two columns are NaN.
    """

    times: pd.DatetimeIndex
    days: np.ndarray
    values: np.ndarray
    cells: tuple
    nonconverged: np.ndarray
    has_elevation: np.ndarray
    columns: tuple[str, ...] = COLUMNS

    # ------------------------------------------------------------------
    @property
    def ncell(self) -> int:
        return len(self.cells)

    @property
    def converged(self) -> bool:
        """True when every sub-step of every cell met the solver tolerance."""
        return not self.nonconverged.any()

    @property
    def nonconverged_cells(self) -> list:
        """Labels of the cells with any sub-step that hit ``mxiter``."""
        return [self.cells[i] for i in np.flatnonzero(self.nonconverged.sum(axis=1))]

    def __len__(self) -> int:
        return len(self.days)

    def column(self, name: str) -> np.ndarray:
        """One column as a ``(time, cell)`` float64 array."""
        try:
            return self.values[:, :, self.columns.index(name)]
        except ValueError:
            raise KeyError(
                f"no column {name!r}; available columns are {list(self.columns)}"
            ) from None

    def __getitem__(self, name: str) -> pd.DataFrame:
        """One column as a table, indexed by time with a column per cell."""
        frame = pd.DataFrame(self.column(name), index=self.times,
                             columns=pd.Index(self.cells, name="cell"))
        frame.index.name = "time"
        return frame

    def cell(self, label) -> Results:
        """One cell's results, exactly as :meth:`Model.run` would return them."""
        try:
            i = self.cells.index(label)
        except ValueError:
            raise KeyError(f"no cell {label!r}") from None
        return Results(
            times=self.times, days=self.days, values=self.values[:, i, :].copy(),
            nonconverged_upper=int(self.nonconverged[i, 0]),
            nonconverged_lower=int(self.nonconverged[i, 1]),
            has_elevation=bool(self.has_elevation[i]),
        )

    def balance_error(self, *, relative: bool = False) -> pd.Series:
        """Each cell's worst water-balance residual; see :meth:`Results.balance_error`."""
        return pd.Series([_balance_error(self.values[:, i, :], relative)
                          for i in range(self.ncell)],
                         index=pd.Index(self.cells, name="cell"), name="balance_error")

    # ------------------------------------------------------------------
    def to_xarray(self, *, cell_levels=None):
        """An xarray Dataset on ``(time, cell)``, one variable per column.

        Tuple cell labels become a MultiIndex with levels named by
        ``cell_levels``; see :func:`lumpyrem.io.to_dataset`.
        """
        from .io import to_dataset
        return to_dataset(self, cell_levels=cell_levels)

    def to_netcdf(self, path, *, cell_levels=None, **kwargs):
        """Write :meth:`to_xarray` to netCDF.  Extra arguments go to xarray."""
        from .io import to_netcdf
        return to_netcdf(self, path, cell_levels=cell_levels, **kwargs)

    def __repr__(self) -> str:
        failed = len(self.nonconverged_cells)
        state = "converged" if failed == 0 else f"{failed} cell(s) NOT converged"
        return (f"<GridResults: {self.ncell} cells x {len(self.days) - 1} output times, "
                f"{self.times[0].date()} to {self.times[-1].date()}, {state}>")
