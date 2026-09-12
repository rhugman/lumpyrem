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

__all__ = ["Results", "ConvergenceWarning", "COLUMNS"]


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
        resid = np.abs(self.values[:, _BALANCE])
        if resid.size == 0:
            return 0.0
        worst = float(resid.max())
        if not relative:
            return worst
        inflow = float(np.abs(self.column("rainfall")).sum()
                       + np.abs(self.column("irrigation")).sum())
        return worst / inflow if inflow > 0.0 else 0.0

    # ------------------------------------------------------------------
    def to_csv(self, path, **kwargs):
        """Write the table as CSV.  Extra arguments go to pandas."""
        return self.to_dataframe().to_csv(path, **kwargs)

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
