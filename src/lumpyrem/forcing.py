"""Forcing input, and the resampling rules that turn it into dense daily arrays.

``lumprem2.f`` gap-fills its three forcing files by three *different* rules,
decided by which file a value happened to be written in: vegetation is linearly
interpolated between listed days, potential evaporation and irrigation are
forward-filled as steps, and rainfall is zero-filled.  The rules are
hydrologically sensible -- rain really does not persist, and a crop factor
really does move smoothly -- but in the Fortran they are a property of the file
format, invisible and unchangeable.

Here they are named, defaulted to the Fortran's choices, and per-variable
overridable.  :data:`FORTRAN_FILL` is the default map; pass ``fill=`` to change
any of it.

Applying one rule to all three variables is the easiest way to get plausible
but wrong answers out of this model, which is why the rule is visible.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

import numpy as np
import pandas as pd

__all__ = ["Fill", "Forcing", "GridForcing", "VARIABLES", "FORTRAN_FILL", "ForcingError"]


class Fill(str, Enum):
    """How a variable's listed values extend to the days between them."""

    #: Linear interpolation between listed days.  The Fortran's vegetation rule.
    INTERPOLATE = "interpolate"
    #: Hold the last listed value until the next one.  The epot/irrigation rule.
    FORWARD = "forward"
    #: Anything not listed is zero.  The rainfall rule.
    ZERO = "zero"
    #: Refuse to fill: every day in the run must carry a listed value.
    DENSE = "dense"


#: The forcing variables, in the order the kernel takes them.
VARIABLES: tuple[str, ...] = (
    "rainfall",
    "pot_evap",
    "crop_factor",
    "veg_gamma",
    "irrigate",
    "gw_irrig_frac",
    "pot_evap_lower",
)

#: What ``lumprem2.f`` does, per variable.  The default.
FORTRAN_FILL: dict[str, Fill] = {
    "rainfall": Fill.ZERO,
    "pot_evap": Fill.FORWARD,
    "crop_factor": Fill.INTERPOLATE,
    "veg_gamma": Fill.INTERPOLATE,
    "irrigate": Fill.FORWARD,
    "gw_irrig_frac": Fill.FORWARD,
    "pot_evap_lower": Fill.FORWARD,
}

#: Used when a variable is not supplied at all.
DEFAULTS: dict[str, float] = {
    "crop_factor": 1.0,
    "irrigate": 0.0,
    "gw_irrig_frac": 0.0,
    "pot_evap_lower": 0.0,
}

_REQUIRED = ("rainfall", "pot_evap", "veg_gamma")

_DAY = pd.Timedelta(days=1)


class ForcingError(ValueError):
    """Forcing could not be resolved onto the requested daily grid."""


def _day_index(times) -> pd.DatetimeIndex:
    """Coerce to a DatetimeIndex and insist every stamp is a whole day."""
    idx = pd.DatetimeIndex(pd.to_datetime(list(times) if not isinstance(
        times, (pd.DatetimeIndex, pd.Series)) else times))
    if idx.tz is not None:
        raise ForcingError("forcing times must be timezone-naive; "
                           "the model's step is a calendar day, not an instant")
    off = idx != idx.normalize()
    if off.any():
        raise ForcingError(
            f"forcing times must fall on midnight -- the model's step is a whole "
            f"day; first offender {idx[off][0]}"
        )
    return idx


# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Forcing:
    """Time series driving the model, at whatever frequency they were observed.

    Build one with :meth:`from_dataframe`, :meth:`from_arrays`,
    :meth:`from_csv` or :meth:`from_dict` rather than calling the constructor.

    Variables, all in depth units per day where they are fluxes:

    ============== ============================================ ==============
    name           meaning                                      default fill
    ============== ============================================ ==============
    rainfall       rainfall rate                                zero
    pot_evap       potential evaporation of the upper store     forward
    crop_factor    multiplies ``pot_evap``                      interpolate
    veg_gamma      evaporation curvature, > 0                   interpolate
    irrigate       1 on days irrigation may top the store up    forward
    gw_irrig_frac  fraction of irrigation drawn from groundwater forward
    pot_evap_lower potential evaporation of the lower store     forward
    ============== ============================================ ==============

    ``rainfall``, ``pot_evap`` and ``veg_gamma`` must be supplied; the rest
    default to a constant.  A value may be a series or a plain number, and a
    number is constant over the whole run regardless of its fill rule.
    """

    frame: pd.DataFrame
    constants: dict[str, float]
    fill: dict[str, Fill]

    # ------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------
    @classmethod
    def from_dataframe(cls, frame: pd.DataFrame, *,
                       fill: Mapping[str, str | Fill] | None = None,
                       **constants: float) -> Forcing:
        """Build from a DataFrame indexed by date.

        Columns are matched to :data:`VARIABLES` by name; any other column is
        rejected rather than ignored, because a misspelled column that is
        silently dropped is a forcing series that silently does nothing.  A
        ``NaN`` means "not listed on this day", and is filled by that
        variable's rule exactly as a missing day would be.

        Args:
            frame: the series, indexed by date.
            fill: per-variable overrides of :data:`FORTRAN_FILL`.
            **constants: variables supplied as a single number instead.
        """
        if not isinstance(frame, pd.DataFrame):
            raise ForcingError(f"expected a DataFrame, got {type(frame).__name__}")
        unknown = [c for c in frame.columns if c not in VARIABLES]
        if unknown:
            raise ForcingError(
                f"unknown forcing column(s) {unknown}; known names are {list(VARIABLES)}"
            )
        clash = sorted(set(constants) & set(frame.columns))
        if clash:
            raise ForcingError(f"{clash} supplied both as a column and as a constant")

        idx = _day_index(frame.index)
        frame = pd.DataFrame(
            {c: np.asarray(frame[c], dtype=np.float64) for c in frame.columns},
            index=idx,
        ).sort_index()
        if frame.index.has_duplicates:
            dup = frame.index[frame.index.duplicated()][0]
            raise ForcingError(f"duplicate forcing date {dup}")

        return cls(frame=frame,
                   constants=_check_constants(constants),
                   fill=_resolve_fill(fill))

    @classmethod
    def from_arrays(cls, start, *, freq: str = "D",
                    fill: Mapping[str, str | Fill] | None = None,
                    **series) -> Forcing:
        """Build from equal-length arrays starting at ``start``.

        Args:
            start: date of the first value.
            freq: spacing of the values; daily by default.
            fill: per-variable overrides of :data:`FORTRAN_FILL`.
            **series: arrays, or plain numbers for constants.
        """
        arrays, constants = {}, {}
        for name, value in series.items():
            if np.ndim(value) == 0:
                constants[name] = value
            else:
                arrays[name] = np.asarray(value, dtype=np.float64)
        lengths = {len(v) for v in arrays.values()}
        if len(lengths) > 1:
            raise ForcingError(
                "arrays have different lengths: "
                + ", ".join(f"{k}={len(v)}" for k, v in sorted(arrays.items()))
            )
        if not arrays:
            raise ForcingError("from_arrays needs at least one array; "
                               "for an all-constant run give rainfall=[...] explicitly")
        n = lengths.pop()
        index = pd.date_range(pd.Timestamp(start), periods=n, freq=freq)
        return cls.from_dataframe(pd.DataFrame(arrays, index=index),
                                  fill=fill, **constants)

    @classmethod
    def from_dict(cls, data: Mapping[str, object], *, index=None,
                  fill: Mapping[str, str | Fill] | None = None) -> Forcing:
        """Build from a mapping of variable name to series or number.

        ``index`` supplies the dates; if omitted, the mapping must carry them
        under a ``"time"``, ``"date"`` or ``"datetime"`` key.
        """
        data = dict(data)
        if index is None:
            for key in ("time", "date", "datetime"):
                if key in data:
                    index = data.pop(key)
                    break
            else:
                raise ForcingError(
                    "no dates: pass index=..., or include a 'time', 'date' or "
                    "'datetime' entry in the mapping"
                )
        arrays = {k: v for k, v in data.items() if np.ndim(v) != 0}
        constants = {k: v for k, v in data.items() if np.ndim(v) == 0}
        frame = pd.DataFrame(
            {k: np.asarray(v, dtype=np.float64) for k, v in arrays.items()},
            index=_day_index(index),
        )
        return cls.from_dataframe(frame, fill=fill, **constants)

    @classmethod
    def from_csv(cls, path, *, time_column: str | int = 0,
                 fill: Mapping[str, str | Fill] | None = None,
                 **read_csv_kwargs) -> Forcing:
        """Read a CSV whose first column is the date and whose headers name
        forcing variables.

        Args:
            path: file to read.
            time_column: name or position of the date column.
            fill: per-variable overrides of :data:`FORTRAN_FILL`.
            **read_csv_kwargs: passed straight to :func:`pandas.read_csv`.
        """
        frame = pd.read_csv(Path(path), **read_csv_kwargs)
        col = frame.columns[time_column] if isinstance(time_column, int) else time_column
        if col not in frame.columns:
            raise ForcingError(f"{path}: no column {col!r}; found {list(frame.columns)}")
        frame = frame.set_index(col)
        frame.index = pd.to_datetime(frame.index)
        return cls.from_dataframe(frame, fill=fill)

    @classmethod
    def coerce(cls, obj, **kwargs) -> Forcing:
        """Accept a Forcing, DataFrame, mapping or CSV path interchangeably."""
        if isinstance(obj, Forcing):
            if kwargs:
                raise ForcingError("options cannot be applied to an existing Forcing")
            return obj
        if isinstance(obj, pd.DataFrame):
            return cls.from_dataframe(obj, **kwargs)
        if isinstance(obj, Mapping):
            return cls.from_dict(obj, **kwargs)
        if isinstance(obj, (str, Path)):
            return cls.from_csv(obj, **kwargs)
        raise ForcingError(
            f"cannot read forcing from {type(obj).__name__}; pass a Forcing, a "
            f"DataFrame, a mapping or a path to a CSV"
        )

    # ------------------------------------------------------------------
    # Span and resampling
    # ------------------------------------------------------------------
    @property
    def span(self) -> tuple[pd.Timestamp, pd.Timestamp]:
        """First and last listed date."""
        if len(self.frame.index) == 0:
            raise ForcingError("forcing carries no dates")
        return self.frame.index[0], self.frame.index[-1]

    def provided(self) -> set[str]:
        """Variables supplied, as either a series or a constant."""
        return set(self.frame.columns) | set(self.constants)

    def to_daily(self, start=None, end=None, *,
                 need: Sequence[str] = _REQUIRED) -> dict[str, np.ndarray]:
        """Resolve onto one value per day over ``[start, end]``.

        Args:
            start: first simulated day; defaults to the first listed date.
            end: last simulated day; defaults to the last listed date.
            need: variables that must resolve to something.  Whatever is
                neither listed nor constant falls back to :data:`DEFAULTS`.

        Returns:
            One float64 array per variable, each of length ``(end - start).days + 1``.
        """
        days = _daily_window(self, start, end, need)
        out: dict[str, np.ndarray] = {}
        for name in VARIABLES:
            if name in self.constants:
                out[name] = np.full(len(days), float(self.constants[name]))
            elif name in self.frame.columns:
                out[name] = _expand(name, self.frame[name], days, self.fill[name])
            else:
                out[name] = np.full(len(days), DEFAULTS[name])
        _validate(out)
        return out


# ---------------------------------------------------------------------------
# Per-cell forcing
# ---------------------------------------------------------------------------

@dataclass(frozen=True, eq=False)
class GridForcing:
    """Forcing that differs from cell to cell, for a :class:`~lumpyrem.ModelGrid`.

    Each variable is one of three things: an ``(ntime, ncell)`` block, one
    column per cell; a series of length ``ntime`` shared by every cell; or a
    plain number.  Variable names, defaults and fill rules are exactly those
    of :class:`Forcing`, applied to each cell's column, and a NaN means "not
    listed on this day" as it does there.

    Build one with :meth:`from_arrays`, or read an xarray Dataset with
    :func:`lumpyrem.io.forcing_from_dataset`.  A grid whose cells all share
    their forcing takes a plain :class:`Forcing` instead.

    ``cells`` optionally labels the columns.  A grid given labelled forcing
    insists the labels match its own cell names, in order, so that forcing
    read from a file cannot be silently applied to the wrong cells.
    """

    index: pd.DatetimeIndex
    ncell: int
    blocks: dict[str, np.ndarray]
    series: dict[str, np.ndarray]
    constants: dict[str, float]
    fill: dict[str, Fill]
    cells: tuple | None = None

    @classmethod
    def from_arrays(cls, time, *, freq: str = "D",
                    fill: Mapping[str, str | Fill] | None = None,
                    cells=None, **variables) -> GridForcing:
        """Build from arrays whose first axis is time.

        Args:
            time: the date of the first row, or a sequence of one date per row.
            freq: spacing of the rows when ``time`` is a single date.
            fill: per-variable overrides of :data:`FORTRAN_FILL`.
            cells: labels for the columns, or ``None``.
            **variables: ``(ntime, ncell)`` blocks, ``(ntime,)`` shared
                series, or numbers.
        """
        unknown = [k for k in variables if k not in VARIABLES]
        if unknown:
            raise ForcingError(
                f"unknown forcing variable(s) {unknown}; known names are {list(VARIABLES)}"
            )
        blocks, series, constants = {}, {}, {}
        for name, value in variables.items():
            arr = np.asarray(value, dtype=np.float64)
            if arr.ndim == 0:
                constants[name] = value
            elif arr.ndim == 1:
                series[name] = arr
            elif arr.ndim == 2:
                blocks[name] = arr
            else:
                raise ForcingError(f"{name}: expected (ntime, ncell), (ntime,) or a "
                                   f"number; got shape {arr.shape}")
        if not blocks:
            raise ForcingError("GridForcing needs at least one (ntime, ncell) block; "
                               "forcing every cell shares is a plain Forcing")
        ncells = {b.shape[1] for b in blocks.values()}
        if len(ncells) > 1:
            raise ForcingError("blocks disagree on the number of cells: " + ", ".join(
                f"{k}={v.shape[1]}" for k, v in sorted(blocks.items())))
        lengths = {len(a) for a in (*blocks.values(), *series.values())}
        if len(lengths) > 1:
            raise ForcingError("arrays disagree on the number of times: " + ", ".join(
                f"{k}={len(v)}" for k, v in sorted({**blocks, **series}.items())))
        ntime = lengths.pop()

        if np.ndim(time) == 0:
            index = _day_index(pd.date_range(pd.Timestamp(time), periods=ntime, freq=freq))
        else:
            index = _day_index(time)
            if len(index) != ntime:
                raise ForcingError(f"{len(index)} dates for {ntime} rows")
        if index.has_duplicates:
            raise ForcingError(f"duplicate forcing date {index[index.duplicated()][0]}")
        ncell = ncells.pop()
        if cells is not None:
            cells = tuple(cells)
            if len(cells) != ncell:
                raise ForcingError(f"{len(cells)} cell labels for {ncell} columns")
        order = np.argsort(index.asi8, kind="stable")
        return cls(index=index[order], ncell=ncell, cells=cells,
                   blocks={k: np.ascontiguousarray(v[order]) for k, v in blocks.items()},
                   series={k: v[order] for k, v in series.items()},
                   constants=_check_constants(constants),
                   fill=_resolve_fill(fill))

    @property
    def span(self) -> tuple[pd.Timestamp, pd.Timestamp]:
        """First and last listed date."""
        if len(self.index) == 0:
            raise ForcingError("forcing carries no dates")
        return self.index[0], self.index[-1]

    def provided(self) -> set[str]:
        """Variables supplied, as a block, a series or a constant."""
        return set(self.blocks) | set(self.series) | set(self.constants)

    def to_daily(self, start=None, end=None, *,
                 need: Sequence[str] = _REQUIRED) -> dict[str, np.ndarray]:
        """Resolve onto one value per day over ``[start, end]``.

        Returns one float64 array per variable: ``(ndays, ncell)`` for a
        block, ``(ndays,)`` for anything every cell shares.
        """
        days = _daily_window(self, start, end, need)
        out: dict[str, np.ndarray] = {}
        for name in VARIABLES:
            if name in self.constants:
                out[name] = np.full(len(days), float(self.constants[name]))
            elif name in self.blocks:
                out[name] = _expand_block(name, self.index, self.blocks[name], days,
                                          self.fill[name])
            elif name in self.series:
                out[name] = _expand(name, pd.Series(self.series[name], index=self.index),
                                    days, self.fill[name])
            else:
                out[name] = np.full(len(days), DEFAULTS[name])
        _validate(out)
        return out


def _daily_window(forcing: Forcing | GridForcing, start, end,
                  need: Sequence[str]) -> pd.DatetimeIndex:
    """The days ``to_daily`` resolves onto, once the request is known to be
    answerable: whole days, in order, with every needed variable available."""
    first, last = forcing.span
    start = first if start is None else pd.Timestamp(start)
    end = last if end is None else pd.Timestamp(end)
    if start != start.normalize() or end != end.normalize():
        raise ForcingError("start and end must fall on midnight")
    if end < start:
        raise ForcingError(f"end {end.date()} precedes start {start.date()}")
    missing = [v for v in need if v not in forcing.provided() and v not in DEFAULTS]
    if missing:
        raise ForcingError(
            f"forcing is missing {missing}; supply each as a series or a number"
        )
    return pd.date_range(start, end, freq="D")


# ---------------------------------------------------------------------------
# Rules
# ---------------------------------------------------------------------------

def _expand(name: str, series: pd.Series, days: pd.DatetimeIndex, rule: Fill) -> np.ndarray:
    """Apply one fill rule to one series, whose NaNs mean "not listed"."""
    listed = series.dropna()
    return _fill_listed(name, listed.index, listed.to_numpy(dtype=np.float64), days, rule)


def _expand_block(name: str, index: pd.DatetimeIndex, block: np.ndarray,
                  days: pd.DatetimeIndex, rule: Fill) -> np.ndarray:
    """Apply one fill rule to each column of an ``(ntime, ncell)`` block.

    Where every cell lists the same days the rule is applied to all columns
    at once -- the same element-wise arithmetic, so each column comes out as
    it would alone.  Cells that list different days are filled one by one.
    """
    nan = np.isnan(block)
    keep = ~nan.all(axis=1)
    if not nan[keep].any():
        return _fill_listed(name, index[keep], block[keep], days, rule)
    out = np.empty((len(days), block.shape[1]))
    for c in range(block.shape[1]):
        out[:, c] = _expand(f"{name} (cell {c})", pd.Series(block[:, c], index=index),
                            days, rule)
    return out


def _fill_listed(name: str, listed: pd.DatetimeIndex, values: np.ndarray,
                 days: pd.DatetimeIndex, rule: Fill) -> np.ndarray:
    """Apply one fill rule.  Day numbers are 1-based, as in the Fortran.

    ``values`` is ``(nlisted,)``, or ``(nlisted, ncell)`` for cells listing the
    same days; the result has ``len(days)`` rows either way.
    """
    if len(listed) == 0:
        raise ForcingError(f"{name}: no values listed")

    # Day numbers relative to the run, so a value listed before the run starts
    # is a negative day and still anchors a forward fill or an interpolation.
    start = days[0]
    listed_days = ((listed - start) // _DAY).to_numpy().astype(np.float64) + 1.0
    target = np.arange(1.0, len(days) + 1.0)
    shape = (len(days),) + values.shape[1:]

    if rule is Fill.ZERO:
        out = np.zeros(shape)
        inside = (listed_days >= 1.0) & (listed_days <= len(days))
        out[listed_days[inside].astype(np.int64) - 1] = values[inside]
        return out

    if rule is Fill.DENSE:
        out = np.full(shape, np.nan)
        inside = (listed_days >= 1.0) & (listed_days <= len(days))
        out[listed_days[inside].astype(np.int64) - 1] = values[inside]
        gap = np.isnan(out).reshape(len(days), -1).any(axis=1)
        if gap.any():
            first = days[int(np.argmax(gap))]
            raise ForcingError(
                f"{name}: fill='dense' requires a value on every day, and "
                f"{first.date()} has none"
            )
        return out

    if rule is Fill.FORWARD:
        if listed_days[0] > 1.0:
            raise ForcingError(
                f"{name}: fill='forward' has nothing to carry forward onto "
                f"{days[0].date()} -- the first listed value is "
                f"{listed[0].date()}. List a value at or before the run "
                f"starts, or choose a different fill rule."
            )
        # searchsorted with side='right' gives, for each day, the index of the
        # last listed entry at or before it.
        pos = np.searchsorted(listed_days, target, side="right") - 1
        return values[pos]

    if rule is Fill.INTERPOLATE:
        # A single listed value is a constant, not an interpolation, so it is
        # held flat wherever it sits.  lumprem2.f has this branch too but
        # cannot reach it for a run longer than a day, because it separately
        # insists the last listed day reach the end of the simulation.
        if len(values) == 1:
            return np.broadcast_to(values[0], shape).copy()
        if listed_days[0] > 1.0 or listed_days[-1] < float(len(days)):
            raise ForcingError(
                f"{name}: fill='interpolate' needs listed values on or outside "
                f"both ends of the run ({days[0].date()} to {days[-1].date()}); "
                f"it has {listed[0].date()} to {listed[-1].date()}. "
                f"Extend the series, or choose a different fill rule."
            )
        # lumprem2.f computes the reciprocal of the gap once and multiplies by
        # it, rather than dividing per day.  Reproduced so that a run driven by
        # sparse forcing still matches the reference bit for bit.
        seg = np.clip(np.searchsorted(listed_days, target, side="right") - 1,
                      0, len(values) - 2)
        d1 = listed_days[seg]
        iddiff = 1.0 / (listed_days[seg + 1] - d1)
        dfac = ((target - d1) * iddiff).reshape((-1,) + (1,) * (values.ndim - 1))
        out = values[seg] + (values[seg + 1] - values[seg]) * dfac
        return out

    raise ForcingError(f"{name}: unknown fill rule {rule!r}")


def _resolve_fill(fill: Mapping[str, str | Fill] | None) -> dict[str, Fill]:
    out = dict(FORTRAN_FILL)
    for name, rule in (fill or {}).items():
        if name not in VARIABLES:
            raise ForcingError(f"fill names {name!r}, which is not a forcing variable; "
                               f"known names are {list(VARIABLES)}")
        try:
            out[name] = Fill(rule)
        except ValueError:
            raise ForcingError(
                f"fill[{name!r}] = {rule!r} is not a rule; "
                f"choose from {[f.value for f in Fill]}"
            ) from None
    return out


def _check_constants(constants: Mapping[str, object]) -> dict[str, float]:
    out = {}
    for name, value in constants.items():
        if name not in VARIABLES:
            raise ForcingError(f"unknown forcing variable {name!r}; "
                               f"known names are {list(VARIABLES)}")
        out[name] = float(value)  # type: ignore[arg-type]
    return out


def _first_bad(arr: np.ndarray, mask: np.ndarray) -> str:
    """Where the first offending value sits, and what it is."""
    if arr.ndim == 1:
        i = int(np.argmax(mask))
        return f"day {i + 1} is {arr[i]!r}"
    i, c = np.unravel_index(int(np.argmax(mask)), mask.shape)
    return f"day {i + 1}, cell {c} is {arr[i, c]!r}"


def _validate(daily: Mapping[str, np.ndarray]) -> None:
    """Reject forcing the physics cannot take, naming the day that broke it.

    Arrays are ``(ndays,)``, or ``(ndays, ncell)`` for per-cell forcing, in
    which case the cell is named too.
    """
    for name in ("rainfall", "pot_evap", "crop_factor", "pot_evap_lower"):
        bad = daily[name] < 0.0
        if bad.any():
            raise ForcingError(f"{name} must be >= 0; {_first_bad(daily[name], bad)}")
    bad = ~(daily["veg_gamma"] > 0.0)
    if bad.any():
        raise ForcingError(
            f"veg_gamma must be > 0 -- the evaporation curve is 0/0 at zero, "
            f"where lumprem2.f returns a NaN without saying so; "
            f"{_first_bad(daily['veg_gamma'], bad)}"
        )
    frac = daily["gw_irrig_frac"]
    bad = (frac < 0.0) | (frac > 1.0)
    if bad.any():
        raise ForcingError(f"gw_irrig_frac must be in [0, 1]; {_first_bad(frac, bad)}")
    code = daily["irrigate"]
    bad = (code != 0.0) & (code != 1.0)
    if bad.any():
        raise ForcingError(
            f"irrigate must be 0 or 1; {_first_bad(code, bad)}. A fractional "
            f"value usually means an interpolating fill rule was applied to it."
        )
    for name, arr in daily.items():
        bad = ~np.isfinite(arr)
        if bad.any():
            raise ForcingError(f"{name} is not finite: {_first_bad(arr, bad)}")
