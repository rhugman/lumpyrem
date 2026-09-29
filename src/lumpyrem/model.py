"""The model objects: parameters in, forcing in, labelled results out.

This is the layer the Fortran never had.  ``lumprem2.f`` takes a fixed-format
control file, integer simulation days and four day-indexed forcing files; this
takes typed parameter objects, real dates and any series pandas can hold, and
hands back a table you can index by time.

:class:`Model` is one site.  :class:`ModelGrid` is many independent ones --
cells of a groundwater model, say -- stepped together by the compiled kernel
in parallel.  Both marshal through the same functions below into
``lumpyrem.engine``, and a grid adds no arithmetic of its own, so a cell of a
grid reproduces the same :class:`Model` run alone exactly.

Time convention: a run is a whole number of days beginning at ``start``.  Day 1
covers ``[start, start + 1 day)``, so an output instant ``t`` marks the *end* of
day ``(t - start).days``.  Row 0 of the results sits at ``start`` itself and
carries the initial state.  Flux columns are therefore totals over the interval
ending at their own timestamp, which is the same convention a MODFLOW stress
period uses.
"""

from __future__ import annotations

import warnings
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field, replace

import numpy as np
import pandas as pd

from . import engine as _engine
from .forcing import Forcing, GridForcing
from .parameters import (
    InitialState,
    LowerStore,
    Solver,
    UpperStore,
    VolumeToElevation,
    _require,
)
from .results import COLUMNS, ConvergenceWarning, GridResults, Results

__all__ = ["Model", "ModelGrid", "ModelError"]

_DAY = pd.Timedelta(days=1)

# Used when a model carries no volume-to-elevation conversion: the kernel always
# computes the two elevation columns, so it is given a conversion that cannot
# fail and the columns are blanked afterwards.
_NO_ELEVATION = VolumeToElevation(offset=0.0, factor1=0.0, factor2=0.0,
                                  power=1.0, datum=0.0)
_ELEV = COLUMNS.index("elevation")
_DWT = COLUMNS.index("depth_to_water")


class ModelError(ValueError):
    """The model specification is inconsistent, or the run was ill-posed."""


@dataclass(frozen=True)
class Model:
    """A LUMPREM site: its stores, its solver settings and its initial state.

    Args:
        upper: the root-zone store.  Required.
        lower: the store below the root zone, or ``None`` for the one-store
            configuration (LUMPREM v1).
        elevation: volume-to-elevation conversion.  Optional; without one the
            ``elevation`` and ``depth_to_water`` columns come back as NaN.
        solver: sub-stepping and iteration settings.
        initial: starting volumes and delay-buffer contents.

    Example:
        >>> m = Model(upper=UpperStore(maxvol=0.3, ks=0.05, m=0.4, l=0.5))
        >>> results = m.run(forcing, times="MS")          # doctest: +SKIP
    """

    upper: UpperStore
    lower: LowerStore | None = None
    elevation: VolumeToElevation | None = None
    solver: Solver = field(default_factory=Solver)
    initial: InitialState = field(default_factory=InitialState)

    def __post_init__(self) -> None:
        for name, expect in (("upper", UpperStore), ("solver", Solver),
                             ("initial", InitialState)):
            value = getattr(self, name)
            if not isinstance(value, expect):
                raise ModelError(
                    f"Model.{name} must be a {expect.__name__}, "
                    f"got {type(value).__name__}"
                )
        if self.lower is not None and not isinstance(self.lower, LowerStore):
            raise ModelError("Model.lower must be a LowerStore or None, "
                             f"got {type(self.lower).__name__}")
        if self.elevation is not None and not isinstance(self.elevation, VolumeToElevation):
            raise ModelError("Model.elevation must be a VolumeToElevation or None, "
                             f"got {type(self.elevation).__name__}")

        init, upper = self.initial, self.upper
        _require(init.vol <= upper.maxvol, "InitialState", "vol", init.vol,
                 f"be <= UpperStore.maxvol ({upper.maxvol!r})")
        if self.lower is None:
            _require(init.vol_lower == 0.0, "InitialState", "vol_lower", init.vol_lower,
                     "be 0 when the model has no lower store")
        else:
            _require(init.vol_lower <= self.lower.capacity, "InitialState", "vol_lower",
                     init.vol_lower,
                     f"be <= LowerStore.maxvol + extravol ({self.lower.capacity!r})")
        if self.elevation is not None and self.elevation.bucket == "lower" \
                and self.lower is None:
            raise ModelError(
                "VolumeToElevation.bucket='lower' needs a lower store; the model "
                "has none. Pass lower=LowerStore(...) or bucket='upper'."
            )

    # ------------------------------------------------------------------
    @property
    def two_store(self) -> bool:
        return self.lower is not None

    def with_(self, **changes) -> Model:
        """A copy with some fields replaced -- handy inside a calibration loop."""
        return replace(self, **changes)

    # ------------------------------------------------------------------
    def run(self, forcing, *, times=None, start=None, end=None, engine: str = "auto",
            **forcing_kwargs) -> Results:
        """Run the model and return labelled results.

        Args:
            forcing: a :class:`~lumpyrem.forcing.Forcing`, a DataFrame, a
                mapping or a path to a CSV.
            times: when to report.  ``None`` reports every day; a pandas
                frequency string such as ``"MS"`` reports on that schedule; a
                sequence of dates reports on exactly those; a sequence of
                integers is taken as 1-based simulation days.
            start: first simulated day.  Defaults to the first forcing date.
            end: last simulated day.  Defaults to the last forcing date.
            engine: ``"auto"`` (the compiled kernel when Numba is installed),
                ``"compiled"`` or ``"python"``.  Results are identical; only
                the speed differs.
            **forcing_kwargs: passed to the Forcing constructor when ``forcing``
                is not already one.

        Returns:
            :class:`~lumpyrem.results.Results`.

        Raises:
            ModelError: the run is ill-posed -- no days, an output time off the
                day grid or outside the run, or non-convergence under
                ``Solver(on_nonconvergence="raise")``.
        """
        forcing = Forcing.coerce(forcing, **forcing_kwargs)
        start, ndays = _run_window(forcing, start, end)
        daily = forcing.to_daily(start, start + (ndays - 1) * _DAY)
        outdays = _output_days(times, start, ndays)
        run = _run_kernel((self,), _kernel_forcing(daily), outdays, self.solver, engine)

        results = Results(
            times=_stamps(start, run.days),
            days=run.days,
            values=run.values[0],
            nonconverged_upper=int(run.nonconverged[0, 0]),
            nonconverged_lower=int(run.nonconverged[0, 1]),
            has_elevation=self.elevation is not None,
        )
        if not results.converged:
            _report_convergence(self.solver, results.nonconverged_upper,
                                results.nonconverged_lower)
        return results


# ---------------------------------------------------------------------------
# Many cells
# ---------------------------------------------------------------------------

@dataclass(frozen=True, eq=False)
class ModelGrid:
    """Independent LUMPREM cells, run together.

    Each cell is a full :class:`Model` -- its own stores, conversion and
    initial state -- and nothing passes between cells.  What they share is
    what the kernel steps them with: one :class:`Solver`, one run window and
    one output schedule.  Forcing is shared too, or per cell through a
    :class:`~lumpyrem.forcing.GridForcing`.

    Build one from parameter arrays with :meth:`from_arrays`, or from a
    sequence of models directly.

    Args:
        cells: one :class:`Model` per cell, all with the same ``solver``.
        names: labels for the ``cell`` axis of the results; ``0 .. N-1`` by
            default.

    Example:
        >>> grid = ModelGrid.from_arrays(
        ...     upper=dict(maxvol=[0.3, 0.4, 0.5], ks=0.05, m=0.4, l=0.5))
        >>> grid.run(forcing, times="MS")["total_rech"]       # doctest: +SKIP
    """

    cells: tuple[Model, ...]
    names: tuple | None = None

    def __post_init__(self) -> None:
        cells = tuple(self.cells)
        if not cells:
            raise ModelError("a ModelGrid needs at least one cell")
        for i, cell in enumerate(cells):
            if not isinstance(cell, Model):
                raise ModelError(f"cell {i} must be a Model, got {type(cell).__name__}")
        solver = cells[0].solver
        for i, cell in enumerate(cells):
            if cell.solver != solver:
                raise ModelError(
                    f"cell {i} has {cell.solver} but cell 0 has {solver}; every cell "
                    f"of a grid is stepped with one solver"
                )
        names = tuple(range(len(cells))) if self.names is None else tuple(self.names)
        if len(names) != len(cells):
            raise ModelError(f"{len(names)} names for {len(cells)} cells")
        if len(set(names)) != len(names):
            raise ModelError("cell names must be unique")
        object.__setattr__(self, "cells", cells)
        object.__setattr__(self, "names", names)

    @classmethod
    def from_arrays(
        cls,
        *,
        upper: Mapping[str, object],
        lower: Mapping[str, object] | None = None,
        elevation: Mapping[str, object] | None = None,
        initial: Mapping[str, object] | None = None,
        solver: Solver | None = None,
        names: Sequence | None = None,
        ncell: int | None = None,
    ) -> ModelGrid:
        """Build from per-cell parameter arrays.

        Each mapping holds the fields of the matching parameter object --
        :class:`UpperStore`, :class:`LowerStore`, :class:`VolumeToElevation`,
        :class:`InitialState` -- as a number (or string) shared by every
        cell, or a 1-D array with one value per cell.  The two initial
        buffers are sequences already, so for them a 1-D sequence is shared
        and a sequence of sequences, one per cell, is not.

        Every cell is validated as a :class:`Model` would be, and an error
        names the cell.

        Args:
            upper, lower, elevation, initial: the parameter fields.
            solver: shared by every cell.
            names: labels for the cells.
            ncell: the number of cells, needed only when every value is
                shared.
        """
        groups = {"upper": upper, "lower": lower, "elevation": elevation,
                  "initial": initial}
        n = _count_cells(groups, names, ncell)
        labels = tuple(range(n)) if names is None else tuple(names)
        makers = {"upper": UpperStore, "lower": LowerStore,
                  "elevation": VolumeToElevation, "initial": InitialState}
        solver = Solver() if solver is None else solver
        cells = []
        for i in range(n):
            try:
                parts = {group: makers[group](**_cell_fields(fields, i))
                         for group, fields in groups.items() if fields is not None}
                cells.append(Model(solver=solver, **parts))
            except ValueError as exc:
                raise ModelError(f"cell {labels[i]!r}: {exc}") from exc
        return cls(cells=tuple(cells), names=labels)

    # ------------------------------------------------------------------
    @property
    def ncell(self) -> int:
        return len(self.cells)

    @property
    def solver(self) -> Solver:
        return self.cells[0].solver

    def __len__(self) -> int:
        return len(self.cells)

    def __getitem__(self, i: int) -> Model:
        """The ``i``-th cell's model, by position."""
        return self.cells[i]

    def __repr__(self) -> str:
        two = sum(c.two_store for c in self.cells)
        return f"<ModelGrid: {self.ncell} cells, {two} with a lower store>"

    # ------------------------------------------------------------------
    def run(self, forcing, *, times=None, start=None, end=None, engine: str = "auto",
            **forcing_kwargs) -> GridResults:
        """Run every cell and return results on ``(time, cell, variable)``.

        Args:
            forcing: shared by every cell -- anything :meth:`Model.run` takes
                -- or a :class:`~lumpyrem.forcing.GridForcing` with one column
                per cell.
            times, start, end, engine: as for :meth:`Model.run`.  The
                compiled engine runs cells on every core Numba is allowed;
                set ``NUMBA_NUM_THREADS`` to use fewer.
            **forcing_kwargs: passed to the Forcing constructor when
                ``forcing`` is not already one.

        Returns:
            :class:`~lumpyrem.results.GridResults`.
        """
        if isinstance(forcing, GridForcing):
            if forcing_kwargs:
                raise ModelError("options cannot be applied to an existing GridForcing")
            if forcing.ncell != self.ncell:
                raise ModelError(f"the forcing has {forcing.ncell} cells and the grid "
                                 f"has {self.ncell}")
            if forcing.cells is not None and forcing.cells != self.names:
                raise ModelError(
                    "the forcing's cell labels do not match the grid's names, in "
                    "order; pass names=forcing.cells when building the grid, or "
                    "reorder the forcing"
                )
        else:
            forcing = Forcing.coerce(forcing, **forcing_kwargs)
        start, ndays = _run_window(forcing, start, end)
        daily = forcing.to_daily(start, start + (ndays - 1) * _DAY)
        outdays = _output_days(times, start, ndays)
        run = _run_kernel(self.cells, _kernel_forcing(daily), outdays, self.solver, engine)

        results = GridResults(
            times=_stamps(start, run.days),
            days=run.days,
            values=run.values.transpose(1, 0, 2),
            cells=self.names,
            nonconverged=run.nonconverged,
            has_elevation=np.array([c.elevation is not None for c in self.cells]),
        )
        if not results.converged:
            failed = results.nonconverged_cells
            upper, lower = (int(v) for v in run.nonconverged.sum(axis=0))
            _report_convergence(self.solver, upper, lower,
                                where=f" across {len(failed)} of {self.ncell} cells "
                                      f"(first: {failed[0]!r})")
        return results


def _count_cells(groups: Mapping[str, Mapping | None], names, ncell) -> int:
    """How many cells the per-cell arrays describe, checked for agreement."""
    lengths = {}
    for group, fields in groups.items():
        for name, value in (fields or {}).items():
            if _per_cell(name, value):
                lengths[f"{group}.{name}"] = len(value)
    if names is not None:
        lengths["names"] = len(names)
    if ncell is not None:
        lengths["ncell"] = int(ncell)
    counts = set(lengths.values())
    if len(counts) > 1:
        raise ModelError("per-cell values disagree on the number of cells: " + ", ".join(
            f"{k}={v}" for k, v in lengths.items()))
    if not counts:
        raise ModelError("every value is shared, so the number of cells is unknown; "
                         "pass ncell= or names=")
    n = counts.pop()
    if n < 1:
        raise ModelError("a ModelGrid needs at least one cell")
    return n


def _per_cell(name: str, value) -> bool:
    """Is this field given per cell, rather than shared?"""
    if name in ("drain_buffer", "macro_buffer"):
        return len(value) > 0 and np.ndim(value[0]) == 1
    if value is None or isinstance(value, str):
        return False
    ndim = np.ndim(value)
    if ndim > 1:
        raise ModelError(f"{name}: expected a number or a 1-D array, got {ndim}-D")
    return ndim == 1


def _cell_fields(fields: Mapping[str, object], i: int) -> dict:
    """Cell ``i``'s values of a group of fields."""
    return {name: (value[i] if _per_cell(name, value) else value)
            for name, value in fields.items()}


# ---------------------------------------------------------------------------
# Marshalling into the kernel -- shared by Model and ModelGrid
# ---------------------------------------------------------------------------

#: Kernel argument name -> forcing variable.
_FORCING_ARGS = {
    "rain": "rainfall",
    "epot": "pot_evap",
    "cropfac": "crop_factor",
    "gamma": "veg_gamma",
    "irrigcode": "irrigate",
    "gwirrigfrac": "gw_irrig_frac",
    "epot_br": "pot_evap_lower",
}


def _run_window(forcing, start, end) -> tuple[pd.Timestamp, int]:
    """The run's first day and its length in days."""
    first, last = forcing.span
    start = pd.Timestamp(first if start is None else start)
    end = pd.Timestamp(last if end is None else end)
    if start != start.normalize() or end != end.normalize():
        raise ModelError("start and end must fall on midnight; the model's "
                         "step is a whole calendar day")
    ndays = int((end - start) // _DAY) + 1
    if ndays < 1:
        raise ModelError(f"the run ends ({end.date()}) before it starts "
                         f"({start.date()})")
    return start, ndays


def _kernel_params(model: Model) -> dict:
    """One model's parameters under the kernel's names."""
    elevation = model.elevation or _NO_ELEVATION
    upper, lower, init = model.upper, model.lower, model.initial
    return dict(
        maxvol=upper.maxvol, irrigvolfrac=upper.irrigvolfrac,
        rdelay=upper.rdelay, mdelay=upper.mdelay,
        ks=upper.ks, m=upper.m, l=upper.l, mflowmax=upper.mflowmax,
        maxvol_br=lower.maxvol if lower else 0.0,
        extravol_br=lower.extravol if lower else 0.0,
        gamma_br=lower.gamma if lower else 0.0,
        ks_br=lower.ks if lower else 0.0,
        m_br=lower.m if lower else 0.0,
        l_br=lower.l if lower else 0.0,
        offset=elevation.offset, factor1=elevation.factor1,
        factor2=elevation.factor2, power=elevation.power,
        datum=elevation.datum, bucket=elevation.bucket,
        elevmin=-1.0e20 if elevation.elevmin is None else elevation.elevmin,
        elevmax=1.0e20 if elevation.elevmax is None else elevation.elevmax,
        vol=init.vol, vol_br=init.vol_lower,
    )


def _kernel_forcing(daily: Mapping[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Daily forcing as the kernel takes it: ``(1, ndays)`` when every cell
    shares a variable, ``(ncell, ndays)`` when each has its own."""
    out = {}
    for arg, name in _FORCING_ARGS.items():
        a = daily[name]
        a = a.reshape(1, -1) if a.ndim == 1 else a.T
        dtype = np.int64 if arg == "irrigcode" else np.float64
        out[arg] = np.ascontiguousarray(a, dtype=dtype)
    return out


def _run_kernel(models: Sequence[Model], forcing, outdays, solver: Solver,
                engine: str) -> _engine.CellRun:
    """Pack the models into rows, run them, and blank what has no conversion."""
    if engine not in _engine.ENGINES:
        raise ModelError(f"engine must be one of {list(_engine.ENGINES)}, got {engine!r}")
    params = np.stack([_engine.param_row(**_kernel_params(m)) for m in models])
    rbuf = _engine.buffer_rows([m.initial.drain_buffer for m in models])
    mbuf = _engine.buffer_rows([m.initial.macro_buffer for m in models])
    run = _engine.run_cells(params, rbuf, mbuf, forcing, outdays, nstep=solver.nstep,
                            mxiter=solver.mxiter, tol=solver.tol, engine=engine)
    for i, model in enumerate(models):
        if model.elevation is None:
            run.values[i, :, _ELEV] = np.nan
            run.values[i, :, _DWT] = np.nan
    return run


def _stamps(start: pd.Timestamp, days: np.ndarray) -> pd.DatetimeIndex:
    return pd.DatetimeIndex([start + int(d) * _DAY for d in days], name="time")


def _report_convergence(solver: Solver, upper: int, lower: int, where: str = "") -> None:
    """Warn or raise, as the solver asks, about sub-steps that hit ``mxiter``."""
    detail = (f"{upper} upper-store and {lower} lower-store sub-step(s){where} hit "
              f"mxiter={solver.mxiter} without meeting tol={solver.tol:g}")
    advice = (
        "The drainage law has an unbounded slope as a store approaches full "
        "whenever m < 1, so raising nstep can make this worse rather than "
        "better. Try a larger m, a smaller ks, or a looser tol."
    )
    policy = solver.on_nonconvergence
    if policy == "raise":
        raise ModelError(f"the run did not converge: {detail}. {advice}")
    if policy == "warn":
        warnings.warn(f"{detail}. Results are returned but are not solutions. "
                      f"{advice}", ConvergenceWarning, stacklevel=3)


# ---------------------------------------------------------------------------
# Output times
# ---------------------------------------------------------------------------

def _output_days(times, start: pd.Timestamp, ndays: int) -> list[int]:
    """Resolve ``times`` to 1-based simulation day numbers."""
    if times is None:
        return list(range(1, ndays + 1))

    if isinstance(times, str):
        stamps = pd.date_range(start, start + ndays * _DAY, freq=times)
        days = sorted({int((t - start) // _DAY)
                       for t in stamps
                       if t == t.normalize() and start < t <= start + ndays * _DAY})
        if not days:
            raise ModelError(
                f"times={times!r} yields no output time inside the run "
                f"({start.date()} to {(start + ndays * _DAY).date()})"
            )
        return days

    if isinstance(times, Mapping) or not isinstance(times, Iterable):
        raise ModelError(
            f"times must be None, a frequency string, or a sequence of dates or "
            f"day numbers; got {type(times).__name__}"
        )

    items = list(times)
    if not items:
        raise ModelError("times is empty; pass None to report every day")

    if all(isinstance(t, (int, np.integer)) and not isinstance(t, bool) for t in items):
        days = [int(t) for t in items]
    else:
        days = []
        for t in pd.DatetimeIndex(pd.to_datetime(items)):
            if t != t.normalize():
                raise ModelError(
                    f"output time {t} is not midnight; the model reports at day "
                    f"boundaries only"
                )
            offset = (t - start) // _DAY
            if t != start + int(offset) * _DAY:
                raise ModelError(f"output time {t} is not a whole number of days "
                                 f"after the start of the run ({start})")
            days.append(int(offset))

    for day, original in zip(days, items, strict=True):
        if not 1 <= day <= ndays:
            where = original if isinstance(original, (int, np.integer)) else \
                f"{original} (day {day})"
            raise ModelError(
                f"output time {where} falls outside the run, which covers days "
                f"1 to {ndays} ({start.date()} to "
                f"{(start + ndays * _DAY).date()})"
            )
    if any(b <= a for a, b in zip(days, days[1:], strict=False)):
        raise ModelError("output times must be strictly increasing")
    return days
