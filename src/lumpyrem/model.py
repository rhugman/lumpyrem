"""The model object: parameters in, forcing in, labelled results out.

This is the layer the Fortran never had.  ``lumprem2.f`` takes a fixed-format
control file, integer simulation days and four day-indexed forcing files; this
takes typed parameter objects, real dates and any series pandas can hold, and
hands back a table you can index by time.

Time convention: a run is a whole number of days beginning at ``start``.  Day 1
covers ``[start, start + 1 day)``, so an output instant ``t`` marks the *end* of
day ``(t - start).days``.  Row 0 of the results sits at ``start`` itself and
carries the initial state.  Flux columns are therefore totals over the interval
ending at their own timestamp, which is the same convention a MODFLOW stress
period uses.
"""

from __future__ import annotations

import warnings
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field, replace

import numpy as np
import pandas as pd

from . import core
from .forcing import Forcing
from .parameters import (
    InitialState,
    LowerStore,
    Solver,
    UpperStore,
    VolumeToElevation,
    _require,
)
from .results import COLUMNS, ConvergenceWarning, Results

__all__ = ["Model", "ModelError"]

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
    def run(self, forcing, *, times=None, start=None, end=None, **forcing_kwargs) -> Results:
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

        daily = forcing.to_daily(start, end)
        outdays = _output_days(times, start, ndays)

        elevation = self.elevation or _NO_ELEVATION
        lower = self.lower
        sim = core.simulate(
            maxvol=self.upper.maxvol,
            irrigvolfrac=self.upper.irrigvolfrac,
            rdelay=self.upper.rdelay,
            mdelay=self.upper.mdelay,
            ks=self.upper.ks, m=self.upper.m, l=self.upper.l,
            mflowmax=self.upper.mflowmax,
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
            vol=self.initial.vol,
            vol_br=self.initial.vol_lower,
            rbuf=self.initial.drain_buffer,
            mbuf=self.initial.macro_buffer,
            nstep=self.solver.nstep, mxiter=self.solver.mxiter, tol=self.solver.tol,
            # Python floats, not NumPy scalars.  Measured across all 209 golden
            # cases the two agree exactly on macOS/arm64, so this is insurance
            # rather than a known trap: `**` and `exp` reach libm by different
            # routes for the two types, and nothing guarantees every platform
            # rounds the last bit the same way on both. One list conversion per
            # run is a cheap price for taking that off the table.
            rain=daily["rainfall"].tolist(),
            epot=daily["pot_evap"].tolist(),
            cropfac=daily["crop_factor"].tolist(),
            gamma=daily["veg_gamma"].tolist(),
            irrigcode=[int(v) for v in daily["irrigate"]],
            gwirrigfrac=daily["gw_irrig_frac"].tolist(),
            epot_br=daily["pot_evap_lower"].tolist() if lower else None,
            outdays=outdays,
        )

        values = sim.values
        if self.elevation is None:
            values[:, _ELEV] = np.nan
            values[:, _DWT] = np.nan

        results = Results(
            times=pd.DatetimeIndex([start + int(d) * _DAY for d in sim.days], name="time"),
            days=sim.days,
            values=values,
            nonconverged_upper=sim.nonconverged_upper,
            nonconverged_lower=sim.nonconverged_lower,
            has_elevation=self.elevation is not None,
        )
        self._report_convergence(results)
        return results

    # ------------------------------------------------------------------
    def _report_convergence(self, results: Results) -> None:
        if results.converged:
            return
        detail = (f"{results.nonconverged_upper} upper-store and "
                  f"{results.nonconverged_lower} lower-store sub-step(s) hit "
                  f"mxiter={self.solver.mxiter} without meeting tol="
                  f"{self.solver.tol:g}")
        advice = (
            "The drainage law has an unbounded slope as a store approaches full "
            "whenever m < 1, so raising nstep can make this worse rather than "
            "better. Try a larger m, a smaller ks, or a looser tol."
        )
        policy = self.solver.on_nonconvergence
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
