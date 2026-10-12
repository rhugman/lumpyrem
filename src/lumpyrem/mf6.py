"""MODFLOW 6 through flopy: recharge cells from DISV, stress periods from TDIS.

A ``ModelGrid`` knows nothing about the groundwater model it feeds.  This
module is the bridge, in four steps:

1. :func:`top_active_cells` reads a flopy DISV package and returns the
   uppermost active cell of every column -- where recharge arrives.
2. :func:`stress_periods` reads a flopy TDIS package (and STO, for which
   periods are steady state) and returns each period's dates, plus the window
   its recharge rate is averaged over.
3. :func:`output_times` turns those windows into the ``times`` a grid run
   reports at, and :func:`period_rates` turns the run's interval totals into
   one rate per stress period per cell -- depth per day, as RCH takes it.
4. :func:`recharge_package` creates or updates the model's RCH package from
   those rates.

The time conventions line up without adjustment.  A stress period covers
``[start, end)``; a results row labelled ``t`` is the total over the days
ending at ``t``, so the row at a period's ``end`` is that period's total.

flopy is imported only when a function here is given flopy objects, so the
rest of the package does not depend on it.
"""

from __future__ import annotations

import os
import re
import shlex
import warnings
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from .results import GridResults

__all__ = [
    "output_times",
    "period_rates",
    "recharge_package",
    "stress_periods",
    "top_active_cells",
    "write_recharge_file",
]

_DAY = pd.Timedelta(days=1)

#: MODFLOW 6's LENBOUNDNAME.
_MAX_BOUNDNAME = 40

#: TDIS ``TIME_UNITS`` -> days.  LUMPREM steps whole days, so periods must
#: come out as whole days whatever the model's unit.
_DAYS_PER_UNIT = {
    "seconds": 1.0 / 86400.0, "minutes": 1.0 / 1440.0, "hours": 1.0 / 24.0,
    "days": 1.0, "years": 365.25,
}


class MF6Error(ValueError):
    """The MODFLOW model and the recharge run do not fit together."""


# ---------------------------------------------------------------------------
# Cells
# ---------------------------------------------------------------------------

def top_active_cells(disv, idomain: np.ndarray | None = None) -> pd.DataFrame:
    """The uppermost active cell of every column of a DISV grid.

    Args:
        disv: a flopy ``ModflowGwfdisv`` package.
        idomain: ``(nlay, ncpl)`` to use instead of the package's own;
            cells with ``idomain > 0`` are active.  Without either, every
            cell is.

    Returns:
        One row per column that has an active cell, in ``icpl`` order, with
        0-based ``layer`` and ``icpl``, the ``node`` number
        (``layer * ncpl + icpl``), the cell-centre ``x`` and ``y`` from
        CELL2D, and the cell ``top``.  Columns with no active cell are left
        out.
    """
    nlay = int(disv.nlay.get_data())
    ncpl = int(disv.ncpl.get_data())
    if idomain is None:
        idomain = disv.idomain.array
    idomain = np.ones((nlay, ncpl), dtype=int) if idomain is None \
        else np.asarray(idomain).reshape(nlay, ncpl)

    active = idomain > 0
    has_cell = active.any(axis=0)
    layer = np.argmax(active, axis=0)
    icpl = np.flatnonzero(has_cell)
    layer = layer[icpl]

    cell2d = disv.cell2d.array
    order = np.argsort(cell2d["icell2d"])
    xc = np.asarray(cell2d["xc"], dtype=float)[order]
    yc = np.asarray(cell2d["yc"], dtype=float)[order]

    top = np.asarray(disv.top.array, dtype=float).reshape(ncpl)
    botm = np.asarray(disv.botm.array, dtype=float).reshape(nlay, ncpl)
    tops = np.vstack([top[None, :], botm[:-1]])
    return pd.DataFrame({
        "layer": layer.astype(np.int64),
        "icpl": icpl.astype(np.int64),
        "node": (layer * ncpl + icpl).astype(np.int64),
        "x": xc[icpl],
        "y": yc[icpl],
        "top": tops[layer, icpl],
    })


# ---------------------------------------------------------------------------
# Stress periods
# ---------------------------------------------------------------------------

def _steady_flags(sto, nper: int) -> np.ndarray:
    """Which periods STO makes steady state.  MF6's rules: a period's setting
    persists until a later one changes it, and a model with no STO package,
    or none for period 1, starts steady state."""
    flags = np.ones(nper, dtype=bool)
    if sto is None:
        return flags
    changes = {}
    for name, value in (("steady_state", True), ("transient", False)):
        data = getattr(sto, name)
        for key, *_ in data.get_active_key_list():
            if data.get_data(key):
                changes[int(key)] = value
    state = True
    for k in range(nper):
        state = changes.get(k, state)
        flags[k] = state
    return flags


def stress_periods(
    tdis,
    sto=None,
    *,
    start=None,
    steady_state: Sequence[int] | None = None,
    steady_state_window=None,
) -> pd.DataFrame:
    """Stress-period dates from a flopy TDIS package.

    Args:
        tdis: a flopy ``ModflowTdis`` package (``sim.tdis``).
        sto: the model's flopy ``ModflowGwfsto``, to read which periods are
            steady state.  Without it, or ``steady_state``, every period is
            transient.
        start: the date the model starts.  Defaults to TDIS
            ``START_DATE_TIME``; one or the other is required.
        steady_state: 0-based periods to treat as steady state, instead of
            reading STO.
        steady_state_window: what a steady-state period's recharge is the
            average of.  A ``(start, end)`` pair applies to every
            steady-state period; a mapping ``{kper: (start, end)}`` sets them
            one by one.  A steady-state period without one is averaged over
            its own dates, like a transient one -- rarely what a pre-
            development steady state means.

    Returns:
        One row per period: ``kper`` (0-based), ``start`` and ``end``
        (``[start, end)``), ``perlen`` in days, ``nstp``, ``tsmult``,
        ``steady``, and ``rate_start``/``rate_end`` -- the window its rate
        is averaged over.

    Raises:
        MF6Error: no start date, or a period boundary off a whole day, since
            a LUMPREM step is one calendar day.
    """
    units = str(tdis.time_units.get_data() or "days").lower()
    if units not in _DAYS_PER_UNIT:
        raise MF6Error(f"TDIS TIME_UNITS {units!r} is not one of {list(_DAYS_PER_UNIT)}")
    if start is None:
        start = tdis.start_date_time.get_data()
        if start is None:
            raise MF6Error("TDIS has no START_DATE_TIME; pass start=")
    start = pd.Timestamp(start)
    if start != start.normalize():
        raise MF6Error(f"the model starts at {start}, not on midnight; LUMPREM steps "
                       f"whole calendar days")

    perioddata = tdis.perioddata.array
    perlen = np.asarray(perioddata["perlen"], dtype=float) * _DAYS_PER_UNIT[units]
    days = np.round(perlen)
    if np.any(np.abs(perlen - days) > 1e-6) or np.any(days < 1):
        k = int(np.flatnonzero((np.abs(perlen - days) > 1e-6) | (days < 1))[0])
        raise MF6Error(f"stress period {k + 1} is {perlen[k]!r} days long; LUMPREM "
                       f"steps whole days, so every period must be a whole number of "
                       f"days")
    nper = len(days)
    edges = start + pd.to_timedelta(np.concatenate([[0.0], np.cumsum(days)]), unit="D")

    if steady_state is not None:
        steady = np.zeros(nper, dtype=bool)
        steady[list(steady_state)] = True
    elif sto is not None:
        steady = _steady_flags(sto, nper)
    else:
        steady = np.zeros(nper, dtype=bool)

    frame = pd.DataFrame({
        "kper": np.arange(nper),
        "start": edges[:-1],
        "end": edges[1:],
        "perlen": days.astype(np.int64),
        "nstp": np.asarray(perioddata["nstp"], dtype=np.int64),
        "tsmult": np.asarray(perioddata["tsmult"], dtype=float),
        "steady": steady,
    })
    frame["rate_start"] = frame["start"]
    frame["rate_end"] = frame["end"]

    if steady_state_window is not None:
        if isinstance(steady_state_window, Mapping):
            windows = dict(steady_state_window)
        else:
            windows = {int(k): steady_state_window for k in frame.index[steady]}
        for k, (a, b) in windows.items():
            a, b = pd.Timestamp(a), pd.Timestamp(b)
            if not (a < b and a == a.normalize() and b == b.normalize()):
                raise MF6Error(f"steady-state window for period {k + 1} must be two "
                               f"midnights in order; got {a} to {b}")
            frame.loc[k, ["rate_start", "rate_end"]] = [a, b]
    return frame


def output_times(periods: pd.DataFrame, run_start=None) -> pd.DatetimeIndex:
    """The ``times`` a grid run must report at for :func:`period_rates`.

    Every rate window's two ends, less the run's own start (row 0 of any
    results).  Pass it straight to ``ModelGrid.run(times=...)``.
    """
    stamps = pd.DatetimeIndex(np.concatenate([periods["rate_start"].to_numpy(),
                                              periods["rate_end"].to_numpy()]))
    stamps = stamps.unique().sort_values()
    if run_start is not None:
        run_start = pd.Timestamp(run_start)
        if stamps[0] < run_start:
            raise MF6Error(f"a rate window starts on {stamps[0].date()}, before the "
                           f"recharge run does ({run_start.date()})")
        stamps = stamps[stamps > run_start]
    return stamps


def period_rates(results: GridResults, periods: pd.DataFrame, column: str = "net_recharge",
                 *, factor: float = 1.0) -> np.ndarray:
    """One recharge rate per stress period per cell.

    The total of ``column`` over each period's ``[rate_start, rate_end)``,
    divided by the window's length in days and multiplied by ``factor`` --
    ``0.001`` turns millimetres into metres, say.

    Args:
        results: a grid run that reported at :func:`output_times`.
        periods: from :func:`stress_periods`.
        column: a flux column of the results.
        factor: unit conversion applied to the rate.

    Returns:
        ``(nper, ncell)`` float64, depth per day.
    """
    values = results.column(column)
    times = results.times
    i0 = times.get_indexer(pd.DatetimeIndex(periods["rate_start"]))
    i1 = times.get_indexer(pd.DatetimeIndex(periods["rate_end"]))
    if (i0 < 0).any() or (i1 < 0).any():
        bad = periods.loc[(i0 < 0) | (i1 < 0)].iloc[0]
        raise MF6Error(
            f"the results do not report at {bad['rate_start'].date()} and "
            f"{bad['rate_end'].date()} (period {int(bad['kper']) + 1}); run the grid "
            f"with times=output_times(periods)"
        )
    length = ((periods["rate_end"] - periods["rate_start"]) / _DAY).to_numpy()
    out = np.empty((len(periods), values.shape[1]))
    for k, (a, b) in enumerate(zip(i0, i1, strict=True)):
        # Row b is the total over the interval ending at times[b]; a single
        # row is taken as it is, so a transient period's rate is exact.
        total = values[b] if b == a + 1 else values[a + 1:b + 1].sum(axis=0)
        out[k] = total / length[k] * factor
    return out


# ---------------------------------------------------------------------------
# The RCH package
# ---------------------------------------------------------------------------

def _cellids(cells, nlay: int, ncpl: int) -> tuple[np.ndarray, np.ndarray]:
    if isinstance(cells, pd.DataFrame):
        layer = cells["layer"].to_numpy(dtype=np.int64)
        icpl = cells["icpl"].to_numpy(dtype=np.int64)
    else:
        arr = np.asarray(cells, dtype=np.int64)
        if arr.ndim != 2 or arr.shape[1] != 2:
            raise MF6Error(f"cells must be a DataFrame with layer and icpl columns or an "
                           f"(ncell, 2) array of 0-based (layer, icpl); got {arr.shape}")
        layer, icpl = arr[:, 0], arr[:, 1]
    if layer.min() < 0 or layer.max() >= nlay or icpl.min() < 0 or icpl.max() >= ncpl:
        raise MF6Error(f"a cell lies outside the grid's {nlay} layers x {ncpl} cells; "
                       f"cell ids are 0-based")
    return layer, icpl


def _rch_filename(gwf, pname: str) -> str:
    """The file the model's name file lists for its RCH package, if any."""
    entries = gwf.name_file.packages.get_data()
    rch = [] if entries is None else [
        (str(fname), str(name)) for ftype, fname, name in entries
        if str(ftype).lower().startswith("rch")]
    for fname, name in rch:
        if name == pname:
            return fname
    if len(rch) == 1:
        return rch[0][0]
    return f"{gwf.name}.rch"


def _quote(token: str) -> str:
    """``token`` as MF6 reads it whole: single-quoted if it holds whitespace."""
    return f"'{token}'" if re.search(r"\s", token) else token


def _boundname_tokens(boundnames: Sequence[str]) -> list[str]:
    """Each boundname as it is written to the file, after checking MF6 can read it."""
    out = []
    for b in boundnames:
        b = str(b)
        if not b.strip():
            raise MF6Error("a boundname is empty")
        if len(b) > _MAX_BOUNDNAME:
            raise MF6Error(f"boundname {b!r} is {len(b)} characters; MODFLOW 6 reads at "
                           f"most {_MAX_BOUNDNAME}")
        if "'" in b:
            raise MF6Error(f"boundname {b!r} contains a single quote, which MODFLOW 6 "
                           f"cannot read inside a name")
        out.append(_quote(b))
    return out


def _register_package(gwf, ftype: str, filename: str, pname: str) -> None:
    """List a package file in the model's name file on disk, or MF6 never reads it.

    An entry for ``filename`` already there is left alone, and one for another
    file under the same ``pname`` is pointed at ``filename``.  Otherwise the
    entry is added at the end of the PACKAGES block.
    """
    sim_path = gwf.simulation.simulation_data.mfpath.get_sim_path()
    path = os.path.join(sim_path, gwf.model_nam_file)
    if not os.path.exists(path):
        raise MF6Error(f"the model's name file {path} does not exist; write the "
                       f"simulation before the RCH file")
    with open(path) as fid:
        lines = fid.read().splitlines()

    upper = [line.strip().upper() for line in lines]
    try:
        begin = next(i for i, u in enumerate(upper) if re.match(r"BEGIN\s+PACKAGES\b", u))
        end = next(i for i in range(begin + 1, len(lines))
                   if re.match(r"END\s+PACKAGES\b", upper[i]))
    except StopIteration:
        raise MF6Error(f"{path} has no BEGIN/END PACKAGES block") from None

    def norm(f):
        return os.path.normcase(os.path.normpath(f))

    replace, others = None, []
    for i in range(begin + 1, end):
        if upper[i][:1] in ("#", "!", "") or upper[i].startswith("//"):
            continue
        tokens = shlex.split(lines[i], comments=True)
        if len(tokens) < 2:
            continue
        kind, fname = tokens[0].upper(), tokens[1]
        name = tokens[2] if len(tokens) > 2 else None
        if kind == ftype.upper() and norm(fname) == norm(filename):
            return
        if name is not None and name.lower() == pname.lower():
            if kind != ftype.upper():
                raise MF6Error(f"the model already has a {kind} package named {pname!r} "
                               f"in {path}; pass another pname")
            replace = i
        elif kind == ftype.upper():
            others.append(fname)

    entry = f"  {ftype}  {_quote(filename)}  {pname}"
    if replace is None:
        lines.insert(end, entry)
    else:
        lines[replace] = entry
    if others:
        warnings.warn(f"{path} also lists {ftype} file(s) {others}; MODFLOW 6 applies "
                      f"those as well as {filename!r}", stacklevel=3)
    with open(path, "w", newline="\n") as fid:
        fid.write("\n".join(lines) + "\n")


def write_recharge_file(
    path,
    layer: np.ndarray,
    icpl: np.ndarray,
    rates: np.ndarray,
    *,
    boundnames: Sequence[str] | None = None,
    save_flows: bool = True,
    print_input: bool = False,
    obs_file: str | None = None,
    fmt: str = "%16.9e",
    header: str = "File generated by lumpyrem.mf6",
) -> str:
    """Write a list-based MF6 RCH file directly, one PERIOD block per row of ``rates``.

    The text flopy would write, without flopy's per-entry formatting: each
    period block is one ``%`` substitution into a template built once, which
    writes fifty thousand cells for five hundred periods in seconds rather
    than minutes.  ``layer`` and ``icpl`` are 0-based; the file is 1-based.
    Boundnames holding spaces are quoted, so MODFLOW reads each as one name.

    This writes the package file only; MODFLOW reads it once the model's name
    file lists it, which :func:`recharge_package` sees to.
    """
    rates = np.asarray(rates, dtype=np.float64)
    ncell = layer.size
    if rates.ndim != 2 or rates.shape[1] != ncell:
        raise MF6Error(f"rates must be (nper, {ncell}); got {rates.shape}")
    if not np.isfinite(rates).all():
        raise MF6Error("rates contain NaN or infinity")
    if boundnames is not None and len(boundnames) != ncell:
        raise MF6Error(f"{len(boundnames)} boundnames for {ncell} cells")

    # The template is %-formatted once per period, so a literal % in a
    # boundname has to be doubled.
    names = [""] * ncell if boundnames is None else \
        [f"  {b}".replace("%", "%%") for b in _boundname_tokens(boundnames)]
    template = "".join(f"  {lay + 1:<3d} {cell + 1:<8d} {fmt}{name}\n"
                       for lay, cell, name in zip(layer.tolist(), icpl.tolist(), names,
                                                  strict=True))
    with open(path, "w", newline="\n") as fid:
        fid.write(f"# {header}\n")
        fid.write("BEGIN options\n")
        if save_flows:
            fid.write("  SAVE_FLOWS\n")
        if print_input:
            fid.write("  PRINT_INPUT\n")
        if boundnames is not None:
            fid.write("  BOUNDNAMES\n")
        if obs_file:
            fid.write(f"  OBS6  FILEIN  {_quote(obs_file)}\n")
        fid.write("END options\n\n")
        fid.write(f"BEGIN dimensions\n  MAXBOUND  {ncell}\nEND dimensions\n\n")
        for k in range(rates.shape[0]):
            fid.write(f"BEGIN period  {k + 1}\n")
            fid.write(template % tuple(rates[k].tolist()))
            fid.write(f"END period  {k + 1}\n\n")
    return str(path)


def recharge_package(
    gwf,
    cells,
    rates: np.ndarray,
    *,
    boundnames: Sequence[str] | None = None,
    pname: str = "rch",
    filename: str | None = None,
    obs_file: str | None = None,
    save_flows: bool = True,
    print_input: bool = False,
    method: str = "flopy",
    write: bool = True,
    **flopy_kwargs,
):
    """Create or update a model's list-based RCH package from per-period rates.

    Everything the package needs that the model already knows is read from
    it: the grid's shape from DISV, the number of stress periods from the
    simulation's TDIS, and the RCH file name from the model's name file.

    Args:
        gwf: a flopy ``ModflowGwf`` with its DISV loaded.
        cells: where recharge is applied -- a DataFrame with 0-based
            ``layer`` and ``icpl`` columns, such as :func:`top_active_cells`
            returns, or an ``(ncell, 2)`` array of the same.
        rates: ``(nper, ncell)`` recharge rates, from :func:`period_rates`.
        boundnames: one name per cell, for RCH observations by name.
        pname: the package name to create or replace.
        filename: the RCH file.  Defaults to the one the name file lists, or
            ``<model>.rch``.
        obs_file: an existing RCH OBS6 file to reference from the options
            block, with either method.  With ``method="flopy"`` the file is
            referenced as it is; pass ``observations=`` instead to have flopy
            write one.
        save_flows, print_input: RCH options.
        method: ``"flopy"`` (the default) replaces the model's RCH package
            object with a new ``ModflowGwfrch`` carrying the rates, so the
            flopy model stays the one source of truth and a later
            ``sim.write_simulation()`` writes it with everything else.
            ``"fast"`` writes the file directly with
            :func:`write_recharge_file` and lists it in the model's name file
            on disk -- the same file, many times faster for tens of thousands
            of cells over hundreds of periods (``scripts/benchmark_rch.py``).
            It needs the simulation already written, and leaves the flopy
            objects alone: a loaded RCH package, if any, no longer matches its
            file, and a later ``sim.write_simulation()`` rewrites the name
            file without the entry -- so write the RCH file last.
        write: with ``method="flopy"``, write the package file and the
            model's name file now.  Pass ``False`` when the simulation is
            written afterwards, which would only write them again.
        **flopy_kwargs: passed to ``ModflowGwfrch`` with ``method="flopy"``,
            ``observations=`` among them.

    Returns:
        The file path (``"fast"``) or the new flopy package (``"flopy"``).
    """
    disv = gwf.get_package("disv")
    if disv is None:
        raise MF6Error("the model has no DISV package loaded")
    nlay, ncpl = int(disv.nlay.get_data()), int(disv.ncpl.get_data())
    layer, icpl = _cellids(cells, nlay, ncpl)
    rates = np.asarray(rates, dtype=np.float64)
    nper = int(gwf.simulation.tdis.nper.get_data())
    if rates.shape != (nper, layer.size):
        raise MF6Error(f"rates must be (nper, ncell) = ({nper}, {layer.size}) for this "
                       f"model; got {rates.shape}")
    filename = filename or _rch_filename(gwf, pname)

    loaded = gwf.get_package(pname)
    if loaded is not None and str(loaded.package_type).lower() != "rch":
        raise MF6Error(f"the model already has a {str(loaded.package_type).upper()}6 "
                       f"package named {pname!r}; pass another pname")

    if method == "fast":
        if loaded is not None:
            warnings.warn(f"the model's loaded {pname!r} package object is not updated by "
                          f"method='fast'; its file is", stacklevel=2)
        path = os.path.join(gwf.simulation.simulation_data.mfpath.get_sim_path(), filename)
        path = write_recharge_file(path, layer, icpl, rates, boundnames=boundnames,
                                   save_flows=save_flows, print_input=print_input,
                                   obs_file=obs_file)
        _register_package(gwf, "RCH6", filename, pname)
        return path
    if method != "flopy":
        raise MF6Error(f"method must be 'fast' or 'flopy', got {method!r}")
    if boundnames is not None:
        _boundname_tokens(boundnames)

    if obs_file and flopy_kwargs.get("observations") is not None:
        raise MF6Error("pass obs_file (an existing OBS6 file) or observations= (one "
                       "flopy writes), not both")
    import flopy

    if loaded is not None:
        gwf.remove_package(loaded)
    cellid = [(int(a), int(b)) for a, b in zip(layer.tolist(), icpl.tolist(), strict=True)]
    names = None if boundnames is None else [str(b) for b in boundnames]
    spd = {}
    for k in range(nper):
        if names is None:
            spd[k] = list(zip(cellid, rates[k].tolist(), strict=True))
        else:
            spd[k] = list(zip(cellid, rates[k].tolist(), names, strict=True))
    package = flopy.mf6.ModflowGwfrch(
        gwf, pname=pname, filename=filename, boundnames=boundnames is not None,
        save_flows=save_flows, print_input=print_input, maxbound=layer.size,
        stress_period_data=spd, **flopy_kwargs,
    )
    if obs_file:
        # ModflowGwfrch takes observations to write, not a file to point at;
        # its OBS6 FILEIN record is set directly.
        record = getattr(package, "_obs_filerecord", None)
        if record is None:
            raise MF6Error("this flopy version cannot reference an existing obs_file "
                           "with method='flopy'; use method='fast' or observations=")
        record.set_data([(obs_file,)])
    if write:
        # The package, the packages it owns (an OBS6 file from observations=),
        # and the name file that lists it: what MODFLOW needs to read it.
        package.write()
        for child in gwf.packagelist:
            if getattr(child, "parent_file", None) is package:
                child.write()
        gwf.name_file.write()
    return package
