"""Running cells: the parameter layout, and which kernel runs it.

Two kernels compute the same numbers.  ``lumpyrem.compiled`` is the Numba
build and steps cells in parallel; ``lumpyrem.core`` is the readable
translation it was written from, run here one cell after another.  They are
held to exact equality (``tests/test_compiled.py``), so choosing between them
is a choice of speed, never of results.  ``engine="auto"`` takes the compiled
kernel whenever Numba imports, and ``core`` otherwise, so the package runs
without Numba installed -- just slower: 6x the Fortran's time per cell on
the benchmark workload (``scripts/benchmark_numba.py``) against the compiled
kernel's 0.9x, and one core instead of all of them.

A NumPy kernel vectorised over cells is deliberately not the fallback.  The
plan measured it at 2-17x slower than one core of Fortran, and it could not
hold the bit-identity gate the other two share: on x86 with AVX-512, NumPy
evaluates ``exp`` with its own SIMD routines rather than the maths library's.

Everything a cell can differ in is packed into one float64 row, so an
``(ncell, NPARAM)`` array carries a whole grid.  Forcing is one array per
variable of shape ``(1, ntime)`` when every cell shares it or
``(ncell, ntime)`` when each has its own, chosen per variable.
"""

from __future__ import annotations

import os
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import cache

import numpy as np

from . import core

__all__ = ["ENGINES", "THREADS_ENV", "CellRun", "compiled_available", "resolve",
           "resolve_threads", "run_cells"]

#: Environment variable giving the default thread count for grid runs, so a
#: PEST worker can be held to one core without changing the script it runs.
THREADS_ENV = "LUMPYREM_THREADS"

# -- The parameter row.
MAXVOL, IRRIGVOLFRAC, RDELAY, MDELAY, KS, M, L, MFLOWMAX = range(8)
MAXVOL_BR, EXTRAVOL_BR, GAMMA_BR, KS_BR, M_BR, L_BR = range(8, 14)
OFFSET, FACTOR1, FACTOR2, POWER, DATUM, NBUCKET, ELEVMIN, ELEVMAX = range(14, 22)
VOL, VOL_BR = range(22, 24)
NPARAM = 24

#: Forcing variables in the order the kernels take them.
FORCING: tuple[str, ...] = (
    "rain", "epot", "cropfac", "gamma", "irrigcode", "gwirrigfrac", "epot_br",
)

ENGINES: tuple[str, ...] = ("auto", "compiled", "python")


def param_row(
    *, maxvol, irrigvolfrac, rdelay, mdelay, ks, m, l, mflowmax,
    maxvol_br=0.0, extravol_br=0.0, gamma_br=0.0, ks_br=0.0, m_br=0.0, l_br=0.0,
    offset, factor1, factor2, power, datum, bucket="upper",
    elevmin=-1.0e20, elevmax=1.0e20, vol, vol_br=0.0,
) -> np.ndarray:
    """Pack one cell's parameters -- ``core.simulate``'s names -- into a row."""
    if bucket not in ("upper", "lower"):
        raise ValueError(f"bucket must be 'upper' or 'lower', got {bucket!r}")
    if bucket == "lower" and maxvol_br <= 0.0:
        raise ValueError("bucket='lower' requires maxvol_br > 0")
    p = np.empty(NPARAM)
    p[MAXVOL], p[IRRIGVOLFRAC], p[RDELAY], p[MDELAY] = maxvol, irrigvolfrac, rdelay, mdelay
    p[KS], p[M], p[L], p[MFLOWMAX] = ks, m, l, mflowmax
    p[MAXVOL_BR], p[EXTRAVOL_BR], p[GAMMA_BR] = maxvol_br, extravol_br, gamma_br
    p[KS_BR], p[M_BR], p[L_BR] = ks_br, m_br, l_br
    p[OFFSET], p[FACTOR1], p[FACTOR2], p[POWER] = offset, factor1, factor2, power
    p[DATUM], p[ELEVMIN], p[ELEVMAX] = datum, elevmin, elevmax
    p[NBUCKET] = 1.0 if bucket == "upper" else 2.0
    p[VOL], p[VOL_BR] = vol, vol_br
    return p


def param_kwargs(p: np.ndarray) -> dict:
    """The inverse of :func:`param_row`, as Python floats."""
    v = p.tolist()
    return dict(
        maxvol=v[MAXVOL], irrigvolfrac=v[IRRIGVOLFRAC], rdelay=v[RDELAY],
        mdelay=v[MDELAY], ks=v[KS], m=v[M], l=v[L], mflowmax=v[MFLOWMAX],
        maxvol_br=v[MAXVOL_BR], extravol_br=v[EXTRAVOL_BR], gamma_br=v[GAMMA_BR],
        ks_br=v[KS_BR], m_br=v[M_BR], l_br=v[L_BR],
        offset=v[OFFSET], factor1=v[FACTOR1], factor2=v[FACTOR2], power=v[POWER],
        datum=v[DATUM], bucket="upper" if v[NBUCKET] == 1.0 else "lower",
        elevmin=v[ELEVMIN], elevmax=v[ELEVMAX], vol=v[VOL], vol_br=v[VOL_BR],
    )


def buffer_rows(buffers: Sequence[Sequence[float]]) -> np.ndarray:
    """Initial delay buffers as 1-based rows (index 0 unused), zero-padded.

    Padding is exact: every kernel sums a buffer from its newest element, and
    trailing zeros added to such a sum leave it unchanged.
    """
    width = max((len(b) for b in buffers), default=0)
    rows = np.zeros((len(buffers), width + 1))
    for i, b in enumerate(buffers):
        rows[i, 1:1 + len(b)] = b
    return rows


def nrows(outdays: np.ndarray, numdays: int) -> int:
    """Day 0 plus every output time up to the first that reaches the end."""
    reached = np.flatnonzero(outdays >= numdays)
    return 1 + (int(reached[0]) + 1 if reached.size else outdays.size)


# ---------------------------------------------------------------------------
# Choosing and running a kernel
# ---------------------------------------------------------------------------

@cache
def compiled_available() -> bool:
    """True when Numba imports and the compiled kernel can run.

    Cached: a failed import is not, and would otherwise be retried on every run.
    """
    try:
        from . import compiled  # noqa: F401
    except ImportError:
        return False
    return True


def resolve(engine: str) -> str:
    """``"compiled"`` or ``"python"``: which kernel ``engine`` asks for here."""
    if engine not in ENGINES:
        raise ValueError(f"engine must be one of {list(ENGINES)}, got {engine!r}")
    if engine == "python":
        return "python"
    if compiled_available():
        return "compiled"
    if engine == "compiled":
        raise ImportError("engine='compiled' needs Numba; install it with "
                          "`pip install lumpyrem[fast]`, or pass engine='auto'")
    return "python"


def resolve_threads(threads: int | None = None,
                    single_threaded: bool | None = None) -> int | None:
    """The thread count a grid run should use; ``None`` means all of them.

    ``single_threaded=True`` is ``threads=1``.  With neither given, the
    ``LUMPYREM_THREADS`` environment variable decides, and unset it is every
    thread Numba has (``NUMBA_NUM_THREADS``, one per logical core by default).

    One thread is the setting for parallel calibration -- PEST++ running as
    many model runs at once as there are cores -- where a grid run that
    spreads over every core competes with all the others for each of them.
    """
    if single_threaded:
        if threads not in (None, 1):
            raise ValueError(f"single_threaded=True contradicts threads={threads}")
        return 1
    if threads is None:
        env = os.environ.get(THREADS_ENV, "").strip()
        if not env:
            return None
        try:
            threads = int(env)
        except ValueError:
            raise ValueError(f"{THREADS_ENV}={env!r} is not an integer") from None
    if isinstance(threads, bool) or int(threads) != threads or threads < 1:
        raise ValueError(f"threads must be an integer >= 1 or None, got {threads!r}")
    return int(threads)


@dataclass(frozen=True)
class CellRun:
    """What the kernel returns for a set of cells sharing one output schedule."""

    days: np.ndarray            # int64, (nrows,), 0 first
    values: np.ndarray          # float64, (ncell, nrows, ncolumn): the columns kept
    nonconverged: np.ndarray    # int64, (ncell, 2): upper, lower


def run_cells(
    params: np.ndarray,
    rbuf: np.ndarray,
    mbuf: np.ndarray,
    forcing: Mapping[str, np.ndarray],
    outdays,
    *,
    nstep: int,
    mxiter: int,
    tol: float,
    zone: np.ndarray | None = None,
    vegclass: np.ndarray | None = None,
    cropfac_table: np.ndarray | None = None,
    gamma_table: np.ndarray | None = None,
    columns: Sequence[int] | None = None,
    engine: str = "auto",
    threads: int | None = None,
) -> CellRun:
    """Run ``params.shape[0]`` independent cells.

    Every cell reads its forcing through a zone: each variable is
    ``(1, ntime)``, shared by every cell, or ``(nzone, ntime)``, and cell
    ``c`` reads row ``zone[c]``.  Without ``zone`` each cell is its own zone,
    so ``(ncell, ntime)`` forcing is one row per cell.  Holding forcing per
    zone sizes it by the number of zones rather than the number of cells,
    which is what makes a grid of tens of thousands of cells over decades of
    days fit in memory.

    Args:
        params: ``(ncell, NPARAM)`` rows from :func:`param_row`.
        rbuf, mbuf: ``(ncell, width + 1)`` initial buffers from
            :func:`buffer_rows`.
        forcing: every name in :data:`FORCING`, each ``(1, ntime)`` or
            ``(nzone, ntime)``.
        outdays: 1-based output days, shared by every cell.
        nstep, mxiter, tol: solver settings, shared by every cell.
        zone: ``(ncell,)`` 0-based row index into the forcing; the identity
            when omitted.
        vegclass: optional ``(1 | nzone, ntime)`` integer vegetation-class
            codes.  With it, cell ``c`` takes its crop factor and gamma on day
            ``t`` from ``cropfac_table[c, k]`` and ``gamma_table[c, k]``,
            ``k = vegclass[zone[c], t]``, and ``forcing["cropfac"]`` and
            ``forcing["gamma"]`` are not read.
        cropfac_table, gamma_table: ``(ncell, nclass)``, required with
            ``vegclass``.
        columns: kernel column indices to keep; all of them by default.
            ``values`` is then ``(ncell, nrows, len(columns))``, and the rest
            are never stored.
        engine: ``"auto"``, ``"compiled"`` or ``"python"``.
        threads: threads for the compiled engine, from
            :func:`resolve_threads`; ``1`` never starts the Numba thread
            pool.  The Python engine always runs on one.
    """
    ncell = params.shape[0]
    assert params.shape == (ncell, NPARAM), params.shape
    assert rbuf.shape[0] == ncell and mbuf.shape[0] == ncell, (rbuf.shape, mbuf.shape)
    assert set(forcing) == set(FORCING), sorted(forcing)
    zone = (np.arange(ncell, dtype=np.int64) if zone is None
            else np.ascontiguousarray(zone, dtype=np.int64))
    assert zone.shape == (ncell,), zone.shape
    if not (vegclass is None) == (cropfac_table is None) == (gamma_table is None):
        raise ValueError("vegclass, cropfac_table and gamma_table go together")
    ntime = forcing["rain"].shape[1]
    columns = (np.arange(len(core.COLUMNS), dtype=np.int64) if columns is None
               else np.ascontiguousarray(columns, dtype=np.int64))
    outdays = np.ascontiguousarray(outdays, dtype=np.int64)
    n = nrows(outdays, ntime)
    days = np.concatenate(([0], np.minimum(outdays[: n - 1], ntime))).astype(np.int64)

    if resolve(engine) == "compiled":
        from . import compiled
        values, nonconv = compiled.simulate_cells(
            params, zone=zone, rbuf=rbuf, mbuf=mbuf, **forcing, vegclass=vegclass,
            cropfac_table=cropfac_table, gamma_table=gamma_table, outdays=outdays,
            nstep=nstep, mxiter=mxiter, tol=tol, columns=columns, threads=threads)
    else:
        values, nonconv = _run_python(params, rbuf, mbuf, zone, forcing, vegclass,
                                      cropfac_table, gamma_table, outdays, nstep, mxiter,
                                      tol, n, columns)
    return CellRun(days=days, values=values, nonconverged=nonconv)


def _run_python(params, rbuf, mbuf, zone, forcing, vegclass, cropfac_table, gamma_table,
                outdays, nstep, mxiter, tol, n, columns):
    """``core.simulate`` once per cell.

    Forcing goes in as lists of Python floats: ``**`` and ``exp`` reach libm by
    different routes for Python and NumPy scalars, and the Phase 2 gate was
    measured on the Python ones.
    """
    ncell = params.shape[0]
    values = np.empty((ncell, n, columns.size))
    nonconv = np.zeros((ncell, 2), dtype=np.int64)
    for c in range(ncell):
        z = int(zone[c])
        row = {name: a[z if a.shape[0] > 1 else 0] for name, a in forcing.items()}
        if vegclass is not None:
            k = vegclass[z if vegclass.shape[0] > 1 else 0]
            row["cropfac"] = cropfac_table[c][k]
            row["gamma"] = gamma_table[c][k]
        # Each cell's own buffer, without the padding to the widest cell, so
        # that core sizes its arrays to this cell alone.  Exact either way.
        sim = core.simulate(
            **param_kwargs(params[c]),
            rbuf=np.trim_zeros(rbuf[c, 1:], "b").tolist(),
            mbuf=np.trim_zeros(mbuf[c, 1:], "b").tolist(),
            nstep=int(nstep), mxiter=int(mxiter), tol=float(tol),
            **{name: a.tolist() for name, a in row.items()}, outdays=outdays.tolist(),
        )
        values[c] = sim.values[:, columns]
        nonconv[c] = (sim.nonconverged_upper, sim.nonconverged_lower)
    return values, nonconv
