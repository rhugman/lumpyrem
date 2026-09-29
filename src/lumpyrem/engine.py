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

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import cache

import numpy as np

from . import core

__all__ = ["ENGINES", "CellRun", "compiled_available", "resolve", "run_cells"]

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


@dataclass(frozen=True)
class CellRun:
    """What the kernel returns for a set of cells sharing one output schedule."""

    days: np.ndarray            # int64, (nrows,), 0 first
    values: np.ndarray          # float64, (ncell, nrows, 26)
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
    engine: str = "auto",
) -> CellRun:
    """Run ``params.shape[0]`` independent cells.

    Args:
        params: ``(ncell, NPARAM)`` rows from :func:`param_row`.
        rbuf, mbuf: ``(ncell, width + 1)`` initial buffers from
            :func:`buffer_rows`.
        forcing: every name in :data:`FORCING`, each ``(1, ntime)`` or
            ``(ncell, ntime)``.
        outdays: 1-based output days, shared by every cell.
        nstep, mxiter, tol: solver settings, shared by every cell.
        engine: ``"auto"``, ``"compiled"`` or ``"python"``.
    """
    ncell = params.shape[0]
    assert params.shape == (ncell, NPARAM), params.shape
    assert rbuf.shape[0] == ncell and mbuf.shape[0] == ncell, (rbuf.shape, mbuf.shape)
    assert set(forcing) == set(FORCING), sorted(forcing)
    ntime = forcing["rain"].shape[1]
    for name in FORCING:
        shape = forcing[name].shape
        assert shape in ((1, ntime), (ncell, ntime)), (name, shape)
    outdays = np.ascontiguousarray(outdays, dtype=np.int64)
    n = nrows(outdays, ntime)
    days = np.concatenate(([0], np.minimum(outdays[: n - 1], ntime))).astype(np.int64)

    if resolve(engine) == "compiled":
        from . import compiled
        values, nonconv = compiled.simulate_cells(
            params, rbuf=rbuf, mbuf=mbuf, **forcing, outdays=outdays,
            nstep=nstep, mxiter=mxiter, tol=tol)
    else:
        values, nonconv = _run_python(params, rbuf, mbuf, forcing, outdays,
                                      nstep, mxiter, tol, n)
    return CellRun(days=days, values=values, nonconverged=nonconv)


def _run_python(params, rbuf, mbuf, forcing, outdays, nstep, mxiter, tol, n):
    """``core.simulate`` once per cell.

    Forcing goes in as lists of Python floats: ``**`` and ``exp`` reach libm by
    different routes for Python and NumPy scalars, and the Phase 2 gate was
    measured on the Python ones.
    """
    ncell = params.shape[0]
    values = np.empty((ncell, n, len(core.COLUMNS)))
    nonconv = np.zeros((ncell, 2), dtype=np.int64)
    for c in range(ncell):
        row = {name: forcing[name][c if forcing[name].shape[0] > 1 else 0].tolist()
               for name in FORCING}
        # Each cell's own buffer, without the padding to the widest cell, so
        # that core sizes its arrays to this cell alone.  Exact either way.
        sim = core.simulate(
            **param_kwargs(params[c]),
            rbuf=np.trim_zeros(rbuf[c, 1:], "b").tolist(),
            mbuf=np.trim_zeros(mbuf[c, 1:], "b").tolist(),
            nstep=int(nstep), mxiter=int(mxiter), tol=float(tol),
            **row, outdays=outdays.tolist(),
        )
        values[c] = sim.values
        nonconv[c] = (sim.nonconverged_upper, sim.nonconverged_lower)
    return values, nonconv
