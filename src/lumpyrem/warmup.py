"""Compile the Numba kernels ahead of time, so a timed or calibration run does not.

Numba compiles a kernel the first time it is called, which costs far more than
a run of a small grid.  lumpyrem's kernels are declared with ``cache=True``, so
the machine code is saved to disk and later processes load it in a fraction of
a second.  Warm the cache once -- after installing or updating lumpyrem, or
after changing Python, Numba or NumPy -- with::

    python -m lumpyrem.warmup              # or the installed lumpyrem-compile

or from Python with :func:`compile_kernels`.  Every entry point a run uses is
called on a tiny synthetic grid, with the argument types a real run passes, so
every specialisation is in the cache: zoned forcing with and without a
vegetation-class calendar, shared and per-cell forcing, a single model, and the
parallel and single-threaded loops.  Keeping this next to the kernels is what
keeps it honest: a kernel whose signature changes is exercised here in the
same commit, and ``tests/test_warmup.py`` checks it.

The cache lives in ``lumpyrem/__pycache__`` (``*.nbi``, ``*.nbc``) unless
``NUMBA_CACHE_DIR`` is set.  It is keyed on the path and modification time of
``compiled.py``, so a run from another copy of the package -- a PEST++ worker
directory holding its own lumpyrem, say -- compiles again unless every copy
points ``NUMBA_CACHE_DIR`` at one shared cache warmed with ``--cache-dir``.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

__all__ = ["WarmupReport", "cache_dir", "compile_kernels", "main"]

#: A fresh process slower than this to make its first calls is taken to have
#: compiled rather than loaded from the cache.
RECOMPILE_SECONDS = 5.0

_PACKAGE_PARENT = str(Path(__file__).resolve().parents[1])


@dataclass(frozen=True)
class WarmupReport:
    """What :func:`compile_kernels` did, in seconds.

    Attributes:
        first: the first calls in this process -- a compile, or a load from
            the cache if it was already warm.
        second: the same calls again, compiled and in memory.
        fresh: the first calls in a new Python process, which sees only what
            is on disk; ``None`` unless asked for.
        cache: the cache directory.
        threads: threads Numba runs a parallel grid on.
    """

    first: float
    second: float
    fresh: float | None
    cache: str
    threads: int

    @property
    def cached(self) -> bool | None:
        """Whether a new process loads the kernels rather than compiling them."""
        return None if self.fresh is None else self.fresh <= RECOMPILE_SECONDS

    def __str__(self) -> str:
        lines = [f"cache: {self.cache}  ({self.threads} threads)",
                 f"first calls (compile, or load from cache): {self.first:6.2f} s",
                 f"second calls (in memory):                  {self.second:6.2f} s"]
        if self.fresh is not None:
            lines.append(f"fresh process, first calls:                {self.fresh:6.2f} s  "
                         + ("(loaded from cache)" if self.cached else
                            "(slow enough to be a recompile: is the cache writable, and "
                            "is NUMBA_CACHE_DIR the same for every run?)"))
        return "\n".join(lines)


def cache_dir() -> str:
    """Where the compiled kernels are cached."""
    env = os.environ.get("NUMBA_CACHE_DIR")
    return env if env else str(Path(__file__).resolve().parent / "__pycache__")


def _exercise() -> None:
    """Call every kernel entry point once, with the types a real run passes."""
    import numpy as np
    import pandas as pd

    import lumpyrem as lr

    days = pd.date_range("2000-01-01", "2000-03-31")
    rng = np.random.default_rng(0)
    rain = pd.DataFrame(rng.gamma(0.5, 6.0, (len(days), 2)), index=days, columns=[1, 2])
    pet = pd.Series(5.0, index=days)
    veg = pd.DataFrame([[0, 1], [1, 0]], index=pd.to_datetime(["2000-01-01", "2000-02-01"]),
                       columns=[1, 2])
    grid = lr.ModelGrid.from_arrays(
        upper=dict(maxvol=[150.0, 170.0, 190.0, 210.0], ks=0.6, m=0.4, l=0.5,
                   mflowmax=1.6, rdelay=10.4, mdelay=2.9, irrigvolfrac=0.5),
        initial=dict(vol=0.5), solver=lr.Solver(nstep=12, mxiter=200, tol=1e-5))
    ncell = grid.ncell

    # Zoned forcing with a vegetation-class calendar, on every thread and on one
    zoned = lr.ZonedForcing.from_frames(
        [1, 2, 2, 1], rainfall=rain, pot_evap=pet, veg_class=veg,
        crop_factor=np.full((ncell, 2), 0.8), veg_gamma=np.full((ncell, 2), 2.0))
    grid.run(zoned, times="MS", columns=["net_recharge"])
    grid.run(zoned, times="MS", columns=["net_recharge"], single_threaded=True)
    # Zoned forcing without a calendar: the same kernels, the other branch
    grid.run(lr.ZonedForcing.from_frames([1, 2, 2, 1], rainfall=rain, pot_evap=pet,
                                         veg_gamma=2.0), times="MS")
    # Shared and per-cell forcing
    shared = lr.Forcing.from_dataframe(pd.DataFrame({"rainfall": rain[1], "pot_evap": pet}),
                                       veg_gamma=2.0)
    grid.run(shared, times="MS")
    per_cell = lr.GridForcing.from_arrays(days, rainfall=rain[[1, 2, 1, 2]].to_numpy(),
                                          pot_evap=pet.to_numpy(), veg_gamma=2.0)
    grid.run(per_cell, times="MS")
    grid.run(per_cell, times="MS", single_threaded=True)
    # A single model
    grid[0].run(shared, times="MS")


def _fresh_process_seconds() -> float:
    """Time the first calls in a new interpreter, importing this same copy."""
    code = ("import time; t = time.perf_counter(); "
            "from lumpyrem.warmup import _exercise; _exercise(); "
            "print(time.perf_counter() - t)")
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(filter(None, [_PACKAGE_PARENT,
                                                      env.get("PYTHONPATH")]))
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                         env=env, cwd=_PACKAGE_PARENT)
    if out.returncode != 0:
        raise RuntimeError(f"the fresh-process check failed:\n{out.stderr}")
    return float(out.stdout.strip().splitlines()[-1])


def compile_kernels(*, check: bool = False, cache: str | os.PathLike | None = None
                    ) -> WarmupReport:
    """Compile every kernel a run uses and save it to the Numba cache.

    Args:
        check: also start a new Python process and time its first calls, to
            confirm the kernels now load from disk rather than compile.
        cache: a cache directory to use, as ``NUMBA_CACHE_DIR``.  Numba reads
            that variable when it is first imported, so this works only
            before anything has imported Numba; otherwise set the variable
            before starting Python.

    Returns:
        A :class:`WarmupReport`.

    Raises:
        RuntimeError: Numba is not installed -- lumpyrem then runs its
            pure-Python kernel and there is nothing to compile -- or ``cache``
            came too late to take effect.
    """
    if cache is not None:
        cache = str(Path(cache).resolve())
        if "numba" in sys.modules and os.environ.get("NUMBA_CACHE_DIR") != cache:
            raise RuntimeError("Numba is already imported, so the cache directory can no "
                               "longer change; set NUMBA_CACHE_DIR before starting Python")
        os.makedirs(cache, exist_ok=True)
        os.environ["NUMBA_CACHE_DIR"] = cache

    from . import engine
    if not engine.compiled_available():
        raise RuntimeError("Numba is not importable, so lumpyrem runs its pure-Python "
                           "kernel and there is nothing to compile; install numba "
                           "(pip install lumpyrem[fast]) for the compiled kernel")
    import numba

    t0 = time.perf_counter()
    _exercise()
    first = time.perf_counter() - t0
    t0 = time.perf_counter()
    _exercise()
    second = time.perf_counter() - t0
    fresh = _fresh_process_seconds() if check else None
    return WarmupReport(first=first, second=second, fresh=fresh, cache=cache_dir(),
                        threads=int(numba.config.NUMBA_NUM_THREADS))


def main(argv=None) -> int:
    """``python -m lumpyrem.warmup`` / ``lumpyrem-compile``.

    Exit status 0 when the kernels are cached, 1 when Numba is missing, 2 when
    a fresh process still appears to compile.
    """
    parser = argparse.ArgumentParser(
        prog="lumpyrem-compile",
        description="Compile lumpyrem's Numba kernels and save them to the cache.")
    parser.add_argument("--cache-dir", metavar="DIR",
                        help="cache directory (sets NUMBA_CACHE_DIR); point every "
                             "PEST++ worker at the same one")
    parser.add_argument("--no-check", action="store_true",
                        help="skip timing a fresh process against the cache")
    args = parser.parse_args(argv)

    from . import __version__
    print(f"lumpyrem {__version__} from {Path(__file__).resolve().parent}")
    try:
        report = compile_kernels(check=not args.no_check, cache=args.cache_dir)
    except RuntimeError as exc:
        print(exc, file=sys.stderr)
        return 1
    print(report)
    return 2 if report.cached is False else 0


if __name__ == "__main__":
    sys.exit(main())
