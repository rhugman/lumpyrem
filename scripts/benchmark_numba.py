"""Phase 3 benchmark: does the Numba kernel reach Fortran speed, and scale?

Three measurements on the same ten-year, five-substep workload the plan's
section 3 used:

1. **Per cell, one core** -- the compiled kernel against ``rechmod2.f`` called
   directly from a small Fortran harness, so neither side pays for file I/O or
   process start-up.  Both report a checksum, which must agree, so the two are
   known to have done the same work.
2. **Pure Python** -- ``lumpyrem.core``, for scale.
3. **Scaling** -- ``compiled.simulate_cells`` over many cells, across thread
   counts.

It also checks the correctness half of the Phase 3 gate on the way: an N-cell
run must reproduce N single runs exactly.

    python scripts/benchmark_numba.py [--years 10] [--cells 2048]

Run it on a quiet machine: the script reports the load average and warns
when other work is likely to have contaminated the timings.
"""

from __future__ import annotations

import argparse
import os
import platform
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO_ROOT / "src"), str(REPO_ROOT / "tests")]

import numba  # noqa: E402

from lumpyrem import compiled, core  # noqa: E402
from oracle.build import FFLAGS, VENDOR_DIR, gfortran  # noqa: E402

HARNESS = """\
program bench
  implicit none
  integer, parameter :: maxdelay = 500
  integer :: numdays, nout, nstep, mxiter, reps, r, imod, iday1, iday2, icall, i
  double precision :: tol, maxvol, irrigvolfrac, ks, l, m, mflowmax, rdelay, mdelay
  double precision :: maxvol_br, extravol_br, ks_br, l_br, m_br, gamma_br
  double precision :: vol0, vol_br0, vol, vol_br
  double precision :: recharge, mrecharge, runoff, evapn, rainfall, potevapn
  double precision :: irrigation, gwithdrawal, drainage_br, overflow_br
  double precision :: evapn_br, potevapn_br, total, best, elapsed
  double precision :: drainsub(maxdelay), macsub(maxdelay)
  double precision, allocatable :: rain(:), epot(:), cropfac(:), gamma(:)
  double precision, allocatable :: gwfrac(:), epot_br(:)
  integer, allocatable :: irrigcode(:), outdays(:)
  integer(8) :: t0, t1, rate

  read(*,*) numdays, nout, nstep, mxiter, tol, reps
  read(*,*) maxvol, irrigvolfrac, ks, l, m, mflowmax, rdelay, mdelay
  read(*,*) maxvol_br, extravol_br, ks_br, l_br, m_br, gamma_br, vol0, vol_br0
  allocate(rain(numdays), epot(numdays), cropfac(numdays), gamma(numdays))
  allocate(gwfrac(numdays), epot_br(numdays), irrigcode(numdays), outdays(nout))
  read(*,*) outdays
  do i = 1, numdays
    read(*,*) rain(i), epot(i), cropfac(i), gamma(i), irrigcode(i), gwfrac(i), epot_br(i)
  end do

  best = huge(best)
  do r = 1, reps
    call system_clock(t0, rate)
    vol = vol0
    vol_br = vol_br0
    drainsub = 0.0d0
    macsub = 0.0d0
    total = 0.0d0
    icall = 0
    iday1 = 0
    do imod = 1, nout
      iday2 = min(outdays(imod), numdays)
      icall = icall + 1
      call rechmod(icall, mxiter, tol, iday2 - iday1, maxdelay, nstep, maxvol, &
        irrigvolfrac, cropfac(iday1 + 1), gamma(iday1 + 1), ks, l, m, mflowmax, &
        rdelay, mdelay, vol, drainsub, macsub, rain(iday1 + 1), epot(iday1 + 1), &
        recharge, mrecharge, runoff, evapn, rainfall, potevapn, irrigation, &
        gwithdrawal, irrigcode(iday1 + 1), gwfrac(iday1 + 1), maxvol_br, &
        extravol_br, ks_br, l_br, m_br, vol_br, drainage_br, overflow_br, &
        epot_br(iday1 + 1), gamma_br, evapn_br, potevapn_br)
      total = total + recharge + mrecharge + drainage_br + overflow_br
      iday1 = iday2
      if (iday2 >= numdays) exit
    end do
    call system_clock(t1)
    elapsed = dble(t1 - t0) / dble(rate)
    best = min(best, elapsed)
  end do
  write(*, '(3es26.17)') best, vol, total
end program bench
"""


def workload(years: int, two_store: bool, seed: int = 20260928) -> dict:
    """A plausible site: seasonal evaporation, storm rainfall, summer irrigation."""
    rng = np.random.default_rng(seed)
    n = 365 * years
    t = np.arange(n)
    season = 0.5 * (1.0 + np.sin(2 * np.pi * t / 365.0))
    wet = rng.random(n) < 0.3
    spec = dict(
        maxvol=0.5, irrigvolfrac=0.3, rdelay=12.5, mdelay=1.5,
        ks=0.9, m=0.35, l=0.55, mflowmax=0.1,
        maxvol_br=0.0, extravol_br=0.0, gamma_br=0.0,
        ks_br=0.0, m_br=0.0, l_br=0.0,
        offset=15.0, factor1=1.0, factor2=1.11, power=2.0, datum=120.0,
        vol=0.15, vol_br=0.0,
        nstep=5, mxiter=500, tol=1e-5,
        rain=np.where(wet, rng.exponential(0.012, n), 0.0),
        epot=0.001 + 0.005 * season,
        cropfac=0.4 + 0.5 * season,
        gamma=np.full(n, 4.5),
        irrigcode=((season > 0.8) & (rng.random(n) < 0.9)).astype(np.int64),
        gwirrigfrac=np.full(n, 0.6),
        epot_br=np.zeros(n),
        outdays=np.arange(30, n + 30, 30).clip(max=n),
    )
    if two_store:
        spec.update(maxvol_br=0.8, extravol_br=0.2, gamma_br=2.0,
                    ks_br=0.05, m_br=0.4, l_br=0.8, vol_br=0.3,
                    epot_br=0.0005 + 0.001 * season)
    spec["outdays"] = np.unique(spec["outdays"])
    return spec


PARAM_KEYS = ("maxvol", "irrigvolfrac", "rdelay", "mdelay", "ks", "m", "l", "mflowmax",
              "maxvol_br", "extravol_br", "gamma_br", "ks_br", "m_br", "l_br",
              "offset", "factor1", "factor2", "power", "datum", "vol", "vol_br")
FORCING_KEYS = ("rain", "epot", "cropfac", "gamma", "irrigcode", "gwirrigfrac", "epot_br")


def checksum(values: np.ndarray) -> float:
    """Same accumulation, in the same order, as the Fortran harness."""
    total = 0.0
    for row in values[1:]:
        total = total + row[10] + row[11] + row[12] + row[13]
    return total


def build_fortran(workdir: Path) -> Path:
    src = workdir / "bench.f90"
    src.write_text(HARNESS)
    exe = workdir / "bench"
    subprocess.run([gfortran(), *FFLAGS, "-o", str(exe), str(src),
                    str(VENDOR_DIR / "rechmod2.f")], check=True, cwd=workdir,
                   capture_output=True)
    return exe


def time_fortran(exe: Path, spec: dict, reps: int) -> tuple[float, float, float]:
    n, out = len(spec["rain"]), spec["outdays"]
    lines = [
        f"{n} {out.size} {spec['nstep']} {spec['mxiter']} {spec['tol']!r} {reps}",
        " ".join(repr(float(spec[k])) for k in
                 ("maxvol", "irrigvolfrac", "ks", "l", "m", "mflowmax", "rdelay", "mdelay")),
        " ".join(repr(float(spec[k])) for k in
                 ("maxvol_br", "extravol_br", "ks_br", "l_br", "m_br", "gamma_br",
                  "vol", "vol_br")),
        " ".join(str(int(d)) for d in out),
    ]
    for i in range(n):
        lines.append(" ".join([
            repr(float(spec["rain"][i])), repr(float(spec["epot"][i])),
            repr(float(spec["cropfac"][i])), repr(float(spec["gamma"][i])),
            str(int(spec["irrigcode"][i])), repr(float(spec["gwirrigfrac"][i])),
            repr(float(spec["epot_br"][i])),
        ]))
    proc = subprocess.run([str(exe)], input="\n".join(lines) + "\n",
                          capture_output=True, text=True, check=True)
    if "itn limit" in proc.stdout:
        raise RuntimeError("benchmark workload does not converge in the reference")
    best, vol, total = (float(x) for x in proc.stdout.split()[-3:])
    return best, vol, total


def time_numba(spec: dict, reps: int) -> tuple[float, np.ndarray]:
    p = compiled.param_row(**{k: spec[k] for k in PARAM_KEYS})
    forcing = compiled._forcing(*(spec[k] for k in FORCING_KEYS))
    outdays = np.ascontiguousarray(spec["outdays"], dtype=np.int64)
    subdim = core.MAXDELAY
    rbuf = np.zeros(subdim + 1)
    mbuf = np.zeros(subdim + 1)
    values = np.empty((compiled._nrows(outdays, forcing[0].size), compiled.NCOLUMN))
    nonconv = np.zeros(2, dtype=np.int64)
    args = (p, rbuf, mbuf, *forcing, outdays, spec["nstep"], spec["mxiter"],
            spec["tol"], subdim, values, nonconv)
    compiled._simulate_cell(*args)          # compile, or load from cache
    best = float("inf")
    for _ in range(reps):
        t0 = time.perf_counter()
        compiled._simulate_cell(*args)
        best = min(best, time.perf_counter() - t0)
    assert nonconv.sum() == 0, nonconv
    return best, values


def time_python(spec: dict) -> float:
    t0 = time.perf_counter()
    core.simulate(bucket="upper", **{k: spec[k] for k in PARAM_KEYS + FORCING_KEYS},
                  nstep=spec["nstep"], mxiter=spec["mxiter"], tol=spec["tol"],
                  outdays=spec["outdays"].tolist())
    return time.perf_counter() - t0


def grid(spec: dict, ncell: int, seed: int = 7) -> np.ndarray:
    """``ncell`` parameter rows, each a +-20% perturbation of the site."""
    rng = np.random.default_rng(seed)
    base = compiled.param_row(**{k: spec[k] for k in PARAM_KEYS})
    rows = np.tile(base, (ncell, 1))
    for idx in (compiled.MAXVOL, compiled.KS, compiled.M, compiled.L, compiled.RDELAY):
        rows[:, idx] *= rng.uniform(0.8, 1.2, ncell)
    rows[:, compiled.VOL] = np.minimum(rows[:, compiled.VOL], rows[:, compiled.MAXVOL])
    return rows


def run_grid(spec: dict, params: np.ndarray):
    return compiled.simulate_cells(
        params, **{k: spec[k] for k in FORCING_KEYS},
        nstep=spec["nstep"], mxiter=spec["mxiter"], tol=spec["tol"],
        outdays=spec["outdays"])


def check_grid_equals_singles(spec: dict) -> None:
    """Phase 3 gate, correctness half: N cells == N single runs, exactly."""
    base = compiled.param_row(**{k: spec[k] for k in PARAM_KEYS})
    same, _ = run_grid(spec, np.tile(base, (64, 1)))
    varied_params = grid(spec, 64)
    varied, _ = run_grid(spec, varied_params)
    single, _ = run_grid(spec, base[None, :])
    assert all(np.array_equal(same[c], single[0]) for c in range(64))
    for c in range(64):
        one, _ = run_grid(spec, varied_params[c:c + 1])
        assert np.array_equal(varied[c], one[0]), c


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--years", type=int, default=10)
    ap.add_argument("--cells", type=int, default=2048)
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--rounds", type=int, default=6)
    args = ap.parse_args()

    ncpu = os.cpu_count() or 1
    print(f"{platform.platform()}  {platform.processor() or platform.machine()}  "
          f"{ncpu} logical CPUs, numba {numba.__version__}, numpy {np.__version__}")
    print(f"workload: {args.years} years x 5 substeps/day, output every 30 days\n")

    load = os.getloadavg()[0]
    if load > ncpu / 4:
        print(f"WARNING: load average {load:.0f} on {ncpu} CPUs -- absolute timings and "
              "scaling are contaminated by other work.\nFortran and Numba are timed in "
              "alternating rounds, so their ratio is the trustworthy number.\n")

    with tempfile.TemporaryDirectory() as tmp:
        exe = build_fortran(Path(tmp))
        for label, two_store in (("one store", False), ("two stores", True)):
            spec = workload(args.years, two_store)
            # Alternate the two, a few reps at a time, and keep each one's best:
            # whatever else the machine is doing lands on both sides equally.
            f_time = n_time = float("inf")
            for _ in range(args.rounds):
                f_best, f_vol, f_total = time_fortran(exe, spec, args.reps)
                n_best, values = time_numba(spec, args.reps)
                f_time, n_time = min(f_time, f_best), min(n_time, n_best)
            n_vol, n_total = float(values[-1, 0]), checksum(values)
            same = (n_vol == f_vol) and (n_total == f_total)
            print(f"[{label}]")
            print(f"  Fortran rechmod, -O2        {f_time * 1e3:8.2f} ms / cell")
            print(f"  Numba, 1 thread             {n_time * 1e3:8.2f} ms / cell   "
                  f"{n_time / f_time:5.2f}x Fortran")
            print(f"  checksum agrees with Fortran: {'yes, exactly' if same else 'NO'}"
                  + ("" if same else f"  (vol {n_vol!r} vs {f_vol!r}, "
                                     f"total {n_total!r} vs {f_total!r})"))
            if not two_store:
                p_time = time_python(spec)
                print(f"  pure Python (lumpyrem.core) {p_time * 1e3:8.0f} ms / cell   "
                      f"{p_time / f_time:5.0f}x Fortran")
            print()

    spec = workload(args.years, two_store=False)
    check_grid_equals_singles(spec)
    print("N-cell run reproduces N single runs exactly: yes (64 identical, 64 varied)\n")

    params = grid(spec, args.cells)
    run_grid(spec, params[:8])                 # compile the parallel entry point
    threads = sorted({t for t in (1, 2, 4, 6, 8, 10, 12, 16, 24, 32, ncpu) if t <= ncpu})
    print(f"scaling, {args.cells} cells (+-20% parameter spread):")
    print("  threads   wall s   ms/cell   speedup   efficiency")
    base_wall = None
    for t in threads:
        numba.set_num_threads(t)
        best = float("inf")
        for _ in range(3):
            t0 = time.perf_counter()
            run_grid(spec, params)
            best = min(best, time.perf_counter() - t0)
        base_wall = base_wall or best
        speedup = base_wall / best
        print(f"  {t:7d}  {best:7.3f}  {best / args.cells * 1e3:8.3f}  "
              f"{speedup:7.2f}x  {speedup / t:10.0%}")


if __name__ == "__main__":
    main()
