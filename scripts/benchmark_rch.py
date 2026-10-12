"""RCH writers: is lumpyrem's own text writer worth keeping over flopy's?

Four ways to put per-period recharge into a MODFLOW 6 model, on one-layer
DISV grids of increasing size:

1. ``fast`` -- ``recharge_package(method="fast")``, lumpyrem's text writer.
2. ``flopy`` -- ``recharge_package(method="flopy")``, a ``ModflowGwfrch``
   written as flopy writes it, every period inline in the RCH file.
3. ``flopy-external`` -- the same package with each period's list in its
   own external text file.
4. ``flopy-binary`` -- each period in an external binary file, flopy's own
   fast path.  flopy cannot write binary lists with boundnames, so this one
   is written without them; the other three carry a boundname per cell.

Each is timed from rates in hand to files on disk.  Then, for every grid up
to ``--run-cells``, MODFLOW 6 runs each model and the recharge each cell
receives in each period is read back from the budget file and compared with
the rates given, so speed is only counted for output MODFLOW agrees with.

    python scripts/benchmark_rch.py [--grids 1000x12 10000x120 50000x550]
                                    [--run-cells 10000] [--mf6 path/to/mf6]
"""

from __future__ import annotations

import argparse
import importlib
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO_ROOT / "src")]

import flopy  # noqa: E402

from lumpyrem import mf6 as lmf6  # noqa: E402

METHODS = ("fast", "flopy", "flopy-external", "flopy-binary")
DX = 100.0


def find_mf6(given: str | None) -> str | None:
    if given:
        return given
    bindir = Path(importlib.import_module("flopy.utils.get_modflow").flopy_appdata_path)
    return (os.environ.get("MF6_EXE") or shutil.which("mf6", path=str(bindir / "bin"))
            or shutil.which("mf6"))


def build(ws: Path, ncpl: int, nper: int, exe: str | None):
    """A one-layer DISV grid of ``ncpl`` square cells, ``nper`` 30-day periods."""
    ncol = int(np.ceil(np.sqrt(ncpl)))
    nrow = int(np.ceil(ncpl / ncol))
    sim = flopy.mf6.MFSimulation(sim_name="b", sim_ws=str(ws), exe_name=exe or "mf6",
                                 verbosity_level=0)
    flopy.mf6.ModflowTdis(sim, nper=nper, time_units="days",
                          perioddata=[(30.0, 1, 1.0)] * nper)
    flopy.mf6.ModflowIms(sim, complexity="simple")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
    iv = np.arange((nrow + 1) * (ncol + 1))
    vertices = list(zip(iv.tolist(), (iv % (ncol + 1) * DX).tolist(),
                        (iv // (ncol + 1) * DX).tolist(), strict=True))
    cell2d = []
    for c in range(ncpl):
        i, j = divmod(c, ncol)
        v = i * (ncol + 1) + j
        cell2d.append((c, (j + 0.5) * DX, (i + 0.5) * DX, 4, v, v + ncol + 1, v + ncol + 2,
                       v + 1))
    flopy.mf6.ModflowGwfdisv(gwf, nlay=1, ncpl=ncpl, nvert=len(vertices), top=10.0,
                             botm=[0.0], vertices=vertices, cell2d=cell2d)
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=10.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=0, ss=1e-5, transient={0: True})
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[((0, 0), 5.0)])
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord="m.cbc",
                           saverecord=[("BUDGET", "LAST")])
    return sim, gwf


def write(method: str, gwf, cells, rates, names) -> float:
    """Seconds from rates in hand to the RCH package on disk and in the name file.

    The rest of the simulation has already been written, so only RCH counts.
    """
    t0 = time.perf_counter()
    if method == "fast":
        lmf6.recharge_package(gwf, cells, rates, boundnames=names, method="fast")
        return time.perf_counter() - t0
    if method == "flopy":
        pkg = lmf6.recharge_package(gwf, cells, rates, boundnames=names, method="flopy",
                                    write=False)
    else:
        # External files, one per period: text, or binary -- which flopy cannot
        # write with boundnames, so the binary variant goes without them.
        binary = method == "flopy-binary"
        cellid = [(0, c) for c in cells[:, 1].tolist()]
        spd = {}
        for k in range(rates.shape[0]):
            rows = (zip(cellid, rates[k].tolist(), strict=True) if binary else
                    zip(cellid, rates[k].tolist(), names, strict=True))
            spd[k] = {"filename": f"m.rch.{k + 1}{'.bin' if binary else '.txt'}",
                      "binary": binary, "data": list(rows)}
        pkg = flopy.mf6.ModflowGwfrch(gwf, pname="rch", filename="m.rch",
                                      boundnames=None if binary else True, save_flows=True,
                                      maxbound=len(cellid), stress_period_data=spd)
    pkg.write()
    gwf.name_file.write()
    return time.perf_counter() - t0


def budget_rch(gwf, ncell: int, nper: int) -> np.ndarray:
    cbc = gwf.output.budget()
    out = np.zeros((nper, ncell))
    for k in range(nper):
        (rec,) = cbc.get_data(text="RCH", kstpkper=(0, k))
        out[k, rec["node"] - 1] = rec["q"]
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--grids", nargs="+", default=["1000x12", "10000x120", "50000x550"],
                        help="NCELLxNPER for each grid")
    parser.add_argument("--run-cells", type=int, default=10000,
                        help="run MODFLOW and check the budget for grids up to this size")
    parser.add_argument("--mf6", help="the MODFLOW 6 executable")
    args = parser.parse_args(argv)
    exe = find_mf6(args.mf6)
    print(f"flopy {flopy.__version__}; MODFLOW 6: {exe or 'not found, budgets not checked'}")
    rng = np.random.default_rng(0)

    print(f"\n{'grid':>12} {'method':>13} {'write s':>9} {'x fast':>7} {'max rel err':>12}")
    for spec in args.grids:
        ncell, nper = (int(x) for x in spec.lower().split("x"))
        cells = np.column_stack([np.zeros(ncell, dtype=np.int64), np.arange(ncell)])
        rates = rng.uniform(1e-5, 5e-3, size=(nper, ncell))
        names = [f"zone {c % 17}" for c in range(ncell)]
        times = {}
        for method in METHODS:
            with tempfile.TemporaryDirectory() as tmp:
                sim, gwf = build(Path(tmp), ncell, nper, exe)
                sim.write_simulation(silent=True)
                times[method] = write(method, gwf, cells, rates, names)
                err = ""
                if exe and ncell <= args.run_cells:
                    ok, _ = sim.run_simulation(silent=True)
                    if not ok:
                        err = "MODFLOW FAILED"
                    else:
                        q = budget_rch(gwf, ncell, nper)
                        # Cell 0 is the constant head, which takes no RCH flow.
                        rel = q[:, 1:] / (rates[:, 1:] * DX * DX) - 1.0
                        err = f"{np.max(np.abs(rel)):.1e}"
            print(f"{spec:>12} {method:>13} {times[method]:9.2f} "
                  f"{times[method] / times['fast']:7.1f} {err:>12}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
