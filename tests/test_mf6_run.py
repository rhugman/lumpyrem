"""lumpyrem.mf6 end to end: run MODFLOW 6 and read the recharge back out.

The unit tests in ``test_mf6.py`` check the RCH file's text.  These run the
file through MODFLOW 6 itself, which is the only way to know the model reads
the package at all (the name file has to list it), reads each row as the
cell and rate intended, and reads the boundnames as written.

The model is two layers of six 100 m cells, with one top cell inactive so
one column's recharge lands on layer 2.  A constant head in layer 2 under an
active top cell takes the water out; recharge is never applied to it, since
MODFLOW reports no boundary flow into a constant-head cell.  The periods are
a one-day steady state averaged over a year of recharge, then three calendar
months, as in ``test_mf6.py``.

MODFLOW 6 is found through ``$MF6_EXE``, else the one ``get-modflow :flopy``
installs, else ``mf6`` on the PATH; the tests skip without it.  flopy's own
folder comes before the PATH because an older ``mf6`` bundled with other
software is often on the PATH.
"""

from __future__ import annotations

import importlib
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from lumpyrem import ModelGrid, Solver, ZonedForcing
from lumpyrem import mf6 as lmf6

flopy = pytest.importorskip("flopy")

pytestmark = pytest.mark.mf6

DX = 100.0
NROW, NCOL = 2, 3
NCPL = NROW * NCOL
#: One zone label per column; two contain spaces, as farm and district names do.
ZONES = ["farm a", "farm a", "river flat", "river flat", "ridge", "farm a"]
OBS_NAME = "farm a"


@pytest.fixture(scope="module")
def mf6_exe():
    exe = os.environ.get("MF6_EXE")
    if not exe:
        bindir = Path(importlib.import_module("flopy.utils.get_modflow")
                      .flopy_appdata_path) / "bin"
        exe = shutil.which("mf6", path=str(bindir)) or shutil.which("mf6")
    if not exe:
        pytest.skip("MODFLOW 6 not found; run get-modflow :flopy, set MF6_EXE or put "
                    "mf6 on the PATH")
    return exe


def _build(ws, exe):
    """The simulation without RCH: TDIS, IMS, DISV, IC, NPF, STO, CHD, OC."""
    sim = flopy.mf6.MFSimulation(sim_name="t", sim_ws=str(ws), exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=4, time_units="days", start_date_time="2000-12-31",
                          perioddata=[(1.0, 1, 1.0), (31.0, 1, 1.0), (28.0, 1, 1.0),
                                      (31.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, complexity="simple", outer_dvclose=1e-8,
                         inner_dvclose=1e-8)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m")

    vertices = [(i * (NCOL + 1) + j, j * DX, i * DX)
                for i in range(NROW + 1) for j in range(NCOL + 1)]
    cell2d = []
    for i in range(NROW):
        for j in range(NCOL):
            v = i * (NCOL + 1) + j
            cell2d.append((i * NCOL + j, (j + 0.5) * DX, (i + 0.5) * DX, 4,
                           v, v + NCOL + 1, v + NCOL + 2, v + 1))
    idomain = np.ones((2, NCPL), dtype=int)
    idomain[0, 4] = 0
    flopy.mf6.ModflowGwfdisv(gwf, length_units="meters", nlay=2, ncpl=NCPL,
                             nvert=len(vertices), top=10.0,
                             botm=[[5.0] * NCPL, [0.0] * NCPL], vertices=vertices,
                             cell2d=cell2d, idomain=idomain)
    flopy.mf6.ModflowGwfic(gwf, strt=8.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=10.0, k33=1.0, save_flows=True)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=0, ss=1e-5, sy=0.1,
                            steady_state={0: True}, transient={1: True})
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[((1, 0), 8.0)])
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord="m.cbc", head_filerecord="m.hds",
                           saverecord=[("BUDGET", "ALL"), ("HEAD", "ALL")])
    return sim, gwf


def _rates(sim, gwf):
    """Per-period rates from a lumpyrem run, distinct for every cell."""
    periods = lmf6.stress_periods(sim.tdis, gwf.get_package("sto"),
                                  steady_state_window=("2000-01-01", "2001-01-01"))
    cells = lmf6.top_active_cells(gwf.get_package("disv"))
    days = pd.date_range("2000-01-01", "2001-03-31")
    season = 1.0 + 0.8 * np.sin(2 * np.pi * days.dayofyear / 365.25)
    rain = pd.DataFrame({"farm a": 1.2 * season, "river flat": 1.0 * season,
                         "ridge": 1.4 * season}, index=days)
    forcing = ZonedForcing.from_frames(ZONES, rainfall=rain,
                                       pot_evap=pd.Series(1.0, index=days), veg_gamma=2.0)
    grid = ModelGrid.from_arrays(
        upper=dict(maxvol=50.0, ks=np.linspace(2.0, 8.0, NCPL), m=1.0, l=0.5),
        initial=dict(vol=10.0), solver=Solver(nstep=4), ncell=NCPL)
    results = grid.run(forcing, times=lmf6.output_times(periods, days[0]), start=days[0],
                       end=days[-1], columns=["net_recharge"])
    rates = lmf6.period_rates(results, periods, factor=0.001)
    assert (rates > 0).all() and len(np.unique(rates)) == rates.size
    return cells, rates


def _cell_areas(gwf, icpl):
    """Plan area of each cell, by the shoelace formula over its CELL2D vertices."""
    disv = gwf.get_package("disv")
    xy = {int(v["iv"]): (v["xv"], v["yv"]) for v in disv.vertices.array}
    rows = {int(r[0]): [int(i) for i in list(r)[4:4 + int(r[3])]] for r in disv.cell2d.array}
    areas = []
    for c in icpl:
        x, y = np.array([xy[i] for i in rows[int(c)]]).T
        areas.append(0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))
    return np.array(areas)


@pytest.mark.parametrize("method", ["fast", "flopy", "flopy-now", "flopy-obs-file",
                                    "fast-loaded", "flopy-loaded"])
def test_modflow_budget_matches_rates(tmp_path, mf6_exe, method):
    sim, gwf = _build(tmp_path, mf6_exe)
    cells, rates = _rates(sim, gwf)
    names = [ZONES[c] for c in cells["icpl"]]
    obs_csv = "m.rch.obs.csv"
    run = sim

    def write_obs_file():
        with open(tmp_path / "m.rch.obs", "w") as fid:
            fid.write(f"BEGIN continuous FILEOUT {obs_csv}\n"
                      f"  rch_obs  rch  '{OBS_NAME}'\nEND continuous\n")

    if method.endswith("-loaded"):
        # A model already set up with RCH, loaded with only the packages a
        # recharge script needs -- the name file still lists the RCH file.
        lmf6.recharge_package(gwf, cells, np.zeros_like(rates), write=False)
        sim.write_simulation(silent=True)
        write_obs_file()
        run = flopy.mf6.MFSimulation.load(sim_ws=str(tmp_path), exe_name=mf6_exe,
                                          load_only=["disv", "sto"], verbosity_level=0)
        loaded = run.get_model("m")
        lmf6.recharge_package(loaded, cells, rates, boundnames=names,
                              method=method.split("-")[0], obs_file="m.rch.obs")
        nam = (tmp_path / "m.nam").read_text().upper().split()
        assert nam.count("RCH6") == 1, "the name file lists RCH more than once"
    elif method == "fast":
        # The documented order: write the simulation, then the RCH file over it.
        sim.write_simulation(silent=True)
        write_obs_file()
        lmf6.recharge_package(gwf, cells, rates, boundnames=names, method="fast",
                              obs_file="m.rch.obs")
    elif method == "flopy-obs-file":
        # A simulation on disk with its own OBS6 file, which flopy points at.
        sim.write_simulation(silent=True)
        write_obs_file()
        lmf6.recharge_package(gwf, cells, rates, boundnames=names, method="flopy",
                              obs_file="m.rch.obs")
    elif method == "flopy":
        # The documented order: build the package object, then write it all.
        lmf6.recharge_package(gwf, cells, rates, boundnames=names, method="flopy",
                              write=False,
                              observations={obs_csv: [("rch_obs", "rch", OBS_NAME)]})
        sim.write_simulation(silent=True)
    else:
        # A simulation already on disk, given its RCH package afterwards.
        sim.write_simulation(silent=True)
        lmf6.recharge_package(gwf, cells, rates, boundnames=names, method="flopy",
                              observations={obs_csv: [("rch_obs", "rch", OBS_NAME)]})

    success, buff = run.run_simulation(silent=True, report=True)
    assert success, "MODFLOW 6 failed:\n" + "\n".join(buff[-20:])

    # Per cell and period: the flow RCH put into the model is rate x area.
    expected = rates * _cell_areas(gwf, cells["icpl"])
    nodes = cells["node"].to_numpy() + 1        # the budget file's nodes are 1-based
    cbc = gwf.output.budget()
    assert any(t.strip() == b"RCH" for t in cbc.get_unique_record_names()), \
        "no RCH flows in the budget file; MODFLOW did not read the package"
    for k in range(rates.shape[0]):
        (rec,) = cbc.get_data(text="RCH", kstpkper=(0, k))
        got = dict(zip(rec["node"].tolist(), rec["q"].tolist(), strict=True))
        assert sorted(got) == sorted(nodes.tolist()), f"period {k + 1}: wrong cells"
        np.testing.assert_allclose([got[n] for n in nodes], expected[k], rtol=1e-5,
                                   err_msg=f"period {k + 1}")

    # The observation sums the cells whose boundname is OBS_NAME, which only
    # works if MODFLOW read each boundname whole.
    obs = pd.read_csv(tmp_path / obs_csv)
    mine = np.array([n == OBS_NAME for n in names])
    np.testing.assert_allclose(obs["RCH_OBS"].to_numpy(), expected[:, mine].sum(axis=1),
                               rtol=1e-5)
