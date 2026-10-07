"""lumpyrem.mf6: DISV cells, TDIS periods, per-period rates and the RCH file."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from lumpyrem import ModelGrid, Solver, ZonedForcing
from lumpyrem import mf6 as lmf6

flopy = pytest.importorskip("flopy")


@pytest.fixture
def sim(tmp_path):
    """Three columns of two layers; column 1's top layer is inactive.

    Period 1 is a one-day steady state, then three calendar months.
    """
    sim = flopy.mf6.MFSimulation(sim_name="t", sim_ws=str(tmp_path))
    flopy.mf6.ModflowTdis(sim, nper=4, time_units="days", start_date_time="2000-12-31",
                          perioddata=[(1.0, 1, 1.0), (31.0, 1, 1.0), (28.0, 1, 1.0),
                                      (31.0, 1, 1.0)])
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
    vertices = [(i, float(x), float(y)) for i, (x, y) in
                enumerate([(0, 0), (1, 0), (2, 0), (3, 0), (0, 1), (1, 1), (2, 1), (3, 1)])]
    cell2d = [(j, j + 0.5, 0.5, 4, j, j + 1, j + 5, j + 4) for j in range(3)]
    flopy.mf6.ModflowGwfdisv(gwf, nlay=2, ncpl=3, nvert=8, top=10.0,
                             botm=[[5.0] * 3, [0.0] * 3], vertices=vertices, cell2d=cell2d,
                             idomain=[[1, 0, 1], [1, 1, 1]])
    sto = flopy.mf6.ModflowGwfsto(gwf, steady_state={0: True}, transient={1: True})
    return sim, gwf, sto


def test_top_active_cells(sim):
    _, gwf, _ = sim
    top = lmf6.top_active_cells(gwf.get_package("disv"))
    assert top["layer"].tolist() == [0, 1, 0]
    assert top["icpl"].tolist() == [0, 1, 2]
    assert top["node"].tolist() == [0, 4, 2]
    assert top["x"].tolist() == [0.5, 1.5, 2.5]
    assert top["top"].tolist() == [10.0, 5.0, 10.0]


def test_stress_periods(sim):
    s, _, sto = sim
    p = lmf6.stress_periods(s.tdis, sto, steady_state_window=("2000-01-01", "2001-01-01"))
    assert p["steady"].tolist() == [True, False, False, False]
    assert p["start"].tolist() == list(pd.to_datetime(
        ["2000-12-31", "2001-01-01", "2001-02-01", "2001-03-01"]))
    assert p["perlen"].tolist() == [1, 31, 28, 31]
    assert p.loc[0, "rate_start"] == pd.Timestamp("2000-01-01")
    assert (p.loc[1:, "rate_start"] == p.loc[1:, "start"]).all()
    times = lmf6.output_times(p, "2000-01-01")
    # The steady state averages over 2000; it needs no row of its own at its
    # own one-day start.
    assert list(times) == list(pd.to_datetime(
        ["2001-01-01", "2001-02-01", "2001-03-01", "2001-04-01"]))


def test_period_rates_and_rch(sim, tmp_path):
    s, gwf, sto = sim
    p = lmf6.stress_periods(s.tdis, sto, steady_state_window=("2000-01-01", "2001-01-01"))
    days = pd.date_range("2000-01-01", "2001-03-31")
    rain = pd.DataFrame({1: 5.0, 2: 0.0}, index=days)
    cells = lmf6.top_active_cells(gwf.get_package("disv"))
    forcing = ZonedForcing.from_frames([1, 1, 2], rainfall=rain,
                                       pot_evap=pd.Series(1.0, index=days), veg_gamma=2.0)
    grid = ModelGrid.from_arrays(upper=dict(maxvol=50.0, ks=2.0, m=0.5, l=0.5),
                                 initial=dict(vol=10.0), solver=Solver(nstep=4), ncell=3)
    res = grid.run(forcing, times=lmf6.output_times(p, days[0]), start=days[0],
                   end=days[-1], columns=["net_recharge"])
    rates = lmf6.period_rates(res, p, factor=0.001)
    assert rates.shape == (4, 3)

    daily = grid.run(forcing, start=days[0], end=days[-1],
                     columns=["net_recharge"])["net_recharge"]
    feb = daily.loc["2001-02-02":"2001-03-01"].sum() / 28 * 0.001
    np.testing.assert_allclose(rates[2], feb.to_numpy(), rtol=1e-12)
    ss = daily.loc["2000-01-02":"2001-01-01"].sum() / 366 * 0.001
    np.testing.assert_allclose(rates[0], ss.to_numpy(), rtol=1e-12)

    path = lmf6.recharge_package(gwf, cells, rates, boundnames=["z1", "z1", "z2"],
                                 filename="m.rch", obs_file="m.rch.obs")
    text = open(path).read()
    assert "BOUNDNAMES" in text and "OBS6  FILEIN  m.rch.obs" in text
    assert text.count("BEGIN period") == 4 and "MAXBOUND  3" in text
    rows = [line.split() for line in text.splitlines()
            if line.startswith("  ") and line.split()[0].isdigit()]
    assert [(r[0], r[1]) for r in rows[:3]] == [("1", "1"), ("2", "2"), ("1", "3")]
    np.testing.assert_allclose([float(r[2]) for r in rows[3:6]], rates[1], rtol=1e-6)

    with pytest.raises(lmf6.MF6Error, match=r"\(4, 3\)"):
        lmf6.recharge_package(gwf, cells, rates[:3])

    pkg = lmf6.recharge_package(gwf, cells, rates, boundnames=["z1", "z1", "z2"],
                                method="flopy", filename="m_flopy.rch")
    assert gwf.get_package("rch") is pkg
    spd = pkg.stress_period_data.get_data(2)
    np.testing.assert_allclose(spd["recharge"], rates[2])
