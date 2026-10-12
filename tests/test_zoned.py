"""ZonedForcing: forcing held per zone is the same run as forcing held per cell.

A zoned grid run marshals differently -- one forcing row per zone, a zone
index per cell, crop factor and gamma looked up from per-cell class tables --
but must compute exactly what a GridForcing run of the same daily forcing
does.  The check is exact equality, on both engines, against a GridForcing
built by expanding the zones and the class calendar by hand.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from lumpyrem import (
    Forcing,
    ForcingError,
    GridForcing,
    ModelError,
    ModelGrid,
    Solver,
    ZonedForcing,
)
from lumpyrem.results import COLUMNS

DAYS = pd.date_range("2000-01-01", "2002-12-31")
ZONES = ["a", "b", "c"]
CELL_ZONE = np.array(["b", "a", "c", "b", "a", "c", "c"])
NCELL = CELL_ZONE.size


def _inputs():
    rng = np.random.default_rng(1)
    wet = rng.random((len(DAYS), len(ZONES))) < 0.3
    rain = pd.DataFrame(rng.gamma(0.3, 8.0, wet.shape) * wet, index=DAYS, columns=ZONES)
    pet = pd.Series(4.0 + 2.0 * np.sin(np.arange(len(DAYS)) / 58.0), index=DAYS)
    years = pd.date_range("2000-01-01", "2003-01-01", freq="YS")
    veg = pd.DataFrame(rng.integers(0, 2, (len(years), len(ZONES))), index=years,
                       columns=ZONES)
    irr = pd.DataFrame(0.0, index=pd.date_range("2000-01-01", "2003-01-01", freq="MS"),
                       columns=ZONES)
    irr.loc[irr.index.month.isin([11, 12, 1]), "b"] = 1.0
    cf = rng.uniform(0.4, 1.0, (NCELL, 2))
    gm = rng.uniform(1.5, 3.0, (NCELL, 2))
    grid = ModelGrid.from_arrays(
        upper=dict(maxvol=rng.uniform(100, 200, NCELL), ks=rng.uniform(0.3, 1.0, NCELL),
                   m=0.4, l=0.5, mflowmax=1.5, irrigvolfrac=0.6,
                   rdelay=rng.uniform(1, 14, NCELL), mdelay=3.0),
        initial=dict(vol=50.0), solver=Solver(nstep=6, mxiter=200, tol=1e-5))
    return rain, pet, veg, irr, cf, gm, grid


def _by_hand(rain, pet, veg, irr, cf, gm):
    """The same forcing as a GridForcing, one column per cell."""
    cls = veg.reindex(DAYS).ffill().astype(int)
    cfb = np.column_stack([cf[i, cls[z].to_numpy()] for i, z in enumerate(CELL_ZONE)])
    gmb = np.column_stack([gm[i, cls[z].to_numpy()] for i, z in enumerate(CELL_ZONE)])
    irrd = irr.reindex(DAYS).ffill()
    return GridForcing.from_arrays(
        DAYS, rainfall=np.column_stack([rain[z] for z in CELL_ZONE]),
        pot_evap=pet.to_numpy(), irrigate=np.column_stack([irrd[z] for z in CELL_ZONE]),
        gw_irrig_frac=0.5, crop_factor=cfb, veg_gamma=gmb,
        fill={"crop_factor": "dense", "veg_gamma": "dense"})


@pytest.mark.parametrize("engine", ["python", "compiled"])
def test_zoned_matches_per_cell(engine):
    if engine == "compiled":
        pytest.importorskip("numba")
    rain, pet, veg, irr, cf, gm, grid = _inputs()
    zoned = ZonedForcing.from_frames(CELL_ZONE, rainfall=rain, pot_evap=pet, irrigate=irr,
                                     gw_irrig_frac=0.5, veg_class=veg, crop_factor=cf,
                                     veg_gamma=gm)
    ref = grid.run(_by_hand(rain, pet, veg, irr, cf, gm), times="MS", engine=engine)
    got = grid.run(zoned, times="MS", engine=engine)
    assert np.array_equal(ref.values, got.values, equal_nan=True)
    assert got["irrigation"].to_numpy().sum() > 0.0

    keep = ["net_recharge", "irrigation", "balance"]
    sub = grid.run(zoned, times="MS", engine=engine, columns=keep)
    assert sub.columns == tuple(keep)
    assert np.array_equal(sub.values, ref.values[:, :, [COLUMNS.index(c) for c in keep]])
    assert sub.balance_error().max() < 1e-9
    assert sub.cell(2).columns == tuple(keep)


def test_zoned_without_classes_matches_per_cell():
    rain, pet, _, _, _, _, grid = _inputs()
    gamma = pd.DataFrame({"a": 2.0, "b": 3.0, "c": 4.0}, index=DAYS[[0, -1]])
    zoned = ZonedForcing.from_frames(CELL_ZONE, rainfall=rain, pot_evap=pet, veg_gamma=gamma)
    ref = GridForcing.from_arrays(
        DAYS, rainfall=np.column_stack([rain[z] for z in CELL_ZONE]),
        pot_evap=pet.to_numpy(),
        veg_gamma=np.column_stack([np.full(len(DAYS), {"a": 2.0, "b": 3.0, "c": 4.0}[z])
                                   for z in CELL_ZONE]))
    a = grid.run(ref, times="MS", engine="python", columns=["net_recharge"])
    b = grid.run(zoned, times="MS", engine="python", columns=["net_recharge"])
    assert np.array_equal(a.values, b.values)


def test_engines_agree():
    pytest.importorskip("numba")
    rain, pet, veg, irr, cf, gm, grid = _inputs()
    zoned = ZonedForcing.from_frames(CELL_ZONE, rainfall=rain, pot_evap=pet, veg_class=veg,
                                     crop_factor=cf, veg_gamma=gm)
    a = grid.run(zoned, times="MS", engine="python", columns=["net_recharge"])
    b = grid.run(zoned, times="MS", engine="compiled", columns=["net_recharge"])
    assert np.array_equal(a.values, b.values)


def test_rejects_bad_input():
    rain, pet, veg, _, cf, gm, grid = _inputs()
    with pytest.raises(ForcingError, match="no column for zone"):
        ZonedForcing.from_frames(CELL_ZONE, rainfall=rain[["a", "b"]], pot_evap=pet,
                                 veg_gamma=2.0)
    with pytest.raises(ForcingError, match="veg_class codes"):
        ZonedForcing.from_frames(CELL_ZONE, rainfall=rain, pot_evap=pet, veg_class=veg + 1,
                                 crop_factor=cf, veg_gamma=gm)
    with pytest.raises(ForcingError, match="per-cell table"):
        ZonedForcing.from_frames(CELL_ZONE, rainfall=rain, pot_evap=pet, veg_class=veg,
                                 crop_factor=cf)
    zoned = ZonedForcing.from_frames(CELL_ZONE[:3], rainfall=rain, pot_evap=pet,
                                     veg_gamma=2.0)
    with pytest.raises(ModelError, match="3 cells"):
        grid.run(zoned)
    with pytest.raises(ModelError, match="unknown output column"):
        grid.run(ZonedForcing.from_frames(CELL_ZONE, rainfall=rain, pot_evap=pet,
                                          veg_gamma=2.0), columns=["recharge"])


@pytest.mark.parametrize("forcing_kind", ["zoned", "per_cell"])
def test_thread_count_does_not_change_results(forcing_kind, monkeypatch):
    pytest.importorskip("numba")
    rain, pet, veg, irr, cf, gm, grid = _inputs()
    forcing = (ZonedForcing.from_frames(CELL_ZONE, rainfall=rain, pot_evap=pet,
                                        veg_class=veg, crop_factor=cf, veg_gamma=gm)
               if forcing_kind == "zoned" else _by_hand(rain, pet, veg, irr, cf, gm))
    ref = grid.run(forcing, times="MS", engine="compiled")
    for kwargs in (dict(single_threaded=True), dict(threads=1), dict(threads=3)):
        got = grid.run(forcing, times="MS", engine="compiled", **kwargs)
        assert np.array_equal(ref.values, got.values, equal_nan=True), kwargs
    monkeypatch.setenv("LUMPYREM_THREADS", "1")
    got = grid.run(forcing, times="MS", engine="compiled")
    assert np.array_equal(ref.values, got.values, equal_nan=True)


@pytest.mark.parametrize("engine", ["python", "compiled"])
def test_shared_forcing_is_one_zone(engine):
    """Shared forcing and a single zone are the same run through the same kernel."""
    if engine == "compiled":
        pytest.importorskip("numba")
    rain, pet, _, _, _, _, grid = _inputs()
    shared = Forcing.from_dataframe(pd.DataFrame({"rainfall": rain["a"], "pot_evap": pet}),
                                    veg_gamma=2.5)
    one_zone = ZonedForcing.from_frames(["a"] * NCELL, rainfall=rain[["a"]], pot_evap=pet,
                                        veg_gamma=2.5)
    a = grid.run(shared, times="MS", engine=engine)
    b = grid.run(one_zone, times="MS", engine=engine)
    assert np.array_equal(a.values, b.values, equal_nan=True)


@pytest.mark.parametrize("engine", ["python", "compiled"])
@pytest.mark.parametrize("forcing_kind", ["shared", "per_cell"])
def test_columns_keeps_only_those_for_every_forcing(engine, forcing_kind):
    """columns= is not a zoned feature: every grid run stores only what it keeps."""
    if engine == "compiled":
        pytest.importorskip("numba")
    rain, pet, veg, irr, cf, gm, grid = _inputs()
    forcing = (_by_hand(rain, pet, veg, irr, cf, gm) if forcing_kind == "per_cell" else
               Forcing.from_dataframe(pd.DataFrame({"rainfall": rain["b"], "pot_evap": pet}),
                                      veg_gamma=2.0))
    full = grid.run(forcing, times="MS", engine=engine)
    keep = ["net_recharge", "elevation", "evap_upper"]
    sub = grid.run(forcing, times="MS", engine=engine, columns=keep)
    assert sub.values.shape == full.values.shape[:2] + (len(keep),)
    assert np.array_equal(sub.values, full.values[:, :, [COLUMNS.index(c) for c in keep]],
                          equal_nan=True)
    assert np.isnan(sub["elevation"].to_numpy()).all()     # no conversion: blanked
