"""xarray and netCDF: the layout, a lossless round trip, and forcing read back in."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from api_compare import EPOCH
from lumpyrem import Forcing, GridForcing, Model, ModelError, ModelGrid, UpperStore
from lumpyrem.results import COLUMNS

xr = pytest.importorskip("xarray")
from lumpyrem import io  # noqa: E402

GOOD = dict(ks=1.0, m=0.6, l=0.5, mflowmax=1.0)


@pytest.fixture
def grid_run():
    rng = np.random.default_rng(5)
    ndays, ncell = 200, 4
    rain = np.where(rng.random((ndays, ncell)) < 0.3, rng.exponential(6.0, (ndays, ncell)),
                    0.0)
    forcing = GridForcing.from_arrays(EPOCH, rainfall=rain, pot_evap=np.full(ndays, 2.5),
                                      veg_gamma=5.0, cells=["n", "e", "s", "w"])
    grid = ModelGrid.from_arrays(
        upper=dict(GOOD, maxvol=[30.0, 40.0, 50.0, 60.0], rdelay=[0.0, 1.5, 10.0, 45.25]),
        names=forcing.cells)
    return grid, forcing, grid.run(forcing, times="MS")


def test_a_grid_becomes_one_variable_per_column_on_time_and_cell(grid_run):
    _, _, res = grid_run
    ds = res.to_xarray()
    assert set(COLUMNS) <= set(ds.data_vars)
    assert ds["total_rech"].dims == ("time", "cell")
    assert list(ds["cell"].values) == ["n", "e", "s", "w"]
    assert (ds["time"].values == res.times.values).all()
    assert (ds["day"].values == res.days).all()
    for name in COLUMNS:
        assert np.array_equal(ds[name].values, res.column(name), equal_nan=True), name
    assert ds["total_rech"].sel(cell="s").values.tolist() == \
        res.cell("s").column("total_rech").tolist()
    assert ds["vol_upper"].attrs["cell_methods"] == "time: point"
    assert ds["nonconverged_upper"].dims == ("cell",)


def test_netcdf_round_trip_is_exact(grid_run, tmp_path):
    _, _, res = grid_run
    path = tmp_path / "grid.nc"
    res.to_netcdf(path)
    with xr.open_dataset(path) as back:
        back = back.load()
    for name in COLUMNS:
        assert np.array_equal(back[name].values, res.column(name), equal_nan=True), name
    assert back["has_elevation"].values.tolist() == res.has_elevation.tolist()


def test_forcing_read_from_netcdf_drives_the_same_run(grid_run, tmp_path):
    grid, forcing, res = grid_run
    days = forcing.index
    ds = xr.Dataset(
        {"rainfall": (("time", "cell"), forcing.blocks["rainfall"]),
         "pot_evap": ("time", np.full(len(days), 2.5)),
         "veg_gamma": ((), 5.0)},
        coords={"time": days, "cell": list(forcing.cells)})
    path = tmp_path / "forcing.nc"
    ds.to_netcdf(path)
    back = io.read_forcing(path)
    assert isinstance(back, GridForcing) and back.cells == ("n", "e", "s", "w")
    again = grid.run(back, times="MS")
    assert np.array_equal(again.values, res.values, equal_nan=True)


def test_labels_read_from_a_file_protect_the_cell_order(grid_run):
    """The teeth of the label check: the same data with the cells reversed."""
    grid, forcing, _ = grid_run
    ds = xr.Dataset(
        {"rainfall": (("time", "cell"), forcing.blocks["rainfall"][:, ::-1]),
         "pot_evap": ("time", np.full(len(forcing.index), 2.5)), "veg_gamma": ((), 5.0)},
        coords={"time": forcing.index, "cell": list(forcing.cells)[::-1]})
    with pytest.raises(ModelError, match="do not match"):
        grid.run(io.forcing_from_dataset(ds))


def test_forcing_without_a_cell_axis_is_a_plain_forcing():
    days = pd.date_range(EPOCH, periods=30)
    ds = xr.Dataset({"rainfall": ("time", np.ones(30)), "pot_evap": ("time", np.ones(30)),
                     "veg_gamma": ((), 4.0)}, coords={"time": days})
    f = io.forcing_from_dataset(ds)
    assert isinstance(f, Forcing)
    Model(upper=UpperStore(maxvol=40.0, **GOOD)).run(f)


def test_forcing_with_unknown_variables_is_refused():
    days = pd.date_range(EPOCH, periods=5)
    ds = xr.Dataset({"rain": ("time", np.ones(5))}, coords={"time": days})
    with pytest.raises(Exception, match="unknown forcing variable"):
        io.forcing_from_dataset(ds)


def test_a_single_run_becomes_a_dataset_on_time():
    days = pd.date_range(EPOCH, periods=30)
    res = Model(upper=UpperStore(maxvol=40.0, **GOOD)).run(
        Forcing.from_arrays(days[0], rainfall=np.ones(30), pot_evap=1.0, veg_gamma=4.0))
    ds = res.to_xarray()
    assert ds["total_rech"].dims == ("time",)
    assert np.array_equal(ds["total_rech"].values, res.column("total_rech"))


def test_a_per_cell_constant_is_broadcast_over_time(grid_run):
    grid, forcing, res = grid_run
    ds = xr.Dataset(
        {"rainfall": (("time", "cell"), forcing.blocks["rainfall"]),
         "pot_evap": ("time", np.full(len(forcing.index), 2.5)),
         "veg_gamma": ((), 5.0), "crop_factor": ("cell", [1.0, 0.9, 0.8, 0.7])},
        coords={"time": forcing.index, "cell": list(forcing.cells)})
    daily = io.forcing_from_dataset(ds).to_daily()
    assert daily["crop_factor"].shape == (len(forcing.index), 4)
    assert (daily["crop_factor"] == [1.0, 0.9, 0.8, 0.7]).all()


@pytest.fixture
def stacked_run():
    """Forcing on (time, y, x), stacked so that ``cell`` is a (y, x) MultiIndex."""
    days = pd.date_range(EPOCH, periods=200, name="time")
    rng = np.random.default_rng(7)
    rain = np.where(rng.random((200, 2, 3)) < 0.3, rng.exponential(6.0, (200, 2, 3)), 0.0)
    ds = xr.Dataset({"rainfall": (("time", "y", "x"), rain),
                     "pot_evap": ("time", np.full(200, 2.5)), "veg_gamma": 5.0},
                    coords={"time": days, "y": [125.0, 375.0], "x": [0.0, 250.0, 500.0]})
    forcing = io.forcing_from_dataset(ds.stack(cell=("y", "x")))
    grid = ModelGrid.from_arrays(upper=dict(GOOD, maxvol=40.0), names=forcing.cells)
    return ds, forcing, grid.run(forcing, times="MS")


def test_tuple_cell_labels_become_a_multiindex_that_unstacks(stacked_run):
    ds, forcing, res = stacked_run
    assert forcing.cells[1] == (125.0, 250.0)
    out = res.to_xarray(cell_levels=("y", "x"))
    assert isinstance(out.indexes["cell"], pd.MultiIndex)
    maps = out.unstack("cell")
    assert maps["total_rech"].dims == ("time", "y", "x")
    assert maps["y"].values.tolist() == ds["y"].values.tolist()
    assert maps["total_rech"].sel(y=375.0, x=0.0).values.tolist() == \
        res.cell((375.0, 0.0)).column("total_rech").tolist()
    assert list(res.to_xarray().indexes["cell"].names) == ["cell_level_0", "cell_level_1"]


def test_a_multiindex_is_written_to_netcdf_as_its_levels(stacked_run, tmp_path):
    _, _, res = stacked_run
    path = tmp_path / "grid.nc"
    res.to_netcdf(path, cell_levels=("y", "x"))
    with xr.open_dataset(path) as back:
        back = back.load()
    assert "cell" not in back.coords and back["y"].dims == ("cell",)
    again = back.set_index(cell=["y", "x"])
    assert again.indexes["cell"].tolist() == list(res.cells)
    assert np.array_equal(again["total_rech"].values, res.column("total_rech"))


def test_cell_levels_must_fit_the_labels(stacked_run, grid_run):
    _, _, res = stacked_run
    with pytest.raises(ValueError, match="does not match"):
        res.to_xarray(cell_levels=("layer", "row", "col"))
    _, _, plain = grid_run
    with pytest.raises(ValueError, match="not tuples"):
        plain.to_xarray(cell_levels=("y", "x"))
