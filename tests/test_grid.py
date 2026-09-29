"""Phase 3 gate, for the API: an N-cell ModelGrid run is N single runs, exactly.

``test_compiled.py`` shows the kernel's cells are independent.  This holds
the layer above it to the same standard: packing N models into parameter
rows, N initial buffers into padded rows, per-cell forcing into the kernel's
layout, and N result tables back out onto ``(time, cell)``.  Any of those can
go wrong by one cell without any number looking implausible, so the check is
exact equality with each cell run alone -- through ``Model.run`` and, as an
independent route, through ``core.simulate`` directly -- on both engines, and
the gate is shown to reject each such slip.

Cells are golden cases that share a forcing length and output schedule, so
every cell carries a genuinely different, reference-checked parameter set.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

from api_compare import EPOCH, forcing_for, model_for
from kernel_compare import run_case_through
from lumpyrem import (
    ConvergenceWarning,
    Forcing,
    GridForcing,
    Model,
    ModelError,
    ModelGrid,
    Solver,
    UpperStore,
    VolumeToElevation,
)
from lumpyrem import engine as lp_engine
from lumpyrem import model as lp_model
from lumpyrem.core import simulate
from lumpyrem.results import COLUMNS

ENGINES = ["python"] + (["compiled"] if lp_engine.compiled_available() else [])


@pytest.fixture(scope="module")
def shared_schedule(frozen_cases):
    """The largest set of golden cases on one forcing length and schedule."""
    groups: dict = {}
    for case in frozen_cases:
        key = (case.numdays, tuple(case.outdays), case.nstep, case.mxiter, case.tol)
        groups.setdefault(key, []).append(case)
    cases = max(groups.values(), key=len)
    assert len(cases) >= 8, "need several cases on one schedule to form a grid"
    assert any(c.two_store for c in cases) and not all(c.two_store for c in cases)
    return cases


def grid_forcing(cases) -> GridForcing:
    """Each case's own forcing, one column per cell."""
    def block(attr, fallback=0.0):
        return np.column_stack([
            np.asarray(getattr(c, attr)) if (attr != "epot_br" or c.two_store)
            else np.full(c.numdays, fallback) for c in cases])
    return GridForcing.from_arrays(
        EPOCH, rainfall=block("rain"), pot_evap=block("epot"),
        crop_factor=block("cropfac"), veg_gamma=block("gamma"),
        irrigate=block("irrigcode"), gw_irrig_frac=block("gwirrigfrac"),
        pot_evap_lower=block("epot_br"))


def own_forcing(case) -> Forcing:
    """What ``grid_forcing`` gives this case's column, as a Forcing of its own."""
    f = forcing_for(case)
    return f if case.two_store else Forcing.from_arrays(
        EPOCH, rainfall=case.rain, pot_evap=case.epot, crop_factor=case.cropfac,
        veg_gamma=case.gamma, irrigate=case.irrigcode, gw_irrig_frac=case.gwirrigfrac,
        pot_evap_lower=np.zeros(case.numdays))


def same_cells(grid_results, singles) -> bool:
    """Exact equality, cell by cell: values, days and convergence counts."""
    for i, one in enumerate(singles):
        got = grid_results.cell(grid_results.cells[i])
        if not (np.array_equal(got.values, one.values, equal_nan=True)
                and np.array_equal(got.days, one.days)
                and (got.nonconverged_upper, got.nonconverged_lower)
                == (one.nonconverged_upper, one.nonconverged_lower)):
            return False
    return True


# ---------------------------------------------------------------------------
# The gate
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("engine", ENGINES)
def test_grid_with_shared_forcing_is_n_single_runs(shared_schedule, engine):
    cases = shared_schedule
    forcing = forcing_for(cases[0])
    grid = ModelGrid([model_for(c) for c in cases])
    got = grid.run(forcing, times=cases[0].outdays, engine=engine)
    singles = [model_for(c).run(forcing, times=c.outdays, engine=engine) for c in cases]
    assert got.values.shape == (len(cases[0].outdays) + 1, len(cases), len(COLUMNS))
    assert same_cells(got, singles)


@pytest.mark.parametrize("engine", ENGINES)
def test_grid_with_per_cell_forcing_is_n_single_runs(shared_schedule, engine):
    cases = shared_schedule
    grid = ModelGrid([model_for(c) for c in cases])
    got = grid.run(grid_forcing(cases), times=cases[0].outdays, engine=engine)
    singles = [model_for(c).run(own_forcing(c), times=c.outdays, engine=engine)
               for c in cases]
    assert same_cells(got, singles)


def test_grid_cells_reproduce_the_kernel_called_directly(shared_schedule):
    """An independent route: no Model, no engine, just core.simulate per case."""
    cases = shared_schedule
    got = ModelGrid([model_for(c) for c in cases]).run(grid_forcing(cases),
                                                       times=cases[0].outdays)
    for i, case in enumerate(cases):
        assert np.array_equal(got.values[:, i, :], run_case_through(simulate, case)), (
            case.name)


def test_grid_cells_reproduce_the_reference(shared_schedule, reference, agreement):
    cases = shared_schedule
    got = ModelGrid([model_for(c) for c in cases]).run(grid_forcing(cases),
                                                       times=cases[0].outdays)
    for i, case in enumerate(cases):
        assert agreement.accepts(got.values[:, i, :], reference[case.name], ulp_budget=0), (
            f"{case.name}: {agreement.describe()}")


def test_engines_agree_on_a_grid(shared_schedule):
    if "compiled" not in ENGINES:
        pytest.skip("numba not installed")
    cases = shared_schedule
    grid = ModelGrid([model_for(c) for c in cases])
    a = grid.run(grid_forcing(cases), times=cases[0].outdays, engine="python")
    b = grid.run(grid_forcing(cases), times=cases[0].outdays, engine="compiled")
    assert np.array_equal(a.values, b.values)
    assert np.array_equal(a.nonconverged, b.nonconverged)


def test_cells_without_elevation_are_blanked_alone(shared_schedule):
    cases = shared_schedule[:4]
    models = [model_for(c) for c in cases]
    models[1] = models[1].with_(elevation=None)
    got = ModelGrid(models).run(forcing_for(cases[0]), times=cases[0].outdays)
    assert got.has_elevation.tolist() == [True, False, True, True]
    assert np.isnan(got.column("elevation")[:, 1]).all()
    assert not np.isnan(got.column("elevation")[:, [0, 2, 3]]).any()
    assert not got.cell(1).has_elevation


# ---------------------------------------------------------------------------
# Does the gate have teeth?
# ---------------------------------------------------------------------------

def _roll(a):
    return np.roll(a, 1, axis=0)


SLIPS = {
    "parameter rows shifted by one cell":
        lambda args, kwargs, run: (((_roll(args[0]),) + args[1:]), kwargs, run),
    "initial buffers shifted by one cell":
        lambda args, kwargs, run: ((args[0], _roll(args[1]), _roll(args[2])) + args[3:],
                                   kwargs, run),
    "per-cell forcing shifted by one cell":
        lambda args, kwargs, run: ((args[:3] + ({k: _roll(v) if v.shape[0] > 1 else v
                                                 for k, v in args[3].items()},) + args[4:]),
                                   kwargs, run),
    "results shifted by one cell":
        lambda args, kwargs, run: (args, kwargs, lambda r: lp_engine.CellRun(
            days=r.days, values=_roll(r.values), nonconverged=r.nonconverged)),
}


@pytest.mark.parametrize("slip", list(SLIPS))
def test_the_gate_rejects_a_cell_slip(slip, shared_schedule, monkeypatch):
    """Each deliberate off-by-one-cell below must fail the comparison."""
    cases = shared_schedule
    # Cells with initial buffers, so a buffer slip has something to move.
    assert sum(any(c.rbuf) or any(c.mbuf) for c in cases) >= 2
    forcing = grid_forcing(cases)
    singles = [model_for(c).run(own_forcing(c), times=c.outdays) for c in cases]
    grid = ModelGrid([model_for(c) for c in cases])
    assert same_cells(grid.run(forcing, times=cases[0].outdays), singles)

    real = lp_engine.run_cells

    def slipped(*args, **kwargs):
        args, kwargs, post = SLIPS[slip](args, kwargs, lambda r: r)
        return post(real(*args, **kwargs))

    monkeypatch.setattr(lp_model._engine, "run_cells", slipped)
    assert not same_cells(grid.run(forcing, times=cases[0].outdays), singles), (
        f"the grid gate did not notice {slip!r}")


# ---------------------------------------------------------------------------
# Building a grid
# ---------------------------------------------------------------------------

GOOD = dict(maxvol=40.0, ks=3.0, m=0.4, l=0.5)


@pytest.fixture
def forcing():
    rng = np.random.default_rng(3)
    n = 120
    return Forcing.from_arrays(
        EPOCH, rainfall=np.where(rng.random(n) < 0.3, rng.exponential(6.0, n), 0.0),
        pot_evap=np.full(n, 2.5), veg_gamma=5.0)


def test_from_arrays_broadcasts_and_matches_explicit_models(forcing):
    maxvol = np.array([30.0, 40.0, 50.0])
    grid = ModelGrid.from_arrays(
        upper=dict(maxvol=maxvol, ks=3.0, m=[0.3, 0.4, 0.5], l=0.5, rdelay=[0.0, 2.5, 30.0]),
        elevation=dict(offset=5.0, factor1=0.2, factor2=0.1, power=0.5, datum=60.0),
        initial=dict(vol=[10.0, 0.0, 25.0], drain_buffer=[(1.0,), (0.0,), (0.5, 0.5)]),
        solver=Solver(nstep=4), names=["a", "b", "c"])
    assert grid.ncell == 3 and grid.names == ("a", "b", "c")
    assert grid[2].upper.rdelay == 30.0 and grid[2].initial.drain_buffer == (0.5, 0.5)
    assert grid[0].upper.maxvol == 30.0 and grid[1].upper.ks == 3.0
    res = grid.run(forcing, times="MS")
    for i, name in enumerate("abc"):
        assert np.array_equal(res.cell(name).values, grid[i].run(forcing, times="MS").values)
    assert list(res["total_rech"].columns) == ["a", "b", "c"]


def test_from_arrays_names_the_bad_cell():
    with pytest.raises(ModelError, match=r"cell 2: UpperStore.m must be > 0") as exc:
        ModelGrid.from_arrays(upper=dict(GOOD, m=[0.4, 0.3, -1.0]))
    assert exc.value.__cause__.field == "m"
    with pytest.raises(ModelError, match="disagree on the number of cells"):
        ModelGrid.from_arrays(upper=dict(GOOD, m=[0.4, 0.3], ks=[1.0, 2.0, 3.0]))
    with pytest.raises(ModelError, match="pass ncell"):
        ModelGrid.from_arrays(upper=GOOD)
    assert ModelGrid.from_arrays(upper=GOOD, ncell=4).ncell == 4


def test_a_grid_shares_one_solver():
    a = Model(upper=UpperStore(**GOOD))
    with pytest.raises(ModelError, match="one solver"):
        ModelGrid([a, a.with_(solver=Solver(nstep=2))])
    with pytest.raises(ModelError, match="unique"):
        ModelGrid([a, a], names=["x", "x"])


def test_labelled_forcing_must_match_the_grid(forcing):
    grid = ModelGrid.from_arrays(upper=GOOD, ncell=2, names=["east", "west"])
    days = pd.date_range(EPOCH, periods=30)
    block = np.full((30, 2), 1.0)
    swapped = GridForcing.from_arrays(days, rainfall=block, pot_evap=2.0, veg_gamma=5.0,
                                      cells=["west", "east"])
    with pytest.raises(ModelError, match="do not match"):
        grid.run(swapped)
    right = GridForcing.from_arrays(days, rainfall=block, pot_evap=2.0, veg_gamma=5.0,
                                    cells=["east", "west"])
    grid.run(right)
    with pytest.raises(ModelError, match="the grid has 3"):
        ModelGrid.from_arrays(upper=GOOD, ncell=3).run(right)


def test_nonconvergence_is_reported_once_for_the_grid(frozen_cases):
    stiff = next(c for c in frozen_cases if c.allow_nonconvergence)
    models = [model_for(stiff, on_nonconvergence="warn")] * 3
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        res = ModelGrid(models).run(forcing_for(stiff), times=stiff.outdays)
    hits = [w for w in caught if issubclass(w.category, ConvergenceWarning)]
    assert len(hits) == 1 and "across 3 of 3 cells" in str(hits[0].message)
    assert res.nonconverged_cells == [0, 1, 2] and not res.converged


def test_an_unknown_engine_is_refused(forcing):
    with pytest.raises(ModelError, match="engine must be one of"):
        Model(upper=UpperStore(**GOOD)).run(forcing, engine="fortran")


def test_without_numba_auto_falls_back_to_core(forcing, monkeypatch):
    model = Model(upper=UpperStore(**GOOD, rdelay=3.5),
                  elevation=VolumeToElevation(offset=5.0, factor1=0.2, factor2=0.1,
                                              power=0.5, datum=60.0))
    want = model.run(forcing, engine="python")
    monkeypatch.setattr(lp_engine, "compiled_available", lambda: False)
    assert lp_engine.resolve("auto") == "python"
    assert np.array_equal(model.run(forcing).values, want.values)
    with pytest.raises(ImportError, match="needs Numba"):
        model.run(forcing, engine="compiled")


def test_the_package_imports_and_runs_without_numba():
    """In a fresh interpreter where ``import numba`` fails, as it does when the
    ``fast`` extra is not installed."""
    import os
    import subprocess
    import sys
    from pathlib import Path

    code = "\n".join([
        "import sys",
        "sys.modules['numba'] = None      # any import of numba now raises",
        "import numpy as np",
        "from lumpyrem import Forcing, ModelGrid, engine",
        "assert not engine.compiled_available()",
        "grid = ModelGrid.from_arrays(upper=dict(maxvol=[30.0, 40.0], ks=1.0, m=0.6, l=0.5))",
        "forcing = Forcing.from_arrays('2001-01-01', rainfall=np.ones(60), pot_evap=1.0,",
        "                              veg_gamma=4.0)",
        "print(grid.run(forcing, times='MS').values.shape)",
    ])
    src = str(Path(__file__).resolve().parents[1] / "src")
    env = {**os.environ, "PYTHONPATH": os.pathsep.join(
        [src, os.environ.get("PYTHONPATH", "")])}
    out = subprocess.run([sys.executable, "-c", code], env=env, check=True,
                         capture_output=True, text=True).stdout
    assert out.split() == ["(3,", "2,", "26)"]
