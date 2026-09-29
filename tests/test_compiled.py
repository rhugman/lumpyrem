"""Phase 3 gate, correctness half: the Numba build of the kernel.

``lumpyrem.compiled`` is the same arithmetic as ``lumpyrem.core`` in a form
Numba can compile.  It is held to the Phase 1 gate against the Fortran
reference (tests/reference.py), to exact equality with ``core`` itself, and
-- the part of the Phase 3 gate that is about arrays -- an N-cell run must
reproduce N single runs exactly.

The speed half of the gate is measured by ``scripts/benchmark_numba.py``
rather than asserted here: timings on shared CI runners are too noisy to gate
on.
"""

from __future__ import annotations

import numpy as np
import pytest

from kernel_compare import run_case_through

compiled = pytest.importorskip("lumpyrem.compiled", reason="numba not installed")

from lumpyrem import core, engine  # noqa: E402
from lumpyrem.core import COLUMNS  # noqa: E402


def test_compiled_reproduces_every_golden_case(frozen_cases, reference, agreement):
    """The Phase 1 gate, applied unchanged to the compiled build."""
    worst, worst_at, exact = 0.0, None, 0
    for case in frozen_cases:
        want = reference[case.name]
        got = run_case_through(compiled.simulate, case)
        assert got.shape == want.shape, case.name
        peak, col = agreement.worst(got, want)
        exact += peak == 0.0
        if peak > worst:
            worst, worst_at = peak, f"{case.name}:{COLUMNS[col]}"
    print(f"\n{agreement.describe()}\n{exact}/{len(frozen_cases)} cases exact; "
          f"worst {worst:.3g} {agreement.unit()}" + (f" at {worst_at}" if worst_at else ""))
    assert worst <= agreement.budget(), (
        f"worst {worst:.3g} {agreement.unit()} at {worst_at}; {agreement.describe()}")


def test_compiled_matches_core_exactly(frozen_cases):
    """Not a budget: the two builds must be indistinguishable, days included."""
    for case in frozen_cases:
        kwargs = _kwargs(case)
        want = core.simulate(**kwargs)
        got = compiled.simulate(**kwargs)
        assert np.array_equal(got.values, want.values), case.name
        assert np.array_equal(got.days, want.days), case.name
        assert (got.nonconverged_upper, got.nonconverged_lower) == (
            want.nonconverged_upper, want.nonconverged_lower), case.name


def test_grid_reproduces_single_runs(frozen_cases):
    """N cells over shared forcing == N single runs, bit for bit.

    Cells are drawn from the golden cases that share one forcing length and
    output schedule, so every row carries a genuinely different parameter set.
    """
    groups: dict = {}
    for case in frozen_cases:
        key = (case.numdays, tuple(case.outdays), case.nstep, case.mxiter, case.tol)
        groups.setdefault(key, []).append(case)
    shared = max(groups.values(), key=len)
    assert len(shared) >= 8, "need several cases on one schedule to form a grid"

    lead = shared[0]
    forcing = dict(rain=lead.rain, epot=lead.epot, cropfac=lead.cropfac,
                   gamma=lead.gamma, irrigcode=lead.irrigcode,
                   gwirrigfrac=lead.gwirrigfrac,
                   epot_br=lead.epot_br if lead.two_store else None)
    solver = dict(nstep=lead.nstep, mxiter=lead.mxiter, tol=lead.tol,
                  outdays=lead.outdays)

    params = np.stack([engine.param_row(**_param_kwargs(c)) for c in shared])
    rbuf = engine.buffer_rows([c.rbuf for c in shared])
    mbuf = engine.buffer_rows([c.mbuf for c in shared])
    grid, nonconv = compiled.simulate_cells(params, rbuf=rbuf, mbuf=mbuf,
                                            **forcing, **solver)

    for i, case in enumerate(shared):
        one = compiled.simulate(**_param_kwargs(case), rbuf=case.rbuf, mbuf=case.mbuf,
                                **forcing, **solver)
        assert np.array_equal(grid[i], one.values), case.name
        assert tuple(nonconv[i]) == (one.nonconverged_upper, one.nonconverged_lower)


def test_grid_of_identical_cells_is_uniform(frozen_cases):
    case = next(c for c in frozen_cases if c.two_store)
    row = engine.param_row(**_param_kwargs(case))
    forcing = dict(rain=case.rain, epot=case.epot, cropfac=case.cropfac,
                   gamma=case.gamma, irrigcode=case.irrigcode,
                   gwirrigfrac=case.gwirrigfrac, epot_br=case.epot_br)
    rbuf = engine.buffer_rows([case.rbuf] * 33)
    mbuf = engine.buffer_rows([case.mbuf] * 33)
    grid, _ = compiled.simulate_cells(np.tile(row, (33, 1)), rbuf=rbuf, mbuf=mbuf,
                                      **forcing, nstep=case.nstep, mxiter=case.mxiter,
                                      tol=case.tol, outdays=case.outdays)
    want = run_case_through(core.simulate, case)
    for cell in grid:
        assert np.array_equal(cell, want)


def _param_kwargs(case) -> dict:
    return dict(
        maxvol=case.maxvol, irrigvolfrac=case.irrigvolfrac,
        rdelay=case.rdelay, mdelay=case.mdelay,
        ks=case.ks, m=case.m, l=case.l, mflowmax=case.mflowmax,
        maxvol_br=case.maxvol_br, extravol_br=case.extravol_br,
        gamma_br=case.gamma_br, ks_br=case.ks_br, m_br=case.m_br, l_br=case.l_br,
        offset=case.offset, factor1=case.factor1, factor2=case.factor2,
        power=case.power, datum=case.datum, bucket=case.bucket,
        elevmin=case.elevmin if case.elevmin is not None else -1.0e20,
        elevmax=case.elevmax if case.elevmax is not None else 1.0e20,
        vol=case.vol, vol_br=case.vol_br,
    )


def _kwargs(case) -> dict:
    return dict(
        **_param_kwargs(case), rbuf=case.rbuf, mbuf=case.mbuf,
        nstep=case.nstep, mxiter=case.mxiter, tol=case.tol,
        rain=case.rain, epot=case.epot, cropfac=case.cropfac, gamma=case.gamma,
        irrigcode=case.irrigcode, gwirrigfrac=case.gwirrigfrac,
        epot_br=(case.epot_br if case.two_store else None),
        outdays=case.outdays,
    )
