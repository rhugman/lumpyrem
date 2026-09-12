"""Phase 1 gate: the Python kernel must reproduce the Fortran reference.

The plan sets the budget at <= 10 ULP on every column of every golden case.
The port currently achieves bit-identity, so the test asserts the budget and
reports the margin -- a regression from 0 to 3 ULP would still pass the gate
but is worth seeing in the log.
"""

from __future__ import annotations

import numpy as np
import pytest

from kernel_compare import golden_values, run_case_through
from lumpyrem.core import COLUMNS, simulate
from ulp import ulp_diff

#: docs/conversion-plan.md, Phase 1 gate.
ULP_BUDGET = 10


def test_kernel_reproduces_every_golden_case(frozen_cases, rundir):
    worst, worst_at = 0.0, None
    exact = 0
    for case in frozen_cases:
        want = golden_values(case.name, rundir / f"k_{case.name}.csv")
        got = run_case_through(simulate, case)
        assert got.shape == want.shape, (
            f"{case.name}: produced {got.shape}, reference is {want.shape}"
        )
        diffs = ulp_diff(got, want)
        peak = float(diffs.max())
        if peak == 0.0:
            exact += 1
        if peak > worst:
            _, col = np.unravel_index(int(diffs.argmax()), diffs.shape)
            worst, worst_at = peak, f"{case.name}:{COLUMNS[col]}"

    print(f"\n{exact}/{len(frozen_cases)} cases bit-identical; "
          f"worst {worst:.0f} ULP" + (f" at {worst_at}" if worst_at else ""))
    assert worst <= ULP_BUDGET, (
        f"worst disagreement {worst:.0f} ULP at {worst_at}, budget {ULP_BUDGET}"
    )


def test_output_days_match_the_reference(frozen_cases, rundir):
    """The call partitioning must line up, not just the numbers."""
    from lumpyrem.core import simulate as sim
    for case in frozen_cases[:40]:
        got = sim(
            maxvol=case.maxvol, irrigvolfrac=case.irrigvolfrac, rdelay=case.rdelay,
            mdelay=case.mdelay, ks=case.ks, m=case.m, l=case.l,
            mflowmax=case.mflowmax, maxvol_br=case.maxvol_br,
            extravol_br=case.extravol_br, gamma_br=case.gamma_br, ks_br=case.ks_br,
            m_br=case.m_br, l_br=case.l_br, offset=case.offset,
            factor1=case.factor1, factor2=case.factor2, power=case.power,
            datum=case.datum, bucket=case.bucket,
            elevmin=case.elevmin if case.elevmin is not None else -1.0e20,
            elevmax=case.elevmax if case.elevmax is not None else 1.0e20,
            vol=case.vol, vol_br=case.vol_br, rbuf=case.rbuf, mbuf=case.mbuf,
            nstep=case.nstep, mxiter=case.mxiter, tol=case.tol, rain=case.rain,
            epot=case.epot, cropfac=case.cropfac, gamma=case.gamma,
            irrigcode=case.irrigcode, gwirrigfrac=case.gwirrigfrac,
            epot_br=(case.epot_br if case.two_store else None),
            outdays=case.outdays,
        )
        expected = [0] + [min(d, case.numdays) for d in case.outdays]
        # The Fortran stops at the first output time that reaches the last day.
        cut = next((i for i, d in enumerate(expected[1:], 1) if d >= case.numdays),
                   len(expected) - 1)
        assert got.days.tolist() == expected[: cut + 1], case.name


def test_nonconvergence_is_reported_not_swallowed(frozen_cases, rundir):
    """The Fortran prints to stdout and carries on; the port counts instead."""
    byname = {c.name: c for c in frozen_cases}
    case = byname["sens_mxiter_lo"]
    assert case.allow_nonconvergence
    got = simulate(
        maxvol=case.maxvol, irrigvolfrac=case.irrigvolfrac, rdelay=case.rdelay,
        mdelay=case.mdelay, ks=case.ks, m=case.m, l=case.l, mflowmax=case.mflowmax,
        offset=case.offset, factor1=case.factor1, factor2=case.factor2,
        power=case.power, datum=case.datum, vol=case.vol, rbuf=case.rbuf,
        mbuf=case.mbuf, nstep=case.nstep, mxiter=case.mxiter, tol=case.tol,
        rain=case.rain, epot=case.epot, cropfac=case.cropfac, gamma=case.gamma,
        irrigcode=case.irrigcode, gwirrigfrac=case.gwirrigfrac, outdays=case.outdays,
    )
    assert not got.converged
    assert got.nonconverged_upper > 0

    # A case that does converge must report so.
    ok = byname["delay_fractional"]
    got_ok = run_case_through(simulate, ok)  # noqa: F841
    from lumpyrem.core import simulate as sim
    res = sim(
        maxvol=ok.maxvol, irrigvolfrac=ok.irrigvolfrac, rdelay=ok.rdelay,
        mdelay=ok.mdelay, ks=ok.ks, m=ok.m, l=ok.l, mflowmax=ok.mflowmax,
        offset=ok.offset, factor1=ok.factor1, factor2=ok.factor2, power=ok.power,
        datum=ok.datum, vol=ok.vol, rbuf=ok.rbuf, mbuf=ok.mbuf, nstep=ok.nstep,
        mxiter=ok.mxiter, tol=ok.tol, rain=ok.rain, epot=ok.epot,
        cropfac=ok.cropfac, gamma=ok.gamma, irrigcode=ok.irrigcode,
        gwirrigfrac=ok.gwirrigfrac, outdays=ok.outdays,
    )
    assert res.converged


@pytest.mark.parametrize("vd,expected", [(0.0, 0.0), (-0.5, 0.0)])
def test_response_functions_at_the_boundaries(vd, expected):
    from lumpyrem.core import drainage, evap
    assert drainage(vd, 1.0, 0.5, 1.0) == expected
    assert evap(vd, 1.0, 1.0, 4.5) == expected


def test_response_functions_saturate():
    from lumpyrem.core import drainage, evap
    assert drainage(1.0, 0.7, 0.5, 1.0) == 0.7
    assert drainage(1.5, 0.7, 0.5, 1.0) == 0.7
    assert evap(1.0, 0.004, 0.8, 4.5) == 0.8 * 0.004
