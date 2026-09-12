"""Checks on the oracle fixtures themselves.

The Phase 0 risk is not that the Fortran is wrong -- it is that the test-only
``.in`` writer is.  A writer that dropped a parameter, or put it in the wrong
slot, would make every golden file wrong in the same direction while every
downstream gate still passed.  These tests attack that directly.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from oracle.cases import SENSITIVITY_PAIRS, structured_cases
from oracle.defects import balance_may_break, buffer_claim_may_be_lost
from oracle.runner import run_case

# lumprem_variables.rec echoes with 1pg14.7, so seven significant figures is
# all that can be asserted from it.
REC_RTOL = 1e-6

# rec name -> Case attribute.  Deliberately absent, because LUMPREM2 does not
# echo them: MFLOWMAX, NSTEP, MXITER, TOL.  Those are covered by
# test_unechoed_parameters_reach_the_model instead.
REC_SCALARS = {
    "MAXVOL": "maxvol", "IRRIGVOLFRAC": "irrigvolfrac",
    "RDELAY": "rdelay", "MDELAY": "mdelay",
    "KS": "ks", "M": "m", "L": "l",
    "OFFSET": "offset", "FACTOR1": "factor1",
    "FACTOR2": "factor2", "POWER": "power",
    "SURFACE": "datum", "VOL": "vol",
    "MAXVOL_BR": "maxvol_br", "EXTRAVOL_BR": "extravol_br",
    "GAMMA_BR": "gamma_br", "KS_BR": "ks_br", "M_BR": "m_br", "L_BR": "l_br",
    "VOL_BR": "vol_br",
}


def _sample(cases, n=40):
    """A spread across the frozen set, including every structured case."""
    structured = [c for c in cases if not c.name.startswith("random_")]
    randoms = [c for c in cases if c.name.startswith("random_")]
    return structured + randoms[:: max(1, len(randoms) // n)]


def test_rec_echo_matches_what_the_writer_intended(oracle_exe, frozen_cases, rundir):
    """Every value LUMPREM2 echoes back must be the one the writer meant to send."""
    for case in _sample(frozen_cases):
        rec = run_case(oracle_exe, case, rundir / f"rec_{case.name}").rec

        for key, attr in REC_SCALARS.items():
            expected = getattr(case, attr)
            if not case.two_store and key.endswith("_BR"):
                assert key not in rec, f"{case.name}: {key} echoed for a one-store case"
                continue
            assert key in rec, f"{case.name}: {key} missing from lumprem_variables.rec"
            got = float(rec[key][0])
            assert got == pytest.approx(expected, rel=REC_RTOL, abs=1e-12), (
                f"{case.name}: {key} echoed as {got}, writer intended {expected}"
            )

        assert int(rec["NUMDAYS"][0]) == case.numdays
        assert int(rec["NOUTDAYS"][0]) == len(case.outdays)
        assert [int(v) for v in rec["OUTDAY"]] == case.outdays
        assert int(rec["NRBUF"][0]) == len(case.rbuf)
        assert int(rec["NMBUF"][0]) == len(case.mbuf)
        assert np.allclose([float(v) for v in rec["RBUF"]], case.rbuf, rtol=REC_RTOL, atol=1e-12)
        assert np.allclose([float(v) for v in rec["MBUF"]], case.mbuf, rtol=REC_RTOL, atol=1e-12)
        if case.two_store:
            assert rec["BUCKET"][0] == case.bucket


@pytest.mark.parametrize("lo_name,hi_name,param", SENSITIVITY_PAIRS,
                         ids=[p[2] for p in SENSITIVITY_PAIRS])
def test_unechoed_parameters_reach_the_model(oracle_exe, frozen_cases, rundir,
                                             lo_name, hi_name, param):
    """MFLOWMAX, NSTEP, MXITER and TOL are not in the .rec echo.

    They are checked the only other way available: two cases differing in that
    one value must produce different output.  If the writer silently dropped
    one, these would come out identical.
    """
    byname = {c.name: c for c in frozen_cases}
    lo, hi = byname.get(lo_name), byname.get(hi_name)
    assert lo is not None and hi is not None, f"sensitivity pair for {param} missing"
    assert getattr(lo, param) != getattr(hi, param)

    a = run_case(oracle_exe, lo, rundir / lo_name).results.values
    b = run_case(oracle_exe, hi, rundir / hi_name).results.values
    assert a.shape == b.shape
    assert not np.array_equal(a, b), (
        f"changing {param} produced identical output; the writer may not be "
        f"passing it through to the model at all"
    )


def _relative_balance_error(case, res) -> float:
    # Scale by the water actually moved, so large-maxvol cases are not judged
    # against an absolute threshold they cannot meet.
    throughput = max(float(np.abs(res.column("rainfall")).sum()), case.maxvol, 1e-30)
    return float(np.abs(res.column("balance")).max()) / throughput


def test_water_balance_closes_where_the_reference_can(oracle_exe, frozen_cases, rundir):
    """The reference's balance column must be at the rounding floor.

    Excluded: cases where one of the two known reference defects can apply.
    Those are asserted directly in the two tests below rather than waived.
    """
    worst, worst_case = 0.0, None
    for case in _sample(frozen_cases):
        if balance_may_break(case):
            continue
        res = run_case(oracle_exe, case, rundir / f"bal_{case.name}").results
        err = _relative_balance_error(case, res)
        if err > worst:
            worst, worst_case = err, case.name
    assert worst < 1e-12, f"balance error {worst:.3e} in case {worst_case}"


def test_defect_irrigation_bypasses_the_store_empties_fixup(oracle_exe, frozen_cases, rundir):
    """Pin the first defect: mass is created, so balance goes negative.

    If this ever starts closing, the reference has changed and every golden
    file built from an irrigated case needs regenerating.
    """
    case = {c.name: c for c in frozen_cases}["irrig_bypasses_empty_fixup"]
    res = run_case(oracle_exe, case, rundir / "defect_irrig").results
    balance = res.column("balance")
    assert balance.min() < -1e-6, (
        "expected the store-empties fix-up to be bypassed and mass created; "
        f"worst balance was {balance.min():.3e}"
    )
    assert (balance <= 1e-12).all(), "this defect creates mass; it cannot lose any"


def test_defect_lower_store_loses_claimed_buffer_water(oracle_exe, frozen_cases, rundir):
    """Pin the second defect: exactly the claimed buffer water goes missing."""
    case = {c.name: c for c in frozen_cases}["two_store_buffer_claim_lost"]
    expected = buffer_claim_may_be_lost(case)
    assert expected > 0, "the pinned case no longer has buffer water beyond irdelay"

    res = run_case(oracle_exe, case, rundir / "defect_claim").results
    balance = res.column("balance")
    # The loss happens once, on the first call, and is exactly the claimed water.
    assert balance[1] == pytest.approx(expected, rel=1e-9), (
        f"expected {expected} of buffer water to vanish on the first call, "
        f"balance was {balance[1]}"
    )
    assert np.abs(balance[2:]).max() < 1e-12, "the loss must not recur after the first call"


# Each structured case exists to make one thing happen.  If it stops happening,
# the case is still green but no longer tests anything.
EXERCISES = {
    "dry_no_rain": "evap_upper",
    "dry_high_epot": "evap_upper",
    "deluge": "runoff",
    "no_macropore": "runoff",
    "all_macropore": "macro_upper",
    "delay_fractional": "vol_drain",
    "delay_long": "vol_drain",
    "buffers_prefilled": "vol_drain",
    "irrig_always": "irrigation",
    "irrig_blocks": "irrigation",
    "irrig_gw_split": "gw_withdrawal",
    "two_store_basic": "drain_lower",
    "two_store_extravol": "overflow_lower",
    "two_store_lower_elev": "elevation",
}


@pytest.mark.parametrize("name,column", sorted(EXERCISES.items()))
def test_structured_case_exercises_its_target(oracle_exe, frozen_cases, rundir, name, column):
    byname = {c.name: c for c in frozen_cases}
    assert name in byname, f"structured case {name} is not in the frozen set"
    res = run_case(oracle_exe, byname[name], rundir / f"ex_{name}").results
    total = float(np.abs(res.column(column)).sum())
    assert total > 0.0, f"case {name} produced nothing in {column}; it tests nothing"


def test_dense_forcing_makes_gap_filling_a_no_op(oracle_exe, rundir):
    """The premise the whole golden set rests on.

    LUMPREM2 fills gaps in its three forcing inputs by three different rules:
    vegetation is interpolated linearly, evaporation and irrigation are
    forward-filled as steps, and rainfall is zero-filled.  Phase 0 writes dense
    daily forcing so none of them fire.  This builds a case whose forcing is
    exactly what those rules would reconstruct from a sparse listing, writes it
    both ways, and requires identical output.
    """
    n = 120
    base = {c.name: c for c in structured_cases()}["two_store_basic"]

    # Piecewise linear vegetation, piecewise constant evaporation and
    # irrigation, isolated rainfall days: reconstructible by the three rules.
    cropfac = [0.4 + 0.6 * (i / (n - 1)) for i in range(n)]
    gamma = [3.0 + 2.0 * (i / (n - 1)) for i in range(n)]
    epot = [0.004 if i < 50 else 0.009 for i in range(n)]
    irrigcode = [0 if i < 70 else 1 for i in range(n)]
    gwirrigfrac = [0.0 if i < 70 else 0.5 for i in range(n)]
    rain_days = {10: 0.05, 11: 0.02, 60: 0.11}
    rain = [rain_days.get(i, 0.0) for i in range(n)]

    dense = replace(base, name="gapfill_dense", numdays=n, outdays=[30, 60, 90, 120],
                    cropfac=cropfac, gamma=gamma, epot=epot, rain=rain,
                    irrigcode=irrigcode, gwirrigfrac=gwirrigfrac,
                    epot_br=[0.001] * n)
    dense_res = run_case(oracle_exe, dense, rundir / "gapfill_dense").results

    # Same case, forcing written only at the points the rules interpolate from.
    sparse_dir = rundir / "gapfill_sparse"
    sparse = replace(dense, name="gapfill_sparse")
    sparse.write(sparse_dir)
    (sparse_dir / "veg.dat").write_text(
        f"1 {cropfac[0]!r} {gamma[0]!r}\n{n} {cropfac[-1]!r} {gamma[-1]!r}\n")
    (sparse_dir / "epot.dat").write_text(f"1 {epot[0]!r}\n51 {epot[-1]!r}\n")
    (sparse_dir / "irrig.dat").write_text("1 0 0.0\n71 1 0.5\n")
    (sparse_dir / "rain.dat").write_text(
        "".join(f"{d + 1} {v!r}\n" for d, v in sorted(rain_days.items())))
    sparse_res = run_case(oracle_exe, sparse, sparse_dir, keep=True).results

    np.testing.assert_array_equal(
        dense_res.values, sparse_res.values,
        err_msg="dense daily forcing is not equivalent to the gap-filled sparse "
                "form; the premise behind the golden set does not hold",
    )
