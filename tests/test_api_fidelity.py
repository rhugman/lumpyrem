"""Phase 2 gate: the object API must reproduce the Phase 1 numbers exactly.

Phase 1 proved the kernel bit-identical to the Fortran.  Everything Phase 2
adds sits *above* that kernel, so the only way it can move a number is by
handing the kernel something subtly different -- a NumPy scalar where the
kernel expected a Python float, a forcing array rebuilt through a resampler
that should have been a no-op, an output day off by one.  This test refuses
all of that: exact equality, not a budget.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from api_compare import EPOCH, dates_for, forcing_for, model_for, run_case_through_api
from kernel_compare import golden_values, run_case_through
from lumpyrem import Forcing
from lumpyrem.core import simulate
from lumpyrem.results import COLUMNS
from ulp import ulp_diff


def test_api_reproduces_every_golden_case(frozen_cases, rundir):
    """The gate.  Object API vs the frozen Fortran output, bit for bit."""
    worst, worst_at = 0.0, None
    for case in frozen_cases:
        want = golden_values(case.name, rundir / f"a_{case.name}.csv")
        got = run_case_through_api(case).values
        assert got.shape == want.shape, (
            f"{case.name}: API produced {got.shape}, reference is {want.shape}"
        )
        diffs = ulp_diff(got, want)
        peak = float(diffs.max())
        if peak > worst:
            _, col = np.unravel_index(int(diffs.argmax()), diffs.shape)
            worst, worst_at = peak, f"{case.name}:{COLUMNS[col]}"

    assert worst == 0.0, (
        f"the object API moved a number: worst {worst:.0f} ULP at {worst_at}. "
        f"Phase 2 adds no arithmetic, so any difference is a marshalling bug."
    )


def test_api_matches_the_phase_1_function_call(frozen_cases):
    """Independently of the goldens: the two ways of driving the kernel agree."""
    for case in frozen_cases:
        direct = run_case_through(simulate, case)
        through_api = run_case_through_api(case).values
        assert ulp_diff(direct, through_api).max() == 0.0, case.name


def test_output_times_land_on_the_right_days(frozen_cases):
    """Day numbers and calendar stamps must describe the same instants."""
    for case in frozen_cases:
        res = run_case_through_api(case)
        assert res.days.tolist() == [0] + list(case.outdays), case.name
        expected = [EPOCH + pd.Timedelta(days=int(d)) for d in res.days]
        assert list(res.times) == expected, case.name
        # Row 0 is the initial state, at the instant the run begins.
        assert res.times[0] == EPOCH
        assert res.values[0, COLUMNS.index("vol_upper")] == case.vol


def test_dates_and_day_numbers_give_identical_results(frozen_cases):
    """Asking for output by date must not differ from asking by day number."""
    for case in frozen_cases[:60]:
        by_number = run_case_through_api(case).values
        by_date = run_case_through_api(case, times=dates_for(case)).values
        assert ulp_diff(by_number, by_date).max() == 0.0, case.name


def test_daily_output_is_a_superset_of_any_schedule(frozen_cases):
    """Output times partition the run, so a sparser schedule must agree on volumes.

    Fluxes are interval totals and so change with the schedule, but the store
    volumes at a given instant cannot depend on how often the run was asked to
    report.
    """
    case = next(c for c in frozen_cases if not c.two_store and not c.allow_nonconvergence)
    daily = run_case_through_api(case, times=None)
    sparse = run_case_through_api(case)
    sparse_days = set(sparse.days.tolist())
    keep = [i for i, d in enumerate(daily.days.tolist()) if d in sparse_days]
    for name in ("vol_upper", "vol_lower", "vol_drain", "vol_macro", "elevation"):
        assert ulp_diff(daily.column(name)[keep], sparse.column(name)).max() == 0.0, name


def test_reported_fluxes_sum_to_the_same_total(frozen_cases):
    """A flux column is an interval total, so the schedule must not change its sum."""
    case = next(c for c in frozen_cases
                if c.two_store and not c.allow_nonconvergence)
    daily = run_case_through_api(case, times=None)
    sparse = run_case_through_api(case)
    for name in ("rainfall", "total_rech", "runoff", "evap_upper", "irrigation"):
        a = float(daily.column(name).sum())
        b = float(sparse.column(name).sum())
        assert a == pytest.approx(b, rel=1e-12, abs=1e-12), name


def test_a_model_without_an_elevation_conversion_says_so(frozen_cases):
    case = next(c for c in frozen_cases if not c.two_store)
    model = model_for(case).with_(elevation=None)
    res = model.run(forcing_for(case), times=case.outdays)
    assert not res.has_elevation
    assert np.isnan(res.column("elevation")).all()
    assert np.isnan(res.column("depth_to_water")).all()
    # Elevation is output-only: blanking it must not touch the physics.
    full = run_case_through_api(case)
    for name in COLUMNS[:-2]:
        assert ulp_diff(res.column(name), full.column(name)).max() == 0.0, name


# ---------------------------------------------------------------------------
# Does the gate have teeth?
# ---------------------------------------------------------------------------

def _perturbations(case):
    """One-thing-wrong versions of a correct translation.

    Phase 1 learned not to trust a green gate without checking it can go red.
    The failure modes available to a marshalling layer are mis-mapped
    parameters, a shifted output schedule and a forcing series plumbed to the
    wrong variable -- so each is committed deliberately here, and the golden
    comparison has to reject all of them.
    """
    model = model_for(case)
    forcing = forcing_for(case)
    upper = model.upper

    yield ("m and l swapped",
           model.with_(upper=replace(upper, m=upper.l, l=upper.m)),
           forcing, case.outdays)
    yield ("rdelay and mdelay swapped",
           model.with_(upper=replace(upper, rdelay=upper.mdelay, mdelay=upper.rdelay)),
           forcing, case.outdays)
    yield ("output days shifted by one",
           model, forcing, [d - 1 for d in case.outdays])
    yield ("evaporation plumbed into the crop factor",
           model,
           Forcing.from_arrays(EPOCH, rainfall=case.rain, pot_evap=case.epot,
                               crop_factor=case.epot, veg_gamma=case.gamma,
                               irrigate=case.irrigcode,
                               gw_irrig_frac=case.gwirrigfrac),
           case.outdays)


def test_the_gate_rejects_a_mis_marshalled_run(frozen_cases, rundir):
    """Every deliberate slip below must be caught by the same comparison."""
    case = next(c for c in frozen_cases
                if not c.two_store and not c.allow_nonconvergence
                and c.rdelay != c.mdelay and c.m != c.l
                and any(v != 1.0 for v in c.cropfac) and c.outdays[0] > 1)
    want = golden_values(case.name, rundir / f"teeth_{case.name}.csv")
    for label, model, forcing, times in _perturbations(case):
        got = model.run(forcing, times=times).values
        assert got.shape != want.shape or ulp_diff(got, want).max() > 0.0, (
            f"{case.name}: the gate did not notice {label!r}"
        )
    # And the same run, correctly marshalled, still passes -- so the rejections
    # above are the perturbations and not the fixture.
    assert ulp_diff(run_case_through_api(case).values, want).max() == 0.0
