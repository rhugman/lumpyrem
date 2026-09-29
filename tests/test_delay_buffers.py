"""The delay buffers, sized to the delay, beyond where the reference can reach.

``lumprem2.f`` holds each delay buffer in a fixed 500-element array, so the
golden set stops at delays of 498 days.  The buffers are now sized per run:
``core`` keeps the Fortran's shifting array at the length the run needs, and
``compiled`` holds a ring of exactly ``int(delay) + 1`` slots.  Within the
reference's range both are held to it (test_compiled.py,
test_compiled_mutations.py).  Beyond it they are two independent
implementations of the same buffer and are held to each other, exactly, and
to the physics: a delayed pulse arrives when it should, and nothing is lost.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from lumpyrem import Forcing, InitialState, LowerStore, Model, Solver, UpperStore
from lumpyrem import engine as lp_engine
from lumpyrem.core import buffer_length

EPOCH = pd.Timestamp("2001-01-03")
DAY = pd.Timedelta(days=1)

needs_numba = pytest.mark.skipif(not lp_engine.compiled_available(),
                                 reason="numba not installed")


def random_case(rng, ndays):
    """A model whose delays and initial buffers may exceed the Fortran's 500."""
    long = rng.random() < 0.5
    rdelay = float(rng.uniform(0, 900) if long else rng.uniform(0, 40))
    mdelay = float(rng.choice([0.0, rng.uniform(0, 30), rng.uniform(400, 900)]))
    if rng.random() < 0.3:
        rdelay = float(int(rdelay))          # whole days: frdelay == 0
    nbuf = int(rng.integers(1, 1200))
    lower = LowerStore(maxvol=60.0, gamma=2.0, ks=0.4, m=0.5, l=0.5, extravol=5.0) \
        if rng.random() < 0.5 else None
    model = Model(
        upper=UpperStore(maxvol=50.0, ks=float(rng.uniform(0.5, 6.0)), m=0.45, l=0.5,
                         mflowmax=1.0, rdelay=rdelay, mdelay=mdelay),
        lower=lower,
        solver=Solver(nstep=3, tol=1e-6, on_nonconvergence="ignore"),
        initial=InitialState(
            vol=20.0, vol_lower=10.0 if lower else 0.0,
            drain_buffer=tuple(rng.exponential(0.05, nbuf)),
            macro_buffer=tuple(rng.exponential(0.05, int(rng.integers(1, 800))))),
    )
    forcing = Forcing.from_arrays(
        EPOCH, rainfall=np.where(rng.random(ndays) < 0.3, rng.exponential(8.0, ndays), 0.0),
        pot_evap=np.full(ndays, 2.0), veg_gamma=4.0, pot_evap_lower=0.5)
    return model, forcing


@needs_numba
@pytest.mark.parametrize("seed", range(12))
def test_the_ring_and_the_shifting_array_agree_beyond_the_reference(seed):
    rng = np.random.default_rng(seed)
    ndays = 1500
    model, forcing = random_case(rng, ndays)
    for times in (None, "MS", sorted(set(rng.integers(1, ndays + 1, 7).tolist()))):
        ring = model.run(forcing, times=times, engine="compiled")
        shift = model.run(forcing, times=times, engine="python")
        assert np.array_equal(ring.values, shift.values, equal_nan=True), (seed, times)
        assert (ring.nonconverged_upper, ring.nonconverged_lower) == (
            shift.nonconverged_upper, shift.nonconverged_lower)


@pytest.mark.parametrize("engine", ["python"] + (
    ["compiled"] if lp_engine.compiled_available() else []))
def test_a_pulse_arrives_after_a_delay_longer_than_the_fortran_allowed(engine):
    """One day of rain and a 600.25-day delay.  Drainage from day 1 sits in the
    buffer until day 601, when 1 - frdelay of it is released and the rest is
    held a day longer; and every drop of the rain is accounted for."""
    ndays = 700
    rain = np.zeros(ndays)
    rain[0] = 30.0
    model = Model(upper=UpperStore(maxvol=50.0, ks=10.0, m=0.5, l=0.5, rdelay=600.25),
                  initial=InitialState(vol=0.0))
    forcing = Forcing.from_arrays(EPOCH, rainfall=rain, pot_evap=np.zeros(ndays),
                                  veg_gamma=4.0)
    res = model.run(forcing, engine=engine)
    released = res.column("total_rech")
    in_buffer = res.column("vol_drain")

    assert released[1:601].max() == 0.0
    assert released[601] == (1.0 - 0.25) * in_buffer[1]     # day 1's drainage, moved
    assert res.balance_error() < 1e-12
    held = res.column("vol_upper")[-1] + in_buffer[-1]
    assert abs(30.0 - held - released.sum() - res.column("runoff").sum()) < 1e-12


def test_buffers_are_sized_to_the_delay_not_to_maxdelay():
    assert buffer_length(0.0, 0.0) == 1
    assert buffer_length(12.5, 1.5) == 13
    assert buffer_length(2.0, 0.0, rbuf=(0.0,) * 40) == 40
    assert buffer_length(900.9, 3.0) == 901
