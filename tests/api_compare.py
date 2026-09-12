"""Express a frozen oracle Case through the Phase 2 object API.

The Phase 2 gate is that this route produces exactly what the Phase 1 function
call produces, so this helper is the one place where the two descriptions of a
run are mapped onto each other.  Anything the API cannot say about a golden
case would show up here as a translation that cannot be written.
"""

from __future__ import annotations

import pandas as pd

from lumpyrem import (
    Forcing,
    InitialState,
    LowerStore,
    Model,
    Solver,
    UpperStore,
    VolumeToElevation,
)

#: Golden cases carry integer simulation days and no calendar.  Any start date
#: works; this one is a Wednesday well away from a leap day, so month-start and
#: week-start schedules over the generated spans stay uneventful.
EPOCH = pd.Timestamp("2001-01-03")


def model_for(case, *, on_nonconvergence: str = "ignore") -> Model:
    return Model(
        upper=UpperStore(
            maxvol=case.maxvol, ks=case.ks, m=case.m, l=case.l,
            mflowmax=case.mflowmax, irrigvolfrac=case.irrigvolfrac,
            rdelay=case.rdelay, mdelay=case.mdelay,
        ),
        lower=(LowerStore(
            maxvol=case.maxvol_br, gamma=case.gamma_br, ks=case.ks_br,
            m=case.m_br, l=case.l_br, extravol=case.extravol_br,
        ) if case.two_store else None),
        elevation=VolumeToElevation(
            offset=case.offset, factor1=case.factor1, factor2=case.factor2,
            power=case.power, datum=case.datum, bucket=case.bucket,
            elevmin=case.elevmin, elevmax=case.elevmax,
        ),
        solver=Solver(nstep=case.nstep, mxiter=case.mxiter, tol=case.tol,
                      on_nonconvergence=on_nonconvergence),
        initial=InitialState(vol=case.vol, vol_lower=case.vol_br,
                             drain_buffer=tuple(case.rbuf),
                             macro_buffer=tuple(case.mbuf)),
    )


def forcing_for(case, *, start: pd.Timestamp = EPOCH) -> Forcing:
    """Dense daily forcing, so no fill rule fires during the fidelity check."""
    series = dict(
        rainfall=case.rain,
        pot_evap=case.epot,
        crop_factor=case.cropfac,
        veg_gamma=case.gamma,
        irrigate=case.irrigcode,
        gw_irrig_frac=case.gwirrigfrac,
    )
    if case.two_store:
        series["pot_evap_lower"] = case.epot_br
    return Forcing.from_arrays(start, **series)


def run_case_through_api(case, *, times=None, start: pd.Timestamp = EPOCH,
                         on_nonconvergence: str = "ignore"):
    return model_for(case, on_nonconvergence=on_nonconvergence).run(
        forcing_for(case, start=start),
        times=case.outdays if times is None else times,
    )


def dates_for(case, *, start: pd.Timestamp = EPOCH) -> pd.DatetimeIndex:
    """The case's output days as calendar dates, for the date-driven path."""
    return pd.DatetimeIndex([start + pd.Timedelta(days=d) for d in case.outdays])
