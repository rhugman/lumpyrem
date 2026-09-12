"""Faithful scalar port of LUMPREM2's ``rechmod`` kernel and its driver.

This module is a literal translation of ``rechmod2.f`` and the simulation loop
of ``lumprem2.f``, and is held to bit-fidelity against the Fortran.  It is
deliberately *not* idiomatic: expression grouping, iteration structure and
branch order all follow the original, because each of them is load-bearing.
Numbered references below point at the fidelity traps in
``docs/conversion-plan.md``.

Nothing here reads or writes files, and nothing here takes an object.  Floats
and sequences in, arrays out.  The designed API arrives in Phase 2.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import exp

import numpy as np

# lumprem2.f: parameter(MAXDELAY=500).  Phase 3 replaces the shifted buffers
# with ring buffers sized to the actual delay; Phase 1 keeps the literal form.
MAXDELAY = 500

# lumprem2.f floors the volume at `1.0e-10` before raising it to `power`, so a
# negative power cannot blow up.  That literal is a *default REAL* -- single
# precision -- in a double precision expression, so the value actually used is
# float32(1e-10) widened, not 1e-10.  With power = -0.5 the difference is
# ~7e-4, which is millions of ULP, so this is not a detail that can be tidied.
ELEVATION_FLOOR = 1.00000001335143196e-10

#: The 26 output columns, in the order lumprem2.f writes them.
COLUMNS: tuple[str, ...] = (
    "vol_upper", "vol_lower", "vol_drain", "vol_macro",
    "del_vol_upper", "del_vol_lower", "del_vol_drain", "del_vol_macro",
    "rainfall", "irrigation", "drain_upper", "macro_upper",
    "drain_lower", "overflow_lower", "total_rech", "gw_withdrawal",
    "net_recharge", "runoff", "pot_evap_upper", "evap_upper",
    "pot_evap_lower", "evap_lower", "gw_pot_evap", "balance",
    "elevation", "depth-to-water",
)


# ---------------------------------------------------------------------------
# Response functions
# ---------------------------------------------------------------------------

def drainage(vd: float, ks: float, m: float, l: float) -> float:
    """Drainage rate as a function of relative storage.

    Note the unbounded slope as ``vd`` approaches 1 for any ``m < 1``: the
    ``rtemp ** m`` term has infinite derivative at ``rtemp == 0``.  That is why
    the Picard iteration below is never formally contractive near a full store.
    """
    if vd <= 0.0:
        return 0.0
    if vd >= 1.0:
        return ks
    rtemp = 1.0 - vd ** (1.0 / m)
    rtemp = rtemp ** m
    rtemp = 1.0 - rtemp
    return ks * (vd ** l) * rtemp * rtemp


def evap(vd: float, epot: float, cropfac: float, gamma: float) -> float:
    """Evaporation rate as a function of relative storage.

    Undefined at ``gamma == 0``: the denominator is then 1 - 2 + 1 = 0 and the
    numerator is 0 too.  The Fortran reader clamps ``gamma_br`` into [0.1, 10]
    but leaves the upper store's gamma unguarded.
    """
    if vd <= 0.0:
        return 0.0
    if vd >= 1.0:
        return cropfac * epot
    rtemp = exp(-gamma * vd)
    return cropfac * epot * (1.0 - rtemp) / (1.0 - 2.0 * exp(-gamma) + rtemp)


# ---------------------------------------------------------------------------
# State and fluxes
# ---------------------------------------------------------------------------

@dataclass
class State:
    """Everything ``rechmod`` carries between calls.

    The delay buffers are 1-based to match the Fortran, so index 0 is unused
    and every subscript in the translation below reads the same as its source.
    """

    vol: float
    vol_br: float = 0.0
    drainsub: list[float] = field(default_factory=lambda: [0.0] * (MAXDELAY + 1))
    macsub: list[float] = field(default_factory=lambda: [0.0] * (MAXDELAY + 1))
    nonconverged_upper: int = 0
    nonconverged_lower: int = 0


@dataclass
class Fluxes:
    """Totals accumulated over one ``rechmod`` call."""

    rainfall: float = 0.0
    recharge: float = 0.0          # drainage released from the upper store
    drainage_br: float = 0.0
    overflow_br: float = 0.0
    mrecharge: float = 0.0
    runoff: float = 0.0
    evapn: float = 0.0
    potevapn: float = 0.0
    evapn_br: float = 0.0
    potevapn_br: float = 0.0
    irrigation: float = 0.0
    gwithdrawal: float = 0.0


# ---------------------------------------------------------------------------
# The kernel
# ---------------------------------------------------------------------------

def rechmod(
    state: State,
    *,
    mxiter: int,
    tol: float,
    nstep: int,
    maxvol: float,
    irrigvolfrac: float,
    ks: float,
    m: float,
    l: float,
    mflowmax: float,
    rdelay: float,
    mdelay: float,
    rain,
    epot,
    cropfac_day,
    gamma_day,
    irrigcode,
    gwirrigfrac,
    maxvol_br: float = 0.0,
    extravol_br: float = 0.0,
    ks_br: float = 0.0,
    m_br: float = 0.0,
    l_br: float = 0.0,
    gamma_br: float = 0.0,
    epot_br=None,
    subdim: int = MAXDELAY,
) -> Fluxes:
    """Advance ``state`` over ``len(rain)`` days.  Direct port of ``rechmod2.f``.

    Recharge to groundwater is ``recharge`` when ``maxvol_br`` is zero, and
    ``drainage_br + overflow_br`` otherwise (trap 12).
    """
    ndays = len(rain)
    vol = state.vol
    vol_br = state.vol_br
    drainsub = state.drainsub
    macsub = state.macsub

    vvol = irrigvolfrac * maxvol

    f = Fluxes()
    tstep = 1.0 / float(nstep)
    cropfac_br = 1.0            # trap 10: hardwired in the source
    recharge_for_br = 0.0
    iflag_br = 1

    # -- Buffer water that should already have been released is claimed now.
    #    (trap 7, first-pass sweep)
    irdelay = int(rdelay)
    frdelay = rdelay - irdelay
    irdelay = irdelay + 1
    imdelay = int(mdelay)
    fmdelay = mdelay - imdelay
    imdelay = imdelay + 1
    if irdelay + 1 <= subdim:
        for i in range(irdelay + 1, subdim + 1):
            f.recharge = f.recharge + drainsub[i]
            drainsub[i] = 0.0
    if maxvol_br > 0.0:
        recharge_for_br = f.recharge
    if imdelay + 1 <= subdim:
        for i in range(imdelay + 1, subdim + 1):
            f.mrecharge = f.mrecharge + macsub[i]
            macsub[i] = 0.0

    for iday in range(ndays):
        train = rain[iday] * tstep
        tepot = epot[iday]
        tepot_br = epot_br[iday] if epot_br is not None else 0.0
        cropfac = cropfac_day[iday]
        gamma = gamma_day[iday]
        ircode = irrigcode[iday]
        gwfrac = gwirrigfrac[iday]

        dayrech = 0.0
        daydrainage_br = 0.0
        dayoverflow_br = 0.0
        daymrech = 0.0
        dayevap = 0.0
        dayevap_br = 0.0
        dayrunoff = 0.0
        dayirrigation = 0.0
        daygwithdrawal = 0.0

        for _istep in range(nstep):
            tvol = vol
            otvol = vol
            odirrig = 0.0
            rtemp1 = 0.0
            rtemp2 = 0.0
            tempvol = vol
            dirrig = 0.0
            for iter_ in range(1, mxiter + 1):
                # trap 1: Crank-Nicolson, not fully implicit.
                vd = (tvol + vol) * 0.5 / maxvol
                if vd > 1.0:
                    vd = 1.0
                elif vd < 0.0:
                    vd = 0.0
                rtemp1 = drainage(vd, ks, m, l) * tstep
                rtemp2 = evap(vd, tepot, cropfac, gamma) * tstep
                tvol = vol - rtemp1 - rtemp2 + train
                if tvol < 0.0:
                    tvol = 0.0
                dirrig = 0.0
                if ircode != 0:
                    if tvol < vvol:
                        dirrig = vvol - tvol
                        tvol = vvol
                tempvol = tvol
                if tvol > maxvol:
                    tvol = maxvol
                # trap 3: the first iteration never tests convergence.
                if iter_ != 1:
                    if ircode == 0:
                        if abs(tvol - otvol) <= tol * maxvol:
                            break       # trap 2: exits carrying this iterate
                    else:
                        # trap 4: irrigation adds a second criterion.
                        if abs(tvol - otvol) <= tol * maxvol:
                            if abs(dirrig - odirrig) <= tol * (dirrig + odirrig) * 0.5:
                                break
                otvol = tvol
                odirrig = dirrig
            else:
                # The Fortran prints a warning here and carries on with
                # whatever iterate it holds.  Counted rather than raised, so
                # Phase 2 can surface it instead of losing it to stdout.
                state.nonconverged_upper += 1

            # trap 5: the store-empties fix-up.  trap 13: irrigation has
            # already raised tempvol above zero whenever it fired, so this is
            # skipped exactly when it is most needed, and mass is created.
            if tempvol == 0.0:
                rtemp3 = rtemp1 + rtemp2
                if vol < rtemp3:
                    if vol == 0.0:
                        rtemp1 = 0.0
                        rtemp2 = 0.0
                    else:
                        if rtemp3 > 0.0:
                            rtemp2 = rtemp2 * vol / rtemp3
                            rtemp1 = rtemp1 * vol / rtemp3

            dayrech = dayrech + rtemp1
            dayevap = dayevap + rtemp2
            dayirrigation = dayirrigation + dirrig
            daygwithdrawal = daygwithdrawal + gwfrac * dirrig
            if tempvol > maxvol:
                rtemp1 = tempvol - maxvol
                # trap 6: strict inequality; the branches differ on equality.
                if rtemp1 < mflowmax * tstep:
                    daymrech = daymrech + rtemp1
                else:
                    daymrech = daymrech + mflowmax * tstep
                    dayrunoff = dayrunoff + rtemp1 - mflowmax * tstep
                tempvol = maxvol
            vol = tempvol

        f.runoff = f.runoff + dayrunoff
        f.evapn = f.evapn + dayevap
        f.irrigation = f.irrigation + dayirrigation
        f.gwithdrawal = f.gwithdrawal + daygwithdrawal

        # -- Delay buffers shift, empty from the bottom, fill from the top.
        #    (trap 7, whole and fractional days)
        f.recharge = f.recharge + drainsub[irdelay]
        f.mrecharge = f.mrecharge + macsub[imdelay]
        # trap 14: iflag_br is initialised to 1 and only cleared inside the
        # branch that never runs, so this overwrites rather than accumulates.
        # The first-pass claim above is therefore never routed to the lower
        # store and simply vanishes.  Reproduced, not repaired.
        if iflag_br == 1:
            recharge_for_br = drainsub[irdelay]
        else:
            recharge_for_br = recharge_for_br + drainsub[irdelay]
            iflag_br = 0
        if irdelay > 1:
            for i in range(irdelay, 1, -1):
                drainsub[i] = drainsub[i - 1]
        if imdelay > 1:
            for i in range(imdelay, 1, -1):
                macsub[i] = macsub[i - 1]
        drainsub[1] = dayrech
        macsub[1] = daymrech
        f.recharge = f.recharge + (1.0 - frdelay) * drainsub[irdelay]
        recharge_for_br = recharge_for_br + (1.0 - frdelay) * drainsub[irdelay]
        drainsub[irdelay] = frdelay * drainsub[irdelay]
        f.mrecharge = f.mrecharge + (1.0 - fmdelay) * macsub[imdelay]
        macsub[imdelay] = fmdelay * macsub[imdelay]
        f.rainfall = f.rainfall + rain[iday]
        f.potevapn = f.potevapn + epot[iday]

        # -- Optionally route the day's recharge through the lower store.
        if maxvol_br > 0.0:
            totvol_br = maxvol_br + extravol_br
            for _istep in range(nstep):
                tvol_br = vol_br
                otvol_br = vol_br
                rtemp1 = 0.0
                rtemp2 = 0.0
                tempvol = vol_br
                for iter_ in range(1, mxiter + 1):
                    vd = (tvol_br + vol_br) * 0.5 / maxvol_br
                    if vd > 1.0:
                        vd = 1.0
                    elif vd < 0.0:
                        vd = 0.0
                    rtemp1 = drainage(vd, ks_br, m_br, l_br) * tstep
                    rtemp2 = evap(vd, tepot_br, cropfac_br, gamma_br) * tstep
                    tvol_br = vol_br - rtemp1 - rtemp2 + recharge_for_br * tstep
                    if tvol_br < 0.0:
                        tvol_br = 0.0
                    tempvol = tvol_br
                    if tvol_br > totvol_br:
                        tvol_br = totvol_br
                    if iter_ != 1:
                        if abs(tvol_br - otvol_br) <= tol * maxvol_br:
                            break
                    otvol_br = tvol_br
                else:
                    state.nonconverged_lower += 1

                if tempvol == 0.0:
                    rtemp3 = rtemp1 + rtemp2
                    if vol_br < rtemp3:
                        if vol_br == 0.0:
                            rtemp1 = 0.0
                            rtemp2 = 0.0
                        else:
                            if rtemp3 > 0.0:
                                rtemp2 = rtemp2 * vol_br / rtemp3
                                rtemp1 = rtemp1 * vol_br / rtemp3
                daydrainage_br = daydrainage_br + rtemp1
                dayevap_br = dayevap_br + rtemp2
                if tempvol > totvol_br:
                    rtemp1 = tempvol - totvol_br
                    dayoverflow_br = dayoverflow_br + rtemp1
                    tempvol = totvol_br
                vol_br = tempvol

            f.drainage_br = f.drainage_br + daydrainage_br
            f.overflow_br = f.overflow_br + dayoverflow_br
            f.evapn_br = f.evapn_br + dayevap_br
            f.potevapn_br = f.potevapn_br + tepot_br

    state.vol = vol
    state.vol_br = vol_br
    return f


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Simulation:
    """Results on the output days, plus day 0 for the initial state."""

    days: np.ndarray            # int64, shape (nout + 1,)
    values: np.ndarray          # float64, shape (nout + 1, 26)
    columns: tuple[str, ...] = COLUMNS
    nonconverged_upper: int = 0
    nonconverged_lower: int = 0

    def column(self, name: str) -> np.ndarray:
        return self.values[:, self.columns.index(name)]

    @property
    def converged(self) -> bool:
        return self.nonconverged_upper == 0 and self.nonconverged_lower == 0


def _elevation(vol, vol_br, maxvol, maxvol_br, nbucket,
               offset, factor1, factor2, power, elevmin, elevmax) -> float:
    """Volume to elevation, with its clamps (trap 8).

    The floor is applied *before* the conversion, so a negative power cannot
    blow up, and the result is then clipped into [elevmin, elevmax].
    """
    dtemp = min(vol, maxvol) if nbucket == 1 else min(vol_br, maxvol_br)
    if dtemp < ELEVATION_FLOOR:
        dtemp = ELEVATION_FLOOR
    elevation = offset + factor1 * dtemp + factor2 * (dtemp ** power)
    if elevation < elevmin:
        elevation = elevmin
    if elevation > elevmax:
        elevation = elevmax
    return elevation


def simulate(
    *,
    # -- upper store
    maxvol: float,
    irrigvolfrac: float,
    rdelay: float,
    mdelay: float,
    ks: float,
    m: float,
    l: float,
    mflowmax: float,
    # -- lower store; maxvol_br == 0 disables it
    maxvol_br: float = 0.0,
    extravol_br: float = 0.0,
    gamma_br: float = 0.0,
    ks_br: float = 0.0,
    m_br: float = 0.0,
    l_br: float = 0.0,
    # -- volume to elevation
    offset: float,
    factor1: float,
    factor2: float,
    power: float,
    datum: float,
    bucket: str = "upper",
    elevmin: float = -1.0e20,
    elevmax: float = 1.0e20,
    # -- initial conditions
    vol: float,
    vol_br: float = 0.0,
    rbuf=(),
    mbuf=(),
    # -- solver
    nstep: int,
    mxiter: int,
    tol: float,
    # -- dense daily forcing, each of length ndays
    rain,
    epot,
    cropfac,
    gamma,
    irrigcode,
    gwirrigfrac,
    epot_br=None,
    # -- output times, 1-based simulation days
    outdays,
    subdim: int = MAXDELAY,
) -> Simulation:
    """Run LUMPREM2 over dense daily forcing.  Port of the ``lumprem2.f`` loop.

    ``rechmod`` is called once per output time, over the days since the last
    one, which is what makes the call partitioning visible in the results: the
    flux columns are per-call totals, not per-day rates.
    """
    numdays = len(rain)
    if bucket not in ("upper", "lower"):
        raise ValueError(f"bucket must be 'upper' or 'lower', got {bucket!r}")
    if bucket == "lower" and maxvol_br <= 0.0:
        raise ValueError("bucket='lower' requires maxvol_br > 0")
    nbucket = 1 if bucket == "upper" else 2

    state = State(vol=float(vol), vol_br=float(vol_br))
    for i, value in enumerate(rbuf, start=1):
        state.drainsub[i] = float(value)
    for i, value in enumerate(mbuf, start=1):
        state.macsub[i] = float(value)

    totd = 0.0
    totm = 0.0
    for i in range(1, subdim + 1):
        totd = totd + state.drainsub[i]
        totm = totm + state.macsub[i]

    kernel = dict(
        mxiter=mxiter, tol=tol, nstep=nstep, maxvol=maxvol,
        irrigvolfrac=irrigvolfrac, ks=ks, m=m, l=l, mflowmax=mflowmax,
        rdelay=rdelay, mdelay=mdelay, maxvol_br=maxvol_br,
        extravol_br=extravol_br, ks_br=ks_br, m_br=m_br, l_br=l_br,
        gamma_br=gamma_br, subdim=subdim,
    )
    elev = dict(
        maxvol=maxvol, maxvol_br=maxvol_br, nbucket=nbucket, offset=offset,
        factor1=factor1, factor2=factor2, power=power,
        elevmin=elevmin, elevmax=elevmax,
    )

    # -- Day 0: the initial state, as lumprem2.f writes it before stepping.
    elevation = _elevation(state.vol, state.vol_br, **elev)
    days = [0]
    rows = [[state.vol, state.vol_br, totd, totm] + [0.0] * 20
            + [elevation, datum - elevation]]

    oldvol = state.vol
    oldvol_br = state.vol_br
    ototd = totd
    ototm = totm
    idayold = 0

    for imod, outday in enumerate(outdays):
        iday1 = 0 if imod == 0 else idayold
        iday2 = outday
        ifin = 0
        if iday2 >= numdays:
            iday2 = numdays
            ifin = 1

        sl = slice(iday1, iday2)
        f = rechmod(
            state,
            rain=rain[sl], epot=epot[sl], cropfac_day=cropfac[sl],
            gamma_day=gamma[sl], irrigcode=irrigcode[sl],
            gwirrigfrac=gwirrigfrac[sl],
            epot_br=(epot_br[sl] if epot_br is not None else None),
            **kernel,
        )

        deltav = state.vol - oldvol
        deltav_br = state.vol_br - oldvol_br
        totd = 0.0
        totm = 0.0
        for i in range(1, subdim + 1):
            totd = totd + state.drainsub[i]
            totm = totm + state.macsub[i]
        if maxvol_br > 0.0:
            allrecharge = f.mrecharge + f.drainage_br + f.overflow_br
        else:
            allrecharge = f.mrecharge + f.recharge
        rbal = (f.rainfall + f.irrigation - deltav - deltav_br - allrecharge
                - f.runoff - f.evapn - (totd - ototd) - (totm - ototm) - f.evapn_br)

        elevation = _elevation(state.vol, state.vol_br, **elev)
        dwt = datum - elevation

        # trap 9: the lower store's evaporation is deliberately not subtracted.
        gw_pot_evap = f.potevapn - f.evapn

        days.append(iday2)
        rows.append([
            state.vol, state.vol_br, totd, totm,
            deltav, deltav_br, totd - ototd, totm - ototm,
            f.rainfall, f.irrigation, f.recharge, f.mrecharge,
            f.drainage_br, f.overflow_br, allrecharge, f.gwithdrawal,
            allrecharge - f.gwithdrawal, f.runoff,
            f.potevapn, f.evapn, f.potevapn_br, f.evapn_br,
            gw_pot_evap, rbal, elevation, dwt,
        ])

        idayold = iday2
        oldvol = state.vol
        oldvol_br = state.vol_br
        ototd = totd
        ototm = totm
        if ifin == 1:
            break

    return Simulation(
        days=np.asarray(days, dtype=np.int64),
        values=np.asarray(rows, dtype=np.float64),
        nonconverged_upper=state.nonconverged_upper,
        nonconverged_lower=state.nonconverged_lower,
    )
