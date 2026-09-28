"""Numba-compiled build of the ``lumpyrem.core`` kernel.

``core`` stays the readable, line-by-line translation of the Fortran and the
fallback when Numba is not installed.  This module is the same arithmetic in a
form Numba can compile: no dataclasses, no ``for``-``else``, parameters packed
into a float64 row so that many cells can be stepped under ``prange``.  It is
held to the same gate as ``core`` -- bit-identity with the golden files -- so
every expression below keeps the grouping and order of its counterpart there.
Trap numbers refer to ``docs/conversion-plan.md``.

Phase 3, first step: the delay buffers are still the literal ``MAXDELAY``
arrays.  Ring buffers sized to the actual delay come after the benchmark has
shown this is worth building on.
"""

from __future__ import annotations

from math import exp

import numpy as np
from numba import njit, prange

from .core import COLUMNS, ELEVATION_FLOOR, MAXDELAY, Simulation

# -- Per-cell parameter row.  Everything a cell can differ in is a float64 so
#    one (ncell, NPARAM) array carries a whole grid.
MAXVOL, IRRIGVOLFRAC, RDELAY, MDELAY, KS, M, L, MFLOWMAX = range(8)
MAXVOL_BR, EXTRAVOL_BR, GAMMA_BR, KS_BR, M_BR, L_BR = range(8, 14)
OFFSET, FACTOR1, FACTOR2, POWER, DATUM, NBUCKET, ELEVMIN, ELEVMAX = range(14, 22)
VOL, VOL_BR = range(22, 24)
NPARAM = 24

# -- Totals accumulated over one rechmod call; mirrors core.Fluxes.
(F_RAINFALL, F_RECHARGE, F_DRAINAGE_BR, F_OVERFLOW_BR, F_MRECHARGE, F_RUNOFF,
 F_EVAPN, F_POTEVAPN, F_EVAPN_BR, F_POTEVAPN_BR, F_IRRIGATION,
 F_GWITHDRAWAL) = range(12)
NFLUX = 12

NCOLUMN = len(COLUMNS)

# error_model="numpy": a zero divisor gives inf/nan, as the Fortran does,
# instead of the per-division check Python semantics would compile in.  It
# changes nothing for finite inputs.
_JIT = dict(cache=True, nogil=True, error_model="numpy")


@njit(**_JIT)
def drainage(vd, ks, m, l):
    if vd <= 0.0:
        return 0.0
    if vd >= 1.0:
        return ks
    rtemp = 1.0 - vd ** (1.0 / m)
    rtemp = rtemp ** m
    rtemp = 1.0 - rtemp
    return ks * (vd ** l) * rtemp * rtemp


@njit(**_JIT)
def evap(vd, epot, cropfac, gamma):
    if vd <= 0.0:
        return 0.0
    if vd >= 1.0:
        return cropfac * epot
    rtemp = exp(-gamma * vd)
    return cropfac * epot * (1.0 - rtemp) / (1.0 - 2.0 * exp(-gamma) + rtemp)


@njit(**_JIT)
def _rechmod(p, state, drainsub, macsub, nonconv, iday1, iday2,
             rain, epot, cropfac_day, gamma_day, irrigcode, gwirrigfrac, epot_br,
             nstep, mxiter, tol, subdim, f):
    """One ``rechmod`` call over days ``[iday1, iday2)``.  See core.rechmod."""
    maxvol = p[MAXVOL]
    ks = p[KS]
    m = p[M]
    l = p[L]
    mflowmax = p[MFLOWMAX]
    maxvol_br = p[MAXVOL_BR]
    extravol_br = p[EXTRAVOL_BR]
    ks_br = p[KS_BR]
    m_br = p[M_BR]
    l_br = p[L_BR]
    gamma_br = p[GAMMA_BR]
    rdelay = p[RDELAY]
    mdelay = p[MDELAY]

    vol = state[0]
    vol_br = state[1]
    vvol = p[IRRIGVOLFRAC] * maxvol

    rainfall = 0.0
    recharge = 0.0
    drainage_br = 0.0
    overflow_br = 0.0
    mrecharge = 0.0
    runoff = 0.0
    evapn = 0.0
    potevapn = 0.0
    evapn_br = 0.0
    potevapn_br = 0.0
    irrigation = 0.0
    gwithdrawal = 0.0

    tstep = 1.0 / float(nstep)
    cropfac_br = 1.0            # trap 10
    recharge_for_br = 0.0
    iflag_br = 1

    # trap 7, first-pass sweep
    irdelay = int(rdelay)
    frdelay = rdelay - irdelay
    irdelay = irdelay + 1
    imdelay = int(mdelay)
    fmdelay = mdelay - imdelay
    imdelay = imdelay + 1
    if irdelay + 1 <= subdim:
        for i in range(irdelay + 1, subdim + 1):
            recharge = recharge + drainsub[i]
            drainsub[i] = 0.0
    if maxvol_br > 0.0:
        recharge_for_br = recharge
    if imdelay + 1 <= subdim:
        for i in range(imdelay + 1, subdim + 1):
            mrecharge = mrecharge + macsub[i]
            macsub[i] = 0.0

    for iday in range(iday1, iday2):
        train = rain[iday] * tstep
        tepot = epot[iday]
        tepot_br = epot_br[iday]
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
            converged = False
            for iter_ in range(1, mxiter + 1):
                vd = (tvol + vol) * 0.5 / maxvol          # trap 1
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
                if iter_ != 1:                            # trap 3
                    if ircode == 0:
                        if abs(tvol - otvol) <= tol * maxvol:
                            converged = True              # trap 2
                            break
                    else:
                        if abs(tvol - otvol) <= tol * maxvol:
                            if abs(dirrig - odirrig) <= tol * (dirrig + odirrig) * 0.5:
                                converged = True
                                break
                otvol = tvol
                odirrig = dirrig
            if not converged:
                nonconv[0] += 1

            if tempvol == 0.0:                            # traps 5 and 13
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
                if rtemp1 < mflowmax * tstep:
                    daymrech = daymrech + rtemp1
                else:
                    daymrech = daymrech + mflowmax * tstep
                    dayrunoff = dayrunoff + rtemp1 - mflowmax * tstep
                tempvol = maxvol
            vol = tempvol

        runoff = runoff + dayrunoff
        evapn = evapn + dayevap
        irrigation = irrigation + dayirrigation
        gwithdrawal = gwithdrawal + daygwithdrawal

        # trap 7, whole and fractional days
        recharge = recharge + drainsub[irdelay]
        mrecharge = mrecharge + macsub[imdelay]
        if iflag_br == 1:                                 # trap 14
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
        recharge = recharge + (1.0 - frdelay) * drainsub[irdelay]
        recharge_for_br = recharge_for_br + (1.0 - frdelay) * drainsub[irdelay]
        drainsub[irdelay] = frdelay * drainsub[irdelay]
        mrecharge = mrecharge + (1.0 - fmdelay) * macsub[imdelay]
        macsub[imdelay] = fmdelay * macsub[imdelay]
        rainfall = rainfall + rain[iday]
        potevapn = potevapn + epot[iday]

        if maxvol_br > 0.0:
            totvol_br = maxvol_br + extravol_br
            for _istep in range(nstep):
                tvol_br = vol_br
                otvol_br = vol_br
                rtemp1 = 0.0
                rtemp2 = 0.0
                tempvol = vol_br
                converged = False
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
                            converged = True
                            break
                    otvol_br = tvol_br
                if not converged:
                    nonconv[1] += 1

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

            drainage_br = drainage_br + daydrainage_br
            overflow_br = overflow_br + dayoverflow_br
            evapn_br = evapn_br + dayevap_br
            potevapn_br = potevapn_br + tepot_br

    state[0] = vol
    state[1] = vol_br
    f[F_RAINFALL] = rainfall
    f[F_RECHARGE] = recharge
    f[F_DRAINAGE_BR] = drainage_br
    f[F_OVERFLOW_BR] = overflow_br
    f[F_MRECHARGE] = mrecharge
    f[F_RUNOFF] = runoff
    f[F_EVAPN] = evapn
    f[F_POTEVAPN] = potevapn
    f[F_EVAPN_BR] = evapn_br
    f[F_POTEVAPN_BR] = potevapn_br
    f[F_IRRIGATION] = irrigation
    f[F_GWITHDRAWAL] = gwithdrawal


@njit(**_JIT)
def _elevation(p, vol, vol_br):
    """Volume to elevation, with its clamps (traps 8 and 15)."""
    dtemp = min(vol, p[MAXVOL]) if p[NBUCKET] == 1.0 else min(vol_br, p[MAXVOL_BR])
    if dtemp < ELEVATION_FLOOR:
        dtemp = ELEVATION_FLOOR
    elevation = p[OFFSET] + p[FACTOR1] * dtemp + p[FACTOR2] * (dtemp ** p[POWER])
    if elevation < p[ELEVMIN]:
        elevation = p[ELEVMIN]
    if elevation > p[ELEVMAX]:
        elevation = p[ELEVMAX]
    return elevation


@njit(**_JIT)
def _simulate_cell(p, rbuf, mbuf, rain, epot, cropfac, gamma, irrigcode,
                   gwirrigfrac, epot_br, outdays, nstep, mxiter, tol, subdim,
                   values, nonconv):
    """One cell, the ``lumprem2.f`` driver loop.  See core.simulate.

    ``values`` has one row per output time actually reached, plus day 0, and
    is filled in place; ``nonconv`` receives the upper and lower counts.
    """
    numdays = rain.shape[0]
    maxvol_br = p[MAXVOL_BR]

    state = np.empty(2)
    state[0] = p[VOL]
    state[1] = p[VOL_BR]
    drainsub = rbuf.copy()
    macsub = mbuf.copy()
    f = np.empty(NFLUX)
    nonconv[0] = 0
    nonconv[1] = 0

    totd = 0.0
    totm = 0.0
    for i in range(1, subdim + 1):
        totd = totd + drainsub[i]
        totm = totm + macsub[i]

    elevation = _elevation(p, state[0], state[1])
    values[0, :] = 0.0
    values[0, 0] = state[0]
    values[0, 1] = state[1]
    values[0, 2] = totd
    values[0, 3] = totm
    values[0, 24] = elevation
    values[0, 25] = p[DATUM] - elevation

    oldvol = state[0]
    oldvol_br = state[1]
    ototd = totd
    ototm = totm
    idayold = 0

    for imod in range(outdays.shape[0]):
        iday1 = 0 if imod == 0 else idayold
        iday2 = outdays[imod]
        ifin = 0
        if iday2 >= numdays:
            iday2 = numdays
            ifin = 1

        _rechmod(p, state, drainsub, macsub, nonconv, iday1, iday2,
                 rain, epot, cropfac, gamma, irrigcode, gwirrigfrac, epot_br,
                 nstep, mxiter, tol, subdim, f)

        deltav = state[0] - oldvol
        deltav_br = state[1] - oldvol_br
        totd = 0.0
        totm = 0.0
        for i in range(1, subdim + 1):
            totd = totd + drainsub[i]
            totm = totm + macsub[i]
        if maxvol_br > 0.0:
            allrecharge = f[F_MRECHARGE] + f[F_DRAINAGE_BR] + f[F_OVERFLOW_BR]
        else:
            allrecharge = f[F_MRECHARGE] + f[F_RECHARGE]
        rbal = (f[F_RAINFALL] + f[F_IRRIGATION] - deltav - deltav_br - allrecharge
                - f[F_RUNOFF] - f[F_EVAPN] - (totd - ototd) - (totm - ototm)
                - f[F_EVAPN_BR])

        elevation = _elevation(p, state[0], state[1])
        row = values[imod + 1]
        row[0] = state[0]
        row[1] = state[1]
        row[2] = totd
        row[3] = totm
        row[4] = deltav
        row[5] = deltav_br
        row[6] = totd - ototd
        row[7] = totm - ototm
        row[8] = f[F_RAINFALL]
        row[9] = f[F_IRRIGATION]
        row[10] = f[F_RECHARGE]
        row[11] = f[F_MRECHARGE]
        row[12] = f[F_DRAINAGE_BR]
        row[13] = f[F_OVERFLOW_BR]
        row[14] = allrecharge
        row[15] = f[F_GWITHDRAWAL]
        row[16] = allrecharge - f[F_GWITHDRAWAL]
        row[17] = f[F_RUNOFF]
        row[18] = f[F_POTEVAPN]
        row[19] = f[F_EVAPN]
        row[20] = f[F_POTEVAPN_BR]
        row[21] = f[F_EVAPN_BR]
        row[22] = f[F_POTEVAPN] - f[F_EVAPN]              # trap 9
        row[23] = rbal
        row[24] = elevation
        row[25] = p[DATUM] - elevation

        idayold = iday2
        oldvol = state[0]
        oldvol_br = state[1]
        ototd = totd
        ototm = totm
        if ifin == 1:
            break


@njit(parallel=True, **_JIT)
def _simulate_cells(params, rbuf, mbuf, rain, epot, cropfac, gamma, irrigcode,
                    gwirrigfrac, epot_br, outdays, nstep, mxiter, tol, subdim,
                    values, nonconv):
    """Independent cells under ``prange``, sharing forcing and solver settings."""
    for c in prange(params.shape[0]):
        _simulate_cell(params[c], rbuf[c], mbuf[c], rain, epot, cropfac, gamma,
                       irrigcode, gwirrigfrac, epot_br, outdays, nstep, mxiter,
                       tol, subdim, values[c], nonconv[c])


# ---------------------------------------------------------------------------
# Python-side marshalling
# ---------------------------------------------------------------------------

def _nrows(outdays: np.ndarray, numdays: int) -> int:
    """Day 0 plus every output time up to the first that reaches the end."""
    reached = np.flatnonzero(outdays >= numdays)
    return 1 + (int(reached[0]) + 1 if reached.size else outdays.size)


def param_row(
    *, maxvol, irrigvolfrac, rdelay, mdelay, ks, m, l, mflowmax,
    maxvol_br=0.0, extravol_br=0.0, gamma_br=0.0, ks_br=0.0, m_br=0.0, l_br=0.0,
    offset, factor1, factor2, power, datum, bucket="upper",
    elevmin=-1.0e20, elevmax=1.0e20, vol, vol_br=0.0,
) -> np.ndarray:
    """Pack one cell's parameters into the row layout the kernel reads."""
    if bucket not in ("upper", "lower"):
        raise ValueError(f"bucket must be 'upper' or 'lower', got {bucket!r}")
    if bucket == "lower" and maxvol_br <= 0.0:
        raise ValueError("bucket='lower' requires maxvol_br > 0")
    p = np.empty(NPARAM)
    p[MAXVOL], p[IRRIGVOLFRAC], p[RDELAY], p[MDELAY] = maxvol, irrigvolfrac, rdelay, mdelay
    p[KS], p[M], p[L], p[MFLOWMAX] = ks, m, l, mflowmax
    p[MAXVOL_BR], p[EXTRAVOL_BR], p[GAMMA_BR] = maxvol_br, extravol_br, gamma_br
    p[KS_BR], p[M_BR], p[L_BR] = ks_br, m_br, l_br
    p[OFFSET], p[FACTOR1], p[FACTOR2], p[POWER] = offset, factor1, factor2, power
    p[DATUM], p[ELEVMIN], p[ELEVMAX] = datum, elevmin, elevmax
    p[NBUCKET] = 1.0 if bucket == "upper" else 2.0
    p[VOL], p[VOL_BR] = vol, vol_br
    return p


def _buffer(values, subdim: int) -> np.ndarray:
    buf = np.zeros(subdim + 1)
    values = np.asarray(values, dtype=np.float64)
    buf[1:1 + values.size] = values
    return buf


def _forcing(rain, epot, cropfac, gamma, irrigcode, gwirrigfrac, epot_br):
    rain = np.ascontiguousarray(rain, dtype=np.float64)
    epot_br = (np.zeros_like(rain) if epot_br is None
               else np.ascontiguousarray(epot_br, dtype=np.float64))
    return (
        rain,
        np.ascontiguousarray(epot, dtype=np.float64),
        np.ascontiguousarray(cropfac, dtype=np.float64),
        np.ascontiguousarray(gamma, dtype=np.float64),
        np.ascontiguousarray(irrigcode, dtype=np.int64),
        np.ascontiguousarray(gwirrigfrac, dtype=np.float64),
        epot_br,
    )


def simulate(
    *, rbuf=(), mbuf=(), nstep, mxiter, tol,
    rain, epot, cropfac, gamma, irrigcode, gwirrigfrac, epot_br=None,
    outdays, subdim: int = MAXDELAY, **params,
) -> Simulation:
    """Drop-in replacement for ``core.simulate``, compiled."""
    p = param_row(**params)
    forcing = _forcing(rain, epot, cropfac, gamma, irrigcode, gwirrigfrac, epot_br)
    outdays = np.ascontiguousarray(outdays, dtype=np.int64)
    nrows = _nrows(outdays, forcing[0].size)
    values = np.empty((nrows, NCOLUMN))
    nonconv = np.zeros(2, dtype=np.int64)
    _simulate_cell(p, _buffer(rbuf, subdim), _buffer(mbuf, subdim), *forcing,
                   outdays, int(nstep), int(mxiter), float(tol), int(subdim),
                   values, nonconv)
    days = np.concatenate(([0], np.minimum(outdays[: nrows - 1], forcing[0].size)))
    return Simulation(
        days=days.astype(np.int64), values=values,
        nonconverged_upper=int(nonconv[0]), nonconverged_lower=int(nonconv[1]),
    )


def simulate_cells(
    params: np.ndarray, *, rbuf=None, mbuf=None, nstep, mxiter, tol,
    rain, epot, cropfac, gamma, irrigcode, gwirrigfrac, epot_br=None,
    outdays, subdim: int = MAXDELAY,
) -> tuple[np.ndarray, np.ndarray]:
    """Run ``params.shape[0]`` independent cells in parallel over shared forcing.

    ``params`` is an ``(ncell, NPARAM)`` array of rows from ``param_row``.
    Returns ``values`` of shape ``(ncell, nrows, 26)`` and the ``(ncell, 2)``
    non-convergence counts.  A benchmarking primitive for Phase 3, not the
    ``ModelGrid`` API.
    """
    params = np.ascontiguousarray(params, dtype=np.float64)
    assert params.ndim == 2 and params.shape[1] == NPARAM, params.shape
    ncell = params.shape[0]
    forcing = _forcing(rain, epot, cropfac, gamma, irrigcode, gwirrigfrac, epot_br)
    outdays = np.ascontiguousarray(outdays, dtype=np.int64)
    nrows = _nrows(outdays, forcing[0].size)
    rbuf = np.zeros((ncell, subdim + 1)) if rbuf is None else np.ascontiguousarray(rbuf)
    mbuf = np.zeros((ncell, subdim + 1)) if mbuf is None else np.ascontiguousarray(mbuf)
    values = np.empty((ncell, nrows, NCOLUMN))
    nonconv = np.zeros((ncell, 2), dtype=np.int64)
    _simulate_cells(params, rbuf, mbuf, *forcing, outdays, int(nstep), int(mxiter),
                    float(tol), int(subdim), values, nonconv)
    return values, nonconv
