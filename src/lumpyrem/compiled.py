"""Numba-compiled build of the ``lumpyrem.core`` kernel.

``core`` stays the readable, line-by-line translation of the Fortran and the
fallback when Numba is not installed.  This module is the same arithmetic in a
form Numba can compile: no dataclasses, no ``for``-``else``, parameters packed
into a float64 row so that many cells can be stepped under ``prange``.  It is
held to the same gate as ``core`` -- bit-identity with the golden files -- so
every expression below keeps the grouping and order of its counterpart there.
Trap numbers refer to ``docs/conversion-plan.md``.

The one structural departure is the delay buffers.  ``lumprem2.f`` holds
each in a fixed ``MAXDELAY = 500`` array and shifts it every day; here each is
a ring of exactly ``int(delay) + 1`` slots -- the elements trap 7 ever reads --
and a day's shift is a move of the ring's head.  Nothing arithmetic changes:

- The shift moves values without touching them, so it has no rounding to
  preserve.  Logical element ``i`` (1-based, newest first, as in the Fortran)
  lives at ``(head + i - 1) % n``.
- Everything past ``int(delay) + 1`` is zero after the first call's sweep,
  and adding zero to a sum that starts at ``0.0`` is exact.  So the Fortran's
  sums over all ``MAXDELAY`` elements equal sums over the ring, taken in the
  same logical order -- which ``_ring_total`` keeps.
- Initial buffers may hold water beyond the delay.  The first call's sweep
  claims it (trap 7), summed in the Fortran's order, before the ring is
  filled; ``_rechmod`` releases that claim on its first call and zero after,
  which is what the sweep over an emptied tail produces.  Trap 14 then loses
  it from the lower store exactly as the reference does.
"""

from __future__ import annotations

from math import exp

import numpy as np
from numba import njit, prange

from .core import COLUMNS, ELEVATION_FLOOR, Simulation
from .engine import (
    DATUM,
    ELEVMAX,
    ELEVMIN,
    EXTRAVOL_BR,
    FACTOR1,
    FACTOR2,
    GAMMA_BR,
    IRRIGVOLFRAC,
    KS,
    KS_BR,
    L_BR,
    M_BR,
    MAXVOL,
    MAXVOL_BR,
    MDELAY,
    MFLOWMAX,
    NBUCKET,
    NPARAM,
    OFFSET,
    POWER,
    RDELAY,
    VOL,
    VOL_BR,
    L,
    M,
    buffer_rows,
    nrows,
    param_row,
)

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
def _tail(head, n):
    """Physical slot of logical element ``n``, the oldest in a ring of ``n``."""
    t = head - 1
    if t < 0:
        t = n - 1
    return t


@njit(**_JIT)
def _ring_total(ring, head):
    """Sum over logical elements 1 to n, in the Fortran's order."""
    n = ring.shape[0]
    total = 0.0
    k = head
    for _ in range(n):
        total = total + ring[k]
        k = k + 1
        if k == n:
            k = 0
    return total


@njit(**_JIT)
def _fill_ring(buf, n, ring):
    """Load an initial 1-based buffer into a ring of ``n``, head at slot 0.

    Returns what sits beyond logical element ``n``: the first call's sweep
    (trap 7), summed in the order the Fortran sums it.
    """
    for i in range(1, n + 1):
        ring[i - 1] = buf[i] if i < buf.shape[0] else 0.0
    claim = 0.0
    for i in range(n + 1, buf.shape[0]):
        claim = claim + buf[i]
    return claim


@njit(**_JIT)
def _rechmod(p, state, drainsub, macsub, heads, claim, nonconv, iday1, iday2,
             rain, epot, cropfac_day, gamma_day, irrigcode, gwirrigfrac, epot_br,
             nstep, mxiter, tol, f):
    """One ``rechmod`` call over days ``[iday1, iday2)``.  See core.rechmod.

    ``drainsub`` and ``macsub`` are rings of ``irdelay`` and ``imdelay``
    slots whose heads are ``heads[0]`` and ``heads[1]``; ``claim`` holds the
    sweep of the initial buffers, released once.
    """
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

    # trap 7, first-pass sweep: taken from the initial buffers by
    # _fill_ring, and nothing on any later call.
    irdelay = int(rdelay)
    frdelay = rdelay - irdelay
    irdelay = irdelay + 1
    imdelay = int(mdelay)
    fmdelay = mdelay - imdelay
    imdelay = imdelay + 1
    recharge = recharge + claim[0]
    claim[0] = 0.0
    if maxvol_br > 0.0:
        recharge_for_br = recharge
    mrecharge = mrecharge + claim[1]
    claim[1] = 0.0
    hd = heads[0]
    hm = heads[1]

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

        # trap 7, whole and fractional days.  The shift is a move of each
        # head onto its tail slot, which then takes the day's inflow as
        # logical element 1.
        td = _tail(hd, irdelay)
        tm = _tail(hm, imdelay)
        recharge = recharge + drainsub[td]
        mrecharge = mrecharge + macsub[tm]
        if iflag_br == 1:                                 # trap 14
            recharge_for_br = drainsub[td]
        else:
            recharge_for_br = recharge_for_br + drainsub[td]
            iflag_br = 0
        hd = td
        hm = tm
        drainsub[hd] = dayrech
        macsub[hm] = daymrech
        td = _tail(hd, irdelay)
        tm = _tail(hm, imdelay)
        recharge = recharge + (1.0 - frdelay) * drainsub[td]
        recharge_for_br = recharge_for_br + (1.0 - frdelay) * drainsub[td]
        drainsub[td] = frdelay * drainsub[td]
        mrecharge = mrecharge + (1.0 - fmdelay) * macsub[tm]
        macsub[tm] = fmdelay * macsub[tm]
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

    heads[0] = hd
    heads[1] = hm
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
                   gwirrigfrac, epot_br, outdays, nstep, mxiter, tol,
                   values, nonconv):
    """One cell, the ``lumprem2.f`` driver loop.  See core.simulate.

    ``rbuf`` and ``mbuf`` are the initial buffers, 1-based (index 0 unused)
    and of any length.  ``values`` has one row per output time actually
    reached, plus day 0, and is filled in place; ``nonconv`` receives the
    upper and lower counts.
    """
    numdays = rain.shape[0]
    maxvol_br = p[MAXVOL_BR]

    state = np.empty(2)
    state[0] = p[VOL]
    state[1] = p[VOL_BR]
    drainsub = np.empty(int(p[RDELAY]) + 1)
    macsub = np.empty(int(p[MDELAY]) + 1)
    heads = np.zeros(2, dtype=np.int64)
    claim = np.empty(2)
    claim[0] = _fill_ring(rbuf, drainsub.shape[0], drainsub)
    claim[1] = _fill_ring(mbuf, macsub.shape[0], macsub)
    f = np.empty(NFLUX)
    nonconv[0] = 0
    nonconv[1] = 0

    # Day 0 sums the initial buffers whole, sweep included, as the Fortran
    # does before its first call.
    totd = 0.0
    for i in range(1, rbuf.shape[0]):
        totd = totd + rbuf[i]
    totm = 0.0
    for i in range(1, mbuf.shape[0]):
        totm = totm + mbuf[i]

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

        _rechmod(p, state, drainsub, macsub, heads, claim, nonconv, iday1, iday2,
                 rain, epot, cropfac, gamma, irrigcode, gwirrigfrac, epot_br,
                 nstep, mxiter, tol, f)

        deltav = state[0] - oldvol
        deltav_br = state[1] - oldvol_br
        totd = _ring_total(drainsub, heads[0])
        totm = _ring_total(macsub, heads[1])
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


@njit(**_JIT)
def _row(a, c):
    """Cell ``c``'s row of a forcing array that is shared (one row) or per cell."""
    return a[c] if a.shape[0] > 1 else a[0]


@njit(parallel=True, **_JIT)
def _simulate_cells(params, rbuf, mbuf, rain, epot, cropfac, gamma, irrigcode,
                    gwirrigfrac, epot_br, outdays, nstep, mxiter, tol,
                    values, nonconv):
    """Independent cells under ``prange``, sharing solver settings and schedule.

    Each forcing array is ``(1, ntime)``, shared by every cell, or
    ``(ncell, ntime)``, one row per cell.
    """
    for c in prange(params.shape[0]):
        _simulate_cell(params[c], rbuf[c], mbuf[c], _row(rain, c), _row(epot, c),
                       _row(cropfac, c), _row(gamma, c), _row(irrigcode, c),
                       _row(gwirrigfrac, c), _row(epot_br, c), outdays, nstep,
                       mxiter, tol, values[c], nonconv[c])


# ---------------------------------------------------------------------------
# Python-side marshalling
# ---------------------------------------------------------------------------

def _forcing(rain, epot, cropfac, gamma, irrigcode, gwirrigfrac, epot_br, ndim):
    """Contiguous float64 arrays (int64 for the irrigation code) of ``ndim``."""
    def fix(a, dtype):
        a = np.asarray(a, dtype=dtype)
        return np.ascontiguousarray(a.reshape(1, -1) if ndim == 2 and a.ndim == 1 else a)

    rain = fix(rain, np.float64)
    epot_br = np.zeros_like(rain[:1] if ndim == 2 else rain) if epot_br is None \
        else fix(epot_br, np.float64)
    out = (rain, fix(epot, np.float64), fix(cropfac, np.float64),
           fix(gamma, np.float64), fix(irrigcode, np.int64),
           fix(gwirrigfrac, np.float64), epot_br)
    for a in out:
        assert a.ndim == ndim and a.shape[-1] == rain.shape[-1], (a.shape, rain.shape)
    return out


def simulate(
    *, rbuf=(), mbuf=(), nstep, mxiter, tol,
    rain, epot, cropfac, gamma, irrigcode, gwirrigfrac, epot_br=None,
    outdays, **params,
) -> Simulation:
    """Drop-in replacement for ``core.simulate``, compiled."""
    p = param_row(**params)
    forcing = _forcing(rain, epot, cropfac, gamma, irrigcode, gwirrigfrac, epot_br, 1)
    outdays = np.ascontiguousarray(outdays, dtype=np.int64)
    n = nrows(outdays, forcing[0].size)
    values = np.empty((n, NCOLUMN))
    nonconv = np.zeros(2, dtype=np.int64)
    _simulate_cell(p, buffer_rows([rbuf])[0], buffer_rows([mbuf])[0], *forcing,
                   outdays, int(nstep), int(mxiter), float(tol), values, nonconv)
    days = np.concatenate(([0], np.minimum(outdays[: n - 1], forcing[0].size)))
    return Simulation(
        days=days.astype(np.int64), values=values,
        nonconverged_upper=int(nonconv[0]), nonconverged_lower=int(nonconv[1]),
    )


def simulate_cells(
    params: np.ndarray, *, rbuf=None, mbuf=None, nstep, mxiter, tol,
    rain, epot, cropfac, gamma, irrigcode, gwirrigfrac, epot_br=None, outdays,
) -> tuple[np.ndarray, np.ndarray]:
    """Run ``params.shape[0]`` independent cells in parallel.

    ``params`` is an ``(ncell, NPARAM)`` array of rows from ``param_row``, and
    ``rbuf``/``mbuf`` are ``(ncell, width + 1)`` rows from ``buffer_rows``
    (empty buffers when omitted).  Each forcing variable is 1-D, or
    ``(1, ntime)``, when every cell shares it, and ``(ncell, ntime)`` when
    each has its own.  Returns ``values`` of shape ``(ncell, nrows, 26)`` and
    the ``(ncell, 2)`` non-convergence counts.  ``engine.run_cells`` is the
    entry point that also falls back to ``core``.
    """
    params = np.ascontiguousarray(params, dtype=np.float64)
    assert params.ndim == 2 and params.shape[1] == NPARAM, params.shape
    ncell = params.shape[0]
    forcing = _forcing(rain, epot, cropfac, gamma, irrigcode, gwirrigfrac, epot_br, 2)
    for a in forcing:
        assert a.shape[0] in (1, ncell), (a.shape, ncell)
    outdays = np.ascontiguousarray(outdays, dtype=np.int64)
    n = nrows(outdays, forcing[0].shape[1])
    empty = np.zeros((ncell, 1))
    rbuf = empty if rbuf is None else np.ascontiguousarray(rbuf, dtype=np.float64)
    mbuf = empty if mbuf is None else np.ascontiguousarray(mbuf, dtype=np.float64)
    assert rbuf.shape[0] == ncell and mbuf.shape[0] == ncell, (rbuf.shape, mbuf.shape)
    values = np.empty((ncell, n, NCOLUMN))
    nonconv = np.zeros((ncell, 2), dtype=np.int64)
    if ncell == 1:
        # One cell gains nothing from a parallel region but its start-up cost,
        # which a calibration loop of single runs would pay every time.
        _simulate_cell(params[0], rbuf[0], mbuf[0], *(a[0] for a in forcing), outdays,
                       int(nstep), int(mxiter), float(tol), values[0], nonconv[0])
    else:
        _simulate_cells(params, rbuf, mbuf, *forcing, outdays, int(nstep), int(mxiter),
                        float(tol), values, nonconv)
    return values, nonconv
