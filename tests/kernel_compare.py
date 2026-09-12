"""Shared helper: run the Python kernel over a frozen case and compare."""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np

from oracle.results import read_csv_output

GOLDEN_DIR = Path(__file__).resolve().parent / "data" / "golden"


def golden_values(name: str, scratch: Path) -> np.ndarray:
    """The reference output for ``name``, as a (nout+1, 26) float64 array."""
    scratch.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(GOLDEN_DIR / f"{name}.csv.gz", "rb") as fh:
        scratch.write_bytes(fh.read())
    return read_csv_output(scratch).values


def run_case_through(simulate, case) -> np.ndarray:
    """Drive ``simulate`` from a frozen oracle Case.

    ``simulate`` is passed in rather than imported so that mutation tests can
    supply a deliberately broken build of the kernel.
    """
    return simulate(
        maxvol=case.maxvol, irrigvolfrac=case.irrigvolfrac,
        rdelay=case.rdelay, mdelay=case.mdelay,
        ks=case.ks, m=case.m, l=case.l, mflowmax=case.mflowmax,
        maxvol_br=case.maxvol_br, extravol_br=case.extravol_br,
        gamma_br=case.gamma_br, ks_br=case.ks_br, m_br=case.m_br, l_br=case.l_br,
        offset=case.offset, factor1=case.factor1, factor2=case.factor2,
        power=case.power, datum=case.datum, bucket=case.bucket,
        elevmin=case.elevmin if case.elevmin is not None else -1.0e20,
        elevmax=case.elevmax if case.elevmax is not None else 1.0e20,
        vol=case.vol, vol_br=case.vol_br, rbuf=case.rbuf, mbuf=case.mbuf,
        nstep=case.nstep, mxiter=case.mxiter, tol=case.tol,
        rain=case.rain, epot=case.epot, cropfac=case.cropfac, gamma=case.gamma,
        irrigcode=case.irrigcode, gwirrigfrac=case.gwirrigfrac,
        epot_br=(case.epot_br if case.two_store else None),
        outdays=case.outdays,
    ).values
