"""Raise NSTEP until the reference converges -- used when freezing cases.

The drainage function ks*vd^l*(1-(1-vd^(1/m))^m)^2 has an unbounded slope as
vd approaches 1 whenever m < 1, so the Picard iteration is never formally
contractive near a full store; whether it converges depends on how close the
iterates come to that neighbourhood.  Rather than encode a stability theory in
the test fixtures, each candidate case is run and its substep count raised
until the reference solves it.

This runs once, when the frozen case set is built.  The calibrated NSTEP is
then part of the committed case specification, so regeneration is
deterministic and never re-calibrates.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from .infile import Case
from .runner import OracleError, run_case

NSTEP_CAP = 256


def calibrate(exe: Path, case: Case, workdir: Path) -> tuple[Case | None, int]:
    """Return (converging case, attempts), or (None, attempts) if it never does.

    Cases that deliberately exercise non-convergence, and cases whose NSTEP is
    the thing under test, are returned untouched.
    """
    if case.allow_nonconvergence or case.name.startswith("sens_nstep"):
        return case, 0

    attempts = 0
    nstep = case.nstep
    while nstep <= NSTEP_CAP:
        attempts += 1
        candidate = replace(case, nstep=nstep)
        try:
            run_case(exe, candidate, workdir)
            return candidate, attempts
        except OracleError as exc:
            if "iteration limit" not in str(exc):
                return None, attempts       # a real error, not stiffness
            nstep *= 2
    return None, attempts
