"""Run the reference binary on a Case and read its output -- test fixture."""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .infile import Case
from .results import Results, read_csv_output, read_rec


class OracleError(RuntimeError):
    """The reference binary failed, or reported a non-convergence warning."""


@dataclass(frozen=True)
class OracleRun:
    results: Results
    rec: dict[str, list[str]]
    stdout: str


def run_case(exe: Path, case: Case, workdir: Path, keep: bool = False,
             rewrite=None) -> OracleRun:
    """Write ``case`` into ``workdir``, run ``exe``, and read the results back.

    ``rewrite``, if given, is called with the working directory after the case
    has been written and before the binary runs.  It exists so that the
    gap-filling tests can replace the dense forcing files with sparse ones and
    see what the reference does with them -- the one thing dense cases are
    built never to exercise.
    """
    workdir = Path(workdir)
    if workdir.exists() and not keep:
        shutil.rmtree(workdir)
    infile = case.write(workdir)
    if rewrite is not None:
        rewrite(workdir)

    proc = subprocess.run(
        [str(Path(exe).resolve()), infile.name, "case.out", "case.csv"],
        cwd=workdir, capture_output=True, text=True, timeout=600,
    )
    if proc.returncode != 0:
        raise OracleError(
            f"case {case.name!r}: exit {proc.returncode}\n{proc.stdout}\n{proc.stderr}"
        )
    # LUMPREM2 warns and carries on when the Picard loop hits MXITER.  Unless
    # the case asked for that, it means the parameters are junk, not a
    # reference -- so reject it.
    hit_limit = "itn limit exceeded" in proc.stdout
    if hit_limit and not case.allow_nonconvergence:
        raise OracleError(f"case {case.name!r}: iteration limit exceeded")
    if case.allow_nonconvergence and not hit_limit:
        raise OracleError(
            f"case {case.name!r}: declared allow_nonconvergence but converged; "
            "it no longer tests what it was written to test"
        )
    # Errors are reported on stdout with a zero exit code in some paths.
    if "*** " in proc.stdout or "Cannot" in proc.stdout:
        raise OracleError(f"case {case.name!r}: reference reported an error\n{proc.stdout}")

    csv = workdir / "case.csv"
    if not csv.exists():
        raise OracleError(f"case {case.name!r}: no CSV output\n{proc.stdout}")

    return OracleRun(
        results=read_csv_output(csv),
        rec=read_rec(workdir / "lumprem_variables.rec"),
        stdout=proc.stdout,
    )
