"""Readers for the reference binary's output -- test fixtures only.

Two files are read back:

``case.csv``
    The full-precision results table.  LUMPREM2 writes the same 26 columns to
    both its fixed-width tabular file and its optional CSV file; only the CSV
    writer is widened to 17 significant figures (see ``build.py``), which
    leaves the tabular file byte-comparable with the shipped 2020 examples.

``lumprem_variables.rec``
    The values LUMPREM2 echoes back after parsing.  Asserting this against
    what the writer *intended* is the guard against a silently buggy ``.in``
    writer -- the failure mode where every golden file is wrong the same way
    and the gates pass regardless.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .infile import COLUMNS


@dataclass(frozen=True)
class Results:
    """Reference output: output days and the 26 columns at each."""

    days: np.ndarray          # int, shape (nout,), first entry is day 0
    values: np.ndarray        # float64, shape (nout, 26)
    columns: tuple[str, ...]

    def column(self, name: str) -> np.ndarray:
        return self.values[:, self.columns.index(name)]


def read_csv_output(path: Path) -> Results:
    lines = [ln.rstrip("\n") for ln in Path(path).read_text().splitlines() if ln.strip()]
    if not lines:
        raise ValueError(f"{path} is empty")

    # The Fortran header format is `format(a,1000(',',a))`, which emits the
    # separator literal once more before it runs out of list items, so the
    # header carries a trailing comma the data rows do not have.
    header = [f for f in lines[0].split(",") if f != ""]
    if tuple(header[1:]) != COLUMNS:
        raise ValueError(f"unexpected columns in {path}: {header}")

    days, rows = [], []
    for ln in lines[1:]:
        fields = ln.split(",")
        if len(fields) != len(COLUMNS) + 1:
            raise ValueError(f"{path}: expected {len(COLUMNS) + 1} fields, got {len(fields)}")
        days.append(int(fields[0]))
        rows.append([float(f) for f in fields[1:]])

    return Results(np.asarray(days, dtype=np.int64),
                   np.asarray(rows, dtype=np.float64), COLUMNS)


def read_rec(path: Path) -> dict[str, list[str]]:
    """Parse ``lumprem_variables.rec`` into name -> list of raw string values.

    Note what is *absent*: MFLOWMAX appears nowhere in this file, and the whole
    `* solution parameters` block (NSTEP, MXITER, TOL) is omitted.  Those four
    cannot be checked this way and are covered instead by asserting that
    varying them changes the output.
    """
    out: dict[str, list[str]] = {}
    pending: str | None = None
    for raw in Path(path).read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("*") or line.startswith("Variables used"):
            continue
        if line.endswith("--->"):
            pending = line[: -len("--->")].strip()
            out[pending] = []
            continue
        parts = line.split()
        if pending is not None and _all_numeric(parts):
            out[pending].extend(parts)
            continue
        pending = None
        if len(parts) >= 2:
            out[parts[0]] = parts[1:]
    return out


def _all_numeric(parts: list[str]) -> bool:
    try:
        for p in parts:
            float(p)
    except ValueError:
        return False
    return bool(parts)
