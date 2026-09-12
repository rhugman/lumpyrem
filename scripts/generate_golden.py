#!/usr/bin/env python3
"""Build the reference oracle and freeze its output as golden files.

This is the Phase 0 deliverable that every later phase is gated against.  It

  1. patches and compiles the vendored LUMPREM2 source (full-precision CSV),
  2. generates structured + randomised cases with dense daily forcing,
  3. raises NSTEP on any case the reference cannot solve, dropping the rest,
  4. writes the surviving case specifications to tests/data/cases.json, and
  5. writes each case's reference output to tests/data/golden/<name>.csv.gz.

Case specifications are frozen so that regeneration never re-calibrates and
never depends on the RNG: given the committed cases.json, the golden files are
a pure function of the Fortran source and the compiler.

Usage:
    python scripts/generate_golden.py            # regenerate cases and goldens
    python scripts/generate_golden.py --check    # verify, write nothing
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import platform
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "tests"))

from oracle.build import FFLAGS, build, gfortran  # noqa: E402
from oracle.calibrate import calibrate  # noqa: E402
from oracle.cases import generate_cases  # noqa: E402
from oracle.infile import dump_cases, load_cases  # noqa: E402
from oracle.runner import OracleError, run_case  # noqa: E402

DATA_DIR = REPO_ROOT / "tests" / "data"
GOLDEN_DIR = DATA_DIR / "golden"
CASES_JSON = DATA_DIR / "cases.json.gz"
MANIFEST = DATA_DIR / "manifest.json"

N_RANDOM = 200


def toolchain() -> dict:
    ver = subprocess.run([gfortran(), "--version"], capture_output=True, text=True)
    return {
        "gfortran": ver.stdout.splitlines()[0].strip() if ver.stdout else "unknown",
        "fflags": list(FFLAGS),
        "platform": platform.platform(),
        "machine": platform.machine(),
    }


def golden_path(name: str) -> Path:
    return GOLDEN_DIR / f"{name}.csv.gz"


def write_golden(name: str, csv_text: str) -> str:
    """Write the reference CSV gzipped; return its sha256 (of the plain text)."""
    with gzip.GzipFile(golden_path(name), "wb", mtime=0) as fh:
        fh.write(csv_text.encode())
    return hashlib.sha256(csv_text.encode()).hexdigest()


def read_golden(name: str) -> str:
    with gzip.open(golden_path(name), "rb") as fh:
        return fh.read().decode()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="regenerate from frozen cases.json and compare; write nothing")
    args = ap.parse_args()

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        exe = build(tmp / "build", version=2)

        if args.check:
            cases = load_cases(CASES_JSON)
            print(f"checking {len(cases)} frozen cases against {exe.name}")
        else:
            candidates = generate_cases(n_random=N_RANDOM)
            print(f"calibrating {len(candidates)} candidate cases ...")
            cases, dropped = [], []
            for case in candidates:
                solved, _ = calibrate(exe, case, tmp / "cal" / case.name)
                (cases if solved else dropped).append(solved or case.name)
            print(f"  kept {len(cases)}, dropped {len(dropped)} the reference cannot solve")
            if dropped:
                print("  dropped:", ", ".join(sorted(dropped))[:300])

        GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
        checksums, mismatches, failures = {}, [], []
        for case in cases:
            try:
                run_case(exe, case, tmp / "run" / case.name)
            except OracleError as exc:
                failures.append(f"{case.name}: {exc}")
                continue
            csv_text = (tmp / "run" / case.name / "case.csv").read_text()
            if args.check:
                if csv_text != read_golden(case.name):
                    mismatches.append(case.name)
                checksums[case.name] = hashlib.sha256(csv_text.encode()).hexdigest()
            else:
                checksums[case.name] = write_golden(case.name, csv_text)

        if failures:
            print(f"\n{len(failures)} case(s) failed to run:")
            for f in failures[:10]:
                print("   ", f)
            return 1

        if args.check:
            if mismatches:
                print(f"\nFAIL: {len(mismatches)} golden file(s) differ: "
                      f"{', '.join(mismatches[:10])}")
                return 1
            print(f"OK: all {len(checksums)} golden files reproduced byte-for-byte")
            return 0

        # Drop goldens for cases that no longer exist, so the tree stays clean.
        for stale in GOLDEN_DIR.glob("*.csv.gz"):
            if stale.stem[: -len(".csv")] not in checksums:
                stale.unlink()

        dump_cases(cases, CASES_JSON)
        MANIFEST.write_text(json.dumps(
            {"toolchain": toolchain(), "n_cases": len(cases), "sha256": checksums},
            indent=1, sort_keys=True) + "\n")

        total = sum(golden_path(n).stat().st_size for n in checksums)
        print(f"wrote {len(checksums)} golden files ({total/1e6:.2f} MB) "
              f"+ {CASES_JSON.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
