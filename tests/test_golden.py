"""The Phase 0 gate: golden files must regenerate reproducibly.

On the platform that produced them the requirement is byte-identity.  Anywhere
else it is a small ULP budget, because ``exp`` and ``**`` come from the system
maths library and the last bit of a transcendental is not standardised across
glibc, Apple's libm and mingw.  Everything else in the kernel is +, -, * and /,
which IEEE 754 pins exactly and which -ffp-contract=off keeps un-fused.
"""

from __future__ import annotations

import gzip
import hashlib
import os
import platform
from pathlib import Path

import numpy as np

from oracle.defects import balance_may_break
from oracle.results import read_csv_output
from oracle.runner import run_case
from ulp import ulp_diff

GOLDEN_DIR = Path(__file__).resolve().parent / "data" / "golden"

# Budget used only when the running platform is not the one that generated the
# golden files.  Provisional until CI has measured the real figure.
CROSS_PLATFORM_ULP = int(os.environ.get("LUMPYREM_GOLDEN_ULP", "16"))


def _golden_text(name: str) -> str:
    with gzip.open(GOLDEN_DIR / f"{name}.csv.gz", "rb") as fh:
        return fh.read().decode()


def _same_platform(manifest: dict) -> bool:
    tc = manifest["toolchain"]
    return (tc["machine"] == platform.machine()
            and tc["platform"].split("-")[0] == platform.platform().split("-")[0])


def test_every_frozen_case_has_a_golden_file(frozen_cases):
    frozen = {c.name for c in frozen_cases}
    stored = {p.name[: -len(".csv.gz")] for p in GOLDEN_DIR.glob("*.csv.gz")}
    assert frozen == stored, (
        f"missing goldens: {sorted(frozen - stored)}; "
        f"orphaned goldens: {sorted(stored - frozen)}"
    )


def test_manifest_checksums_match_stored_goldens(frozen_cases, manifest):
    recorded = manifest["sha256"]
    assert set(recorded) == {c.name for c in frozen_cases}
    for name, digest in recorded.items():
        actual = hashlib.sha256(_golden_text(name).encode()).hexdigest()
        assert actual == digest, f"{name}: stored golden does not match its recorded sha256"


def test_golden_files_regenerate(oracle_exe, frozen_cases, manifest, rundir):
    """Rebuild the reference from source and reproduce every golden file."""
    strict = _same_platform(manifest)
    budget = 0 if strict else CROSS_PLATFORM_ULP

    worst_ulp, worst_where, byte_mismatches = 0.0, None, []
    for case in frozen_cases:
        produced = (rundir / "golden" / case.name / "case.csv")
        run_case(oracle_exe, case, rundir / "golden" / case.name)
        text = produced.read_text()
        expected_text = _golden_text(case.name)
        if text == expected_text:
            continue
        byte_mismatches.append(case.name)

        got = read_csv_output(produced)
        want = read_csv_output_from_text(expected_text, rundir / f"_{case.name}.csv")
        assert got.days.tolist() == want.days.tolist(), f"{case.name}: output days differ"
        diffs = ulp_diff(got.values, want.values)
        if diffs.max() > worst_ulp:
            row, col = np.unravel_index(int(diffs.argmax()), diffs.shape)
            worst_ulp, worst_where = float(diffs.max()), f"{case.name}:{want.columns[col]}"

    if strict:
        assert not byte_mismatches, (
            f"{len(byte_mismatches)} golden file(s) changed on the platform that "
            f"generated them: {byte_mismatches[:8]}"
        )
    else:
        assert worst_ulp <= budget, (
            f"worst cross-platform disagreement {worst_ulp:.0f} ULP at {worst_where}, "
            f"budget {budget}. Raise LUMPYREM_GOLDEN_ULP deliberately, or "
            f"regenerate the goldens on this platform."
        )
        if byte_mismatches:
            print(f"\n{len(byte_mismatches)}/{len(frozen_cases)} cases differ from the "
                  f"golden platform; worst {worst_ulp:.0f} ULP at {worst_where}")


def read_csv_output_from_text(text: str, scratch: Path):
    scratch.parent.mkdir(parents=True, exist_ok=True)
    scratch.write_text(text)
    return read_csv_output(scratch)


def test_goldens_carry_a_closed_water_balance(frozen_cases, rundir):
    """A golden file with an unexplained broken balance is a broken reference.

    Cases exposed to one of the two known reference defects are exempt; those
    are asserted directly in tests/test_oracle_selfcheck.py.
    """
    for case in frozen_cases:
        if balance_may_break(case):
            continue
        res = read_csv_output_from_text(_golden_text(case.name),
                                        rundir / f"bal_{case.name}.csv")
        throughput = max(float(np.abs(res.column("rainfall")).sum()), case.maxvol, 1e-30)
        err = float(np.abs(res.column("balance")).max()) / throughput
        assert err < 1e-12, f"{case.name}: golden balance error {err:.3e}"
