"""The Phase 0 gate: golden files must regenerate reproducibly.

On the platform that produced them the requirement is byte-identity.  Anywhere
else ``exp`` and ``**`` come from a different maths library, whose last bit
differs, so the requirement is the tolerance regime of tests/reference.py: a
difference relative to the size of the water flows in the case.  ULP was the
first guess and the first CI run retired it -- a near-zero residual such as
``balance`` flips sign under a last-bit change, which is 1e18 ULP and
physically nothing.  Everything else in the kernel is +, -, * and /, which
IEEE 754 pins exactly and which -ffp-contract=off keeps un-fused.
"""

from __future__ import annotations

import gzip
import hashlib
import platform
from pathlib import Path

import numpy as np

from oracle.defects import balance_may_break
from oracle.results import read_csv_output
from reference import CROSS_RTOL, scaled_diff
from ulp import ulp_diff

GOLDEN_DIR = Path(__file__).resolve().parent / "data" / "golden"


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


def test_golden_files_regenerate(local_runs, frozen_cases, manifest, rundir):
    """Rebuild the reference from source and reproduce every golden file."""
    strict = _same_platform(manifest)

    worst, worst_where, worst_ulp, byte_mismatches = 0.0, None, 0.0, []
    for case in frozen_cases:
        produced = local_runs[case.name]
        text = produced.read_text()
        expected_text = _golden_text(case.name)
        if text == expected_text:
            continue
        byte_mismatches.append(case.name)

        got = read_csv_output(produced)
        want = read_csv_output_from_text(expected_text, rundir / f"_{case.name}.csv")
        assert got.days.tolist() == want.days.tolist(), f"{case.name}: output days differ"
        worst_ulp = max(worst_ulp, float(ulp_diff(got.values, want.values).max()))
        diffs = scaled_diff(got.values, want.values)
        if diffs.max() > worst:
            _, col = np.unravel_index(int(diffs.argmax()), diffs.shape)
            worst, worst_where = float(diffs.max()), f"{case.name}:{want.columns[col]}"

    if strict:
        assert not byte_mismatches, (
            f"{len(byte_mismatches)} golden file(s) changed on the platform that "
            f"generated them: {byte_mismatches[:8]}"
        )
    else:
        assert worst <= CROSS_RTOL, (
            f"worst cross-platform disagreement {worst:.3g} of case scale at "
            f"{worst_where}, tolerance {CROSS_RTOL:g}. Raise LUMPYREM_CROSS_RTOL "
            f"deliberately, or regenerate the goldens on this platform."
        )
        if byte_mismatches:
            print(f"\n{len(byte_mismatches)}/{len(frozen_cases)} cases differ from the "
                  f"golden platform; worst {worst:.3g} of case scale at {worst_where} "
                  f"({worst_ulp:.3g} ULP)")


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
