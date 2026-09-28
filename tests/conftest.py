"""Shared fixtures: one compiled oracle and one frozen case set per session."""

from __future__ import annotations

import json
import os
import shutil
import tempfile
from pathlib import Path

import pytest

from oracle.build import build
from oracle.infile import load_cases

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "tests" / "data"


def pytest_configure(config):
    config.addinivalue_line("markers", "oracle: requires gfortran to build the reference")


def _have_fortran() -> bool:
    """The same test oracle.build.gfortran() applies: $FC, else gfortran on PATH."""
    return bool(os.environ.get("FC")) or shutil.which("gfortran") is not None


@pytest.fixture(scope="session")
def oracle_exe(tmp_path_factory) -> Path:
    if not _have_fortran():
        pytest.skip("gfortran not available")
    return build(tmp_path_factory.mktemp("oracle_build"), version=2)


@pytest.fixture(scope="session")
def frozen_cases():
    if not (DATA_DIR / "cases.json.gz").exists():
        pytest.skip("no frozen cases; run scripts/generate_golden.py")
    return load_cases(DATA_DIR / "cases.json.gz")


@pytest.fixture(scope="session")
def manifest() -> dict:
    path = DATA_DIR / "manifest.json"
    if not path.exists():
        pytest.skip("no manifest; run scripts/generate_golden.py")
    return json.loads(path.read_text())


@pytest.fixture(scope="session")
def rundir():
    with tempfile.TemporaryDirectory() as tmp:
        yield Path(tmp)


# ---------------------------------------------------------------------------
# What the port is compared against; see tests/reference.py.
# ---------------------------------------------------------------------------

def _golden_platform(manifest: dict) -> bool:
    import platform
    tc = manifest["toolchain"]
    return (tc["machine"] == platform.machine()
            and tc["platform"].split("-")[0] == platform.platform().split("-")[0])


@pytest.fixture(scope="session")
def agreement(tmp_path_factory, manifest):
    """Exact where Python and the Fortran share exp and **; a tolerance elsewhere."""
    from oracle.build import gfortran
    from reference import Agreement, probe_maths
    if not _have_fortran():
        if _golden_platform(manifest):
            return Agreement(True, "no gfortran; comparing against goldens on "
                                   "the platform that generated them")
        pytest.skip("gfortran not available and not on the golden platform")
    same, detail = probe_maths(gfortran(), tmp_path_factory.mktemp("probe"))
    return Agreement(same, detail)


@pytest.fixture(scope="session")
def local_runs(oracle_exe, frozen_cases, tmp_path_factory) -> dict:
    """Every frozen case run through the reference built here, once per session."""
    from oracle.runner import run_case
    root = tmp_path_factory.mktemp("reference")
    runs = {}
    for case in frozen_cases:
        run_case(oracle_exe, case, root / case.name)
        runs[case.name] = root / case.name / "case.csv"
    return runs


@pytest.fixture(scope="session")
def reference(frozen_cases, manifest, tmp_path_factory, request) -> dict:
    """Reference output per case: the Fortran built on this machine.

    Without gfortran, on the platform that generated them, the committed
    goldens stand in -- they are the same numbers there, byte for byte.
    """
    from oracle.results import read_csv_output
    if not _have_fortran():
        if not _golden_platform(manifest):
            pytest.skip("gfortran not available and not on the golden platform")
        from kernel_compare import golden_values
        scratch = tmp_path_factory.mktemp("goldens")
        return {c.name: golden_values(c.name, scratch / f"{c.name}.csv")
                for c in frozen_cases}
    runs = request.getfixturevalue("local_runs")
    return {name: read_csv_output(path).values for name, path in runs.items()}
