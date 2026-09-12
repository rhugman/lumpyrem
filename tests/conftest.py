"""Shared fixtures: one compiled oracle and one frozen case set per session."""

from __future__ import annotations

import json
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


@pytest.fixture(scope="session")
def oracle_exe(tmp_path_factory) -> Path:
    if shutil.which("gfortran") is None:
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
