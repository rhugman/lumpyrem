"""lumpyrem.warmup: every kernel compiles, and a new process loads it from the cache."""

from __future__ import annotations

import pytest

from lumpyrem import warmup

pytest.importorskip("numba")


def test_compile_kernels_runs_every_entry_point():
    report = warmup.compile_kernels()
    assert report.fresh is None and report.cached is None
    assert report.second <= report.first + 1.0
    assert report.threads >= 1
    assert report.cache == warmup.cache_dir()


def test_fresh_process_loads_from_cache():
    report = warmup.compile_kernels(check=True)
    assert report.fresh is not None
    assert report.cached, str(report)


def test_main_reports_and_exits_zero(capsys):
    assert warmup.main(["--no-check"]) == 0
    out = capsys.readouterr().out
    assert "first calls" in out and "cache:" in out


def test_cache_dir_cannot_change_after_numba_import(tmp_path):
    import numba  # noqa: F401 -- the point is that it is already imported
    with pytest.raises(RuntimeError, match="NUMBA_CACHE_DIR"):
        warmup.compile_kernels(cache=tmp_path / "elsewhere")
