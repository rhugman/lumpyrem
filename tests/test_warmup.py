"""lumpyrem.warmup: every kernel compiles, and a new process loads it from the cache.

Nothing here times anything: on a shared CI runner the clock says more about
the neighbours than about the cache.  Numba's own cache statistics say what
each kernel did; the timings are only reported.
"""

from __future__ import annotations

import pytest

from lumpyrem import warmup

pytest.importorskip("numba")


def test_fresh_process_loads_every_kernel_from_the_cache():
    report = warmup.compile_kernels(check=True)
    assert report.kernels, "the fresh process called no compiled kernel"
    assert not report.recompiled, str(report)
    assert report.cached


def test_main_reports_and_exits_zero(capsys):
    assert warmup.main(["--no-check"]) == 0
    out = capsys.readouterr().out
    assert "first calls" in out and "cache:" in out


def test_cache_dir_cannot_change_after_numba_import(tmp_path):
    import numba  # noqa: F401 -- the point is that it is already imported
    with pytest.raises(RuntimeError, match="NUMBA_CACHE_DIR"):
        warmup.compile_kernels(cache=tmp_path / "elsewhere")
