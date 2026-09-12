"""Compile the vendored LUMPREM Fortran as a full-precision reference oracle.

The shipped binaries are Windows-only and the tabular writer emits 7
significant figures, which is far too coarse to gate a port against.  This
module rebuilds from the vendored source with two surgical edits that widen
the CSV writer to 17 significant figures, then compiles with gfortran.

Nothing here is part of the ``lumpyrem`` package -- it exists only to drive
the reference implementation during testing.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
VENDOR_DIR = REPO_ROOT / "vendor" / "fortran"

# Strict IEEE semantics per operation.  Without -ffp-contract=off gfortran is
# free to fuse a*b+c into an FMA, which lands differently on arm64 and x86 and
# would make the golden files platform-dependent.
FFLAGS = ("-O2", "-std=legacy", "-ffp-contract=off", "-fno-fast-math")


@dataclass(frozen=True)
class Patch:
    """A single exactly-once source edit, so a silent miss cannot happen."""

    old: str
    new: str
    why: str


# The CSV writer is widened rather than the fixed-width tabular writer, so the
# tabular output stays byte-comparable with the shipped 2020 example files.
#
# 1pe26.17 is 18 significant figures in a guaranteed round-trippable form;
# 1pg14.7 was 7.  26 fields at width 27 overflow the 600-character line
# buffer the writer formats into, so the buffer is widened to match.
PATCHES_LUMPREM2 = (
    Patch(
        old="        character*600 cline\n",
        new="        character*1200 cline\n",
        why="widen line buffer to hold 26 full-precision CSV fields",
    ),
    Patch(
        old="1262      format(i6,26(',',1pg14.7))\n",
        new="1262      format(i6,26(',',1pe26.17))\n",
        why="CSV writer: 7 -> 18 significant figures",
    ),
    # Bug fix, not a precision change.  lumprem2.f blanks the keyword "upper"
    # out of the line buffer before parsing ELEVMIN, but the matching
    # statement for "lower" sits after an unconditional `go to 9998` and never
    # runs.  The word survives, is parsed as ELEVMIN, and every run that asks
    # for lower-store elevation dies with "Cannot read value for ELEVMIN".
    # Reviving the dead statement is unambiguously what the author intended,
    # and it cannot change any input that works today -- the only inputs it
    # affects are ones that currently abort.  Without it the oracle could not
    # cover the nbucket=2 branch of the volume-to-elevation conversion.
    Patch(
        old="          nbucket=2\n        end if\n        if(nbucket.eq.0)nbucket=1\n",
        new="          nbucket=2\n          cline(nn:nn+4)=' '\n        end if\n"
            "        if(nbucket.eq.0)nbucket=1\n",
        why="revive dead statement so bucket='lower' is reachable",
    ),
)

# LUMPREM v1 has no CSV writer, so there is nothing to widen.  It is built
# unpatched and used only to cross-check the stale shipped example output; the
# golden files all come from v2.
PATCHES_LUMPREM1: tuple[Patch, ...] = ()


def _apply(text: str, patches: tuple[Patch, ...], name: str) -> str:
    for patch in patches:
        count = text.count(patch.old)
        if count != 1:
            raise RuntimeError(
                f"{name}: expected exactly 1 occurrence of {patch.old!r} "
                f"({patch.why}), found {count}. The vendored source has "
                f"changed; update tests/oracle/build.py."
            )
        text = text.replace(patch.old, patch.new)
    return text


def gfortran() -> str:
    exe = os.environ.get("FC") or shutil.which("gfortran")
    if exe is None:
        raise RuntimeError(
            "gfortran not found. Install it (brew install gcc / apt install "
            "gfortran / choco install mingw) or set $FC."
        )
    return exe


def build(build_dir: Path, version: int = 2) -> Path:
    """Patch and compile LUMPREM ``version``; return the executable path."""
    if version == 2:
        driver, kernel, patches = "lumprem2.f", "rechmod2.f", PATCHES_LUMPREM2
    elif version == 1:
        driver, kernel, patches = "lumprem.f", "rechmod.f", PATCHES_LUMPREM1
    else:
        raise ValueError(f"version must be 1 or 2, got {version}")

    build_dir = Path(build_dir)
    src_dir = build_dir / f"src_v{version}"
    src_dir.mkdir(parents=True, exist_ok=True)

    driver_text = _apply((VENDOR_DIR / driver).read_text(), patches, driver)
    (src_dir / driver).write_text(driver_text)
    shutil.copy2(VENDOR_DIR / kernel, src_dir / kernel)

    exe = build_dir / (f"lumprem{version}.exe" if sys.platform == "win32" else f"lumprem{version}")
    if exe.exists():
        exe.unlink()

    result = subprocess.run(
        [gfortran(), *FFLAGS, "-o", str(exe), driver, kernel],
        cwd=src_dir,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0 or not exe.exists():
        raise RuntimeError(
            f"gfortran failed for LUMPREM v{version}:\n{result.stdout}\n{result.stderr}"
        )
    return exe
