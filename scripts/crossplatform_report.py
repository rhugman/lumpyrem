"""Measure how the reference and the port disagree on *this* platform.

The golden files were generated on macOS/arm64.  Elsewhere the last bit of
``exp`` and ``**`` differs, and a ULP budget was set as a guess until CI could
measure the real figure.  This script is the measurement.  For every frozen
case it produces four outputs --

    G  the committed golden file
    F  the Fortran reference, rebuilt and rerun here
    P  lumpyrem.core, run here
    C  lumpyrem.compiled, run here (if Numba is installed)

-- and reports, per comparison and per column, how far apart they are: in ULP,
absolutely, and relative to the size of that column in that case.  The last
one matters because a residual such as ``balance`` sits at ~1e-17 and flips
sign under a last-bit change, which is 1e18 ULP and physically nothing.

    python scripts/crossplatform_report.py [--json report.json]

Writes Markdown to stdout, and to the GitHub step summary when run in Actions.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO_ROOT / "src"), str(REPO_ROOT / "tests")]

from kernel_compare import golden_values, run_case_through  # noqa: E402
from lumpyrem import core  # noqa: E402
from oracle.build import build, gfortran  # noqa: E402
from oracle.infile import load_cases  # noqa: E402
from oracle.runner import run_case  # noqa: E402
from ulp import ulp_diff  # noqa: E402

try:
    from lumpyrem import compiled
except ImportError:
    compiled = None

COMPARISONS = (
    ("F vs G", "local Fortran against the committed goldens", "F", "G"),
    ("P vs F", "Python kernel against local Fortran", "P", "F"),
    ("P vs G", "Python kernel against the committed goldens", "P", "G"),
    ("C vs P", "compiled kernel against the Python kernel", "C", "P"),
)
BUDGET = 16


def _scale(values: np.ndarray) -> np.ndarray:
    """Per-column magnitude of each case: the largest |value| it takes."""
    return np.maximum(np.abs(values).max(axis=0), np.finfo(float).tiny)


def compare(a: np.ndarray, b: np.ndarray) -> dict:
    ulp = ulp_diff(a, b)
    absd = np.abs(a - b)
    absd[np.isnan(a) & np.isnan(b)] = 0.0
    rel = absd / _scale(b)
    return dict(ulp=ulp, abs=absd, rel=rel)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--json", type=Path)
    args = ap.parse_args()

    cases = load_cases(REPO_ROOT / "tests" / "data" / "cases.json.gz")
    stats = {name: {"exact": 0, "cases": 0, "cols": {}} for name, *_ in COMPARISONS}
    worst_rows = {name: [] for name, *_ in COMPARISONS}

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        exe = build(tmp / "build", version=2)
        for case in cases:
            out = {
                "G": golden_values(case.name, tmp / "g" / f"{case.name}.csv"),
                "F": run_case(exe, case, tmp / "run" / case.name).results.values,
                "P": run_case_through(core.simulate, case),
            }
            if compiled is not None:
                out["C"] = run_case_through(compiled.simulate, case)
            for name, _, x, y in COMPARISONS:
                if x not in out or y not in out:
                    continue
                if out[x].shape != out[y].shape:
                    raise SystemExit(f"{case.name}: {x} {out[x].shape} vs {y} {out[y].shape}")
                d = compare(out[x], out[y])
                s = stats[name]
                s["cases"] += 1
                s["exact"] += int(d["ulp"].max() == 0.0)
                for j, col in enumerate(core.COLUMNS):
                    c = s["cols"].setdefault(col, dict(ulp=0.0, abs=0.0, rel=0.0,
                                                       over=0, cells=0))
                    c["ulp"] = max(c["ulp"], float(d["ulp"][:, j].max()))
                    c["abs"] = max(c["abs"], float(d["abs"][:, j].max()))
                    c["rel"] = max(c["rel"], float(d["rel"][:, j].max()))
                    c["over"] += int((d["ulp"][:, j] > BUDGET).sum())
                    c["cells"] += d["ulp"].shape[0]
                i, j = np.unravel_index(int(d["rel"].argmax()), d["rel"].shape)
                worst_rows[name].append((float(d["rel"][i, j]), case.name,
                                         core.COLUMNS[j], int(i),
                                         float(out[x][i, j]), float(out[y][i, j])))

    lines = [
        f"## Cross-platform report: {platform.system()} {platform.machine()}",
        "",
        f"`{platform.platform()}`, Python {platform.python_version()}, "
        f"numpy {np.__version__}, compiler `{gfortran()}`, "
        f"compiled kernel {'on' if compiled else 'off'}",
        "",
    ]
    for name, what, *_ in COMPARISONS:
        s = stats[name]
        if not s["cases"]:
            continue
        lines += [f"### {name}: {what}", "",
                  f"{s['exact']}/{s['cases']} cases bit-identical.", ""]
        differing = {k: v for k, v in s["cols"].items() if v["ulp"] > 0}
        if differing:
            lines += ["| column | max ULP | cells > 16 ULP | max abs | max rel to column |",
                      "|---|---:|---:|---:|---:|"]
            for col, v in sorted(differing.items(), key=lambda kv: -kv[1]["rel"]):
                lines.append(f"| {col} | {v['ulp']:.3g} | {v['over']}/{v['cells']} | "
                             f"{v['abs']:.3g} | {v['rel']:.3g} |")
            lines += ["", "Worst cases, relative to column magnitude:", ""]
            for rel, cname, col, row, a, b in sorted(worst_rows[name], reverse=True)[:5]:
                lines.append(f"- `{cname}` `{col}` row {row}: {a!r} vs {b!r} "
                             f"(rel {rel:.3g})")
            lines.append("")

    text = "\n".join(lines)
    print(text)
    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary, "a", encoding="utf-8") as fh:
            fh.write(text + "\n")
    if args.json:
        args.json.write_text(json.dumps(
            {"platform": platform.platform(), "stats": stats}, indent=1))


if __name__ == "__main__":
    main()
