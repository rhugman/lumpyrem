# lumpyrem

A Python port of the [LUMPREM](https://pesthomepage.org) lumped-parameter
recharge model, adding array-based forcing and output, a modern API, and
compiled-speed execution.

**Status: Phase 0 complete.** The reference oracle and its test scaffolding are
in place. The model kernel arrives in Phase 1. See
[docs/conversion-plan.md](docs/conversion-plan.md) for the full plan.

## What Phase 0 provides

A regression oracle: the original Fortran, rebuilt from vendored source with a
widened output writer, driven over 207 frozen cases whose results are committed
as golden files. Every later phase is gated on reproducing them.

| | |
|---|---|
| `vendor/fortran/` | LUMPREM and LUMPREM2 sources, unmodified |
| `tests/oracle/` | Build, legacy `.in` writer, case generator, result readers |
| `tests/data/cases.json.gz` | 207 frozen case specifications |
| `tests/data/golden/` | Reference output for each, full precision, gzipped |
| `scripts/generate_golden.py` | Rebuilds and refreezes the golden set |

The legacy `.in` writer exists only because the reference binary reads nothing
else. It lives under `tests/` and is never imported by the package.

## Getting started

Requires Python 3.10+ and `gfortran` (for the oracle only).

```bash
pip install -e ".[dev]"
pytest
```

To verify the golden files regenerate from a clean checkout:

```bash
python scripts/generate_golden.py --check
```

## Reproducibility

The Fortran is compiled with `-ffp-contract=off`, so the compiler cannot fuse
`a*b+c` into an FMA and produce different last bits on different architectures.
Golden files then reproduce **byte-for-byte** on the platform that generated
them.

Across platforms a small ULP budget applies, because `exp` and `**` come from
the system maths library and the last bit of a transcendental is not
standardised between glibc, Apple's libm and mingw. Everything else in the
kernel is `+ - * /`, which IEEE 754 pins exactly. The budget is
`LUMPYREM_GOLDEN_ULP`, currently 16 and provisional until CI has measured the
real figure on all three platforms.

## Two defects in the reference

The golden files record what LUMPREM2 *does*, which in two situations is not
what it intends — its own water-balance column does not close. Both are pinned
by dedicated cases and documented in `tests/oracle/defects.py`. Phase 1 must
decide explicitly whether the port reproduces them or departs from them.

## Licence

MIT. The vendored Fortran is the work of Watermark Numerical Computing and is
included for reference and testing.
