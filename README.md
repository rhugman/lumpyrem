# lumpyrem

A Python port of the [LUMPREM](https://pesthomepage.org) lumped-parameter
recharge model, adding array-based forcing and output, a modern API, and
compiled-speed execution.

**Status: Phase 3 in progress.** The reference oracle is in place, the model
kernel reproduces it bit-for-bit, and the Python API is built on top — real
dates, typed parameters, explicit resampling. A Numba build of the kernel is
bit-identical too, matches Fortran speed per cell and scales across cores
(`pip install -e ".[fast]"`, `scripts/benchmark_numba.py`); the array API
that uses it comes next. See
[docs/conversion-plan.md](docs/conversion-plan.md) for the full plan.

## Using it

```python
import lumpyrem as lr

forcing = lr.Forcing.from_csv("climate.csv")        # or a DataFrame, dict, arrays
model = lr.Model(
    upper=lr.UpperStore(maxvol=300.0, ks=3.0, m=0.4, l=0.5, mflowmax=1.5, rdelay=2.0),
    solver=lr.Solver(nstep=5),
)
results = model.run(forcing, times="MS")            # month-start output
results.df["total_rech"]                            # a pandas Series
results.balance_error()                             # how well the run closes
```

Forcing comes from a DataFrame, a NumPy array, a CSV or a dict, at whatever
frequency it was observed. Output times are a frequency string, a list of dates
or explicit day numbers. Row *k* of the results covers the interval *ending* at
its own timestamp, which is the convention a MODFLOW stress period uses.

A run that fails to converge warns instead of pretending; `Solver` takes
`on_nonconvergence="raise"` or `"ignore"` if you would rather it did something
else.

### Gap-filling is a choice, not a file format

`lumprem2.f` fills gaps in its forcing by three different rules, decided by
which file a value happened to be written in. The rules are sensible; their
invisibility is not. Here they are named, defaulted to the Fortran's choices,
and per-variable overridable:

| variable | default rule |
|---|---|
| `rainfall` | zero-fill — rain does not persist |
| `pot_evap`, `irrigate`, `gw_irrig_frac`, `pot_evap_lower` | forward-fill as steps |
| `crop_factor`, `veg_gamma` | linear interpolation |

```python
lr.Forcing.from_dataframe(df, fill={"rainfall": "forward"})
```

## What is here

`lumpyrem.core` is a faithful scalar port of LUMPREM2's `rechmod` kernel and its
simulation loop. It reproduces the Fortran **bit-for-bit — 0 ULP on every column
of all 209 reference cases**, against a gate that allows 10. Everything above it
— `parameters`, `forcing`, `model`, `results` — is designed rather than
translated, and is held to reproducing those same numbers exactly.

Fidelity is not just asserted. `tests/test_kernel_mutations.py` reintroduces
each known trap in turn and requires the golden set to reject it, which is what
found the two demotions and the one real coverage gap recorded in the plan; the
API layer is checked the same way, by mis-marshalling a run on purpose.

## The reference oracle

A regression oracle: the original Fortran, rebuilt from vendored source with a
widened output writer, driven over 209 frozen cases whose results are committed
as golden files. Every later phase is gated on reproducing them.

| | |
|---|---|
| `vendor/fortran/` | LUMPREM and LUMPREM2 sources, unmodified |
| `tests/oracle/` | Build, legacy `.in` writer, case generator, result readers |
| `tests/data/cases.json.gz` | 209 frozen case specifications |
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

Across platforms `exp` and `**` come from different maths libraries — glibc,
Apple's libm, mingw — whose last bit differs. Everything else in the kernel is
`+ - * /`, which IEEE 754 pins exactly. So the tests measure, on each machine,
whether Python and the Fortran share `exp` and `**`:

- **where they do** — macOS, Linux and Windows (UCRT64) in CI — the port is
  held to ≤ 10 ULP against the Fortran built on that machine, and reproduces it
  bit for bit on all three;
- **where they do not**, and for the committed goldens off the platform that
  made them, the check is a difference relative to the size of what each
  column is computed from, `LUMPYREM_CROSS_RTOL` = `1e-9`. The measured worst
  case is 7.7e-11, Windows against the macOS goldens.

`scripts/crossplatform_report.py` prints the full comparison for the machine
it runs on; CI runs it on all three platforms.

## Two defects in the reference

The golden files record what LUMPREM2 *does*, which in two situations is not
what it intends — its own water-balance column does not close. Both are pinned
by dedicated cases and documented in `tests/oracle/defects.py`. The port
**reproduces** them: the fidelity gate is bit-identity and both defects are in
the golden files, so a mass-conserving alternative belongs above the kernel as
an opt-in rather than silently inside it. `Results.balance_error()` is the
place to notice one.

## Licence

MIT. The vendored Fortran is the work of Watermark Numerical Computing and is
included for reference and testing.
