"""Does the golden set constrain the compiled kernel, trap by trap?

``test_kernel_mutations.py`` asks this of ``core``.  ``compiled`` is a second
copy of the same arithmetic, held so far only by the golden set and by exact
equality with ``core`` -- which would not notice a trap that neither the
goldens nor ``core`` exercise differently.  So the same mistakes are
reintroduced here, into ``compiled`` itself, and each must move its output.

It carries a second survey the core one cannot: the ring buffers.  They are
the one place ``compiled`` departs structurally from the Fortran, and their
correctness rests on an argument -- that summing the ring in logical order,
and claiming the initial buffers' tail once, reproduce the Fortran's sums
over all ``MAXDELAY`` elements.  Each way of getting that wrong is committed
below and must be rejected.

As in the core survey, "rejected" is measured against a fresh build of the
unmutated module on the same machine, and that build is checked against the
installed one, so the survey is about the mutation alone.  Every build is
compiled from source without Numba's cache.
"""

from __future__ import annotations

import importlib.util
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest

from kernel_compare import run_case_through
from ulp import ulp_diff

pytest.importorskip("numba", reason="numba not installed")

SOURCE = Path(__file__).resolve().parents[1] / "src" / "lumpyrem" / "compiled.py"

UNCACHED = ('_JIT = dict(cache=True, nogil=True, error_model="numpy")',
            '_JIT = dict(cache=False, nogil=True, error_model="numpy")')


@dataclass(frozen=True)
class Mutation:
    trap: str
    what: str
    old: str
    new: str
    inert: str | None = None


MUTATIONS = (
    Mutation(
        "trap 1", "fully implicit instead of Crank-Nicolson",
        "vd = (tvol + vol) * 0.5 / maxvol",
        "vd = tvol / maxvol",
    ),
    Mutation(
        "trap 2", "recompute the state after convergence instead of exiting with it",
        "                            converged = True              # trap 2\n"
        "                            break",
        "                            converged = True\n"
        "                            vd = (tvol + vol) * 0.5 / maxvol\n"
        "                            rtemp1 = drainage(vd, ks, m, l) * tstep\n"
        "                            break",
    ),
    Mutation(
        "trap 3", "test convergence on the first iteration too",
        "                if iter_ != 1:                            # trap 3",
        "                if True:",
    ),
    Mutation(
        "trap 4", "drop the irrigation convergence criterion",
        "                            if abs(dirrig - odirrig) <= tol * (dirrig + odirrig) * 0.5:",
        "                            if True:",
        inert="see tests/test_kernel_mutations.py: once irrigation fires, every "
              "later iterate is bit-identical, so exiting early changes nothing.",
    ),
    Mutation(
        "trap 5", "drop the store-empties fix-up",
        "            if tempvol == 0.0:                            # traps 5 and 13",
        "            if False:",
    ),
    Mutation(
        "trap 6", "<= instead of < in the macropore split",
        "                if rtemp1 < mflowmax * tstep:",
        "                if rtemp1 <= mflowmax * tstep:",
        inert="see tests/test_kernel_mutations.py: both branches add the same "
              "amounts at equality.",
    ),
    Mutation(
        "trap 7", "release the whole tail element instead of (1 - frdelay)",
        "        recharge = recharge + (1.0 - frdelay) * drainsub[td]",
        "        recharge = recharge + drainsub[td]",
    ),
    Mutation(
        "trap 8", "floor the volume at a true 1e-10",
        "    if dtemp < ELEVATION_FLOOR:\n        dtemp = ELEVATION_FLOOR",
        "    if dtemp < 1e-10:\n        dtemp = 1e-10",
    ),
    Mutation(
        "trap 9", "subtract lower-store evaporation from gw_pot_evap",
        "        row[22] = f[F_POTEVAPN] - f[F_EVAPN]              # trap 9",
        "        row[22] = f[F_POTEVAPN] - f[F_EVAPN] - f[F_EVAPN_BR]",
    ),
    Mutation(
        "trap 10", "expose the lower store's crop factor instead of hardwiring 1.0",
        "    cropfac_br = 1.0            # trap 10",
        "    cropfac_br = 0.9",
    ),
    Mutation(
        "trap 12", "delay drainage on its way out of a two-store model as well",
        "            allrecharge = f[F_MRECHARGE] + f[F_DRAINAGE_BR] + f[F_OVERFLOW_BR]",
        "            allrecharge = f[F_MRECHARGE] + f[F_RECHARGE] + f[F_OVERFLOW_BR]",
    ),
    Mutation(
        "trap 14", "repair recharge_for_br so the lower store gets claimed water",
        "            recharge_for_br = drainsub[td]",
        "            recharge_for_br = recharge_for_br + drainsub[td]",
    ),
    # -- The ring buffers.
    Mutation(
        "ring: order", "sum the ring oldest-first instead of in the Fortran's order",
        "    k = head\n    for _ in range(n):\n        total = total + ring[k]\n"
        "        k = k + 1\n        if k == n:\n            k = 0",
        "    k = head - 1 if head > 0 else n - 1\n    for _ in range(n):\n"
        "        total = total + ring[k]\n        k = k - 1\n        if k < 0:\n"
        "            k = n - 1",
    ),
    Mutation(
        "ring: claim", "drop the first call's sweep of the initial buffers",
        "    recharge = recharge + claim[0]\n",
        "",
    ),
    Mutation(
        "ring: claim twice", "release the sweep on every call, not just the first",
        "    claim[0] = 0.0\n",
        "",
    ),
    Mutation(
        "ring: day 0", "total the initial buffer over the ring only, not whole",
        "    for i in range(1, rbuf.shape[0]):\n        totd = totd + rbuf[i]",
        "    for i in range(1, min(rbuf.shape[0], drainsub.shape[0] + 1)):\n"
        "        totd = totd + rbuf[i]",
    ),
    Mutation(
        "ring: size", "size the ring ceil(delay) + 1, as the plan first sketched",
        "    drainsub = np.empty(int(p[RDELAY]) + 1)",
        "    drainsub = np.empty(int(np.ceil(p[RDELAY])) + 1)",
    ),
    Mutation(
        "ring: rotate", "forget to move the head, so the ring never shifts",
        "        hd = td\n",
        "",
    ),
    Mutation(
        "ring: persist", "restart each call with the head at slot 0",
        "    hd = heads[0]\n",
        "    hd = 0\n",
    ),
)


def _build(mutation: Mutation | None, tmp_path: Path):
    """Compile ``compiled.py``, mutated or not, as a module of the package."""
    source = SOURCE.read_text()
    assert source.count(UNCACHED[0]) == 1, "the JIT options line has changed"
    source = source.replace(*UNCACHED)
    name = "lumpyrem._compiled_pristine"
    if mutation is not None:
        assert source.count(mutation.old) == 1, (
            f"{mutation.trap}: anchor no longer matches compiled.py; the mutation "
            f"test has drifted from the code it is meant to guard"
        )
        source = source.replace(mutation.old, mutation.new)
        slug = "".join(ch if ch.isalnum() else "_" for ch in mutation.trap)
        name = f"lumpyrem._compiled_{slug}"
    path = tmp_path / f"{name.rsplit('.', 1)[1]}.py"
    path.write_text(source)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
    finally:
        del sys.modules[name]
    return module


def _outputs(module, cases) -> dict:
    out = {}
    for case in cases:
        try:
            out[case.name] = run_case_through(module.simulate, case)
        except Exception:           # a mutation may crash the kernel; that is a rejection
            out[case.name] = None
    return out


@pytest.fixture(scope="module")
def pristine(frozen_cases, tmp_path_factory) -> dict:
    module = _build(None, tmp_path_factory.mktemp("pristine"))
    return _outputs(module, frozen_cases)


def test_the_pristine_build_is_the_installed_kernel(frozen_cases, pristine):
    """Guard the harness: a fresh, uncached build is the kernel under test."""
    from lumpyrem import compiled
    for case in frozen_cases:
        assert np.array_equal(pristine[case.name],
                              run_case_through(compiled.simulate, case)), case.name


@pytest.mark.parametrize("mutation", MUTATIONS, ids=[m.trap for m in MUTATIONS])
def test_golden_set_rejects_each_mutation(mutation, frozen_cases, pristine, tmp_path):
    got = _outputs(_build(mutation, tmp_path), frozen_cases)
    worst = max(_worst_ulp(g, pristine[k]) for k, g in got.items())

    if mutation.inert:
        assert worst == 0.0, (
            f"{mutation.trap} was documented as unable to change results, but it "
            f"moved the output by {worst:.0f} ULP: {mutation.inert}")
        return
    assert worst > 10, (
        f"{mutation.trap} ({mutation.what}) changed nothing across "
        f"{len(frozen_cases)} golden cases, so compiled could get it wrong and pass")


def _worst_ulp(got, want) -> float:
    if got is None or got.shape != want.shape:
        return float("inf")
    return float(ulp_diff(got, want).max())
