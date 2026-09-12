"""Does the golden set actually constrain each fidelity trap?

A gate that passes is only meaningful if it would fail when something breaks.
Each mutation below deliberately reintroduces one of the mistakes the plan
warns about; the golden set must reject it.

Two entries are marked as provably inert.  Those are not gaps in the test set
-- they are traps that turn out not to be traps, and the proof is recorded with
them so nobody re-adds the worry later.
"""

from __future__ import annotations

import importlib.util
import sys
from dataclasses import dataclass
from pathlib import Path

import pytest

from kernel_compare import golden_values, run_case_through
from ulp import ulp_diff

KERNEL_SOURCE = Path(__file__).resolve().parents[1] / "src" / "lumpyrem" / "core.py"


@dataclass(frozen=True)
class Mutation:
    trap: str
    what: str
    old: str
    new: str
    inert: str | None = None      # why it provably cannot change results


MUTATIONS = (
    Mutation(
        "trap 1", "fully implicit instead of Crank-Nicolson",
        "vd = (tvol + vol) * 0.5 / maxvol",
        "vd = tvol / maxvol",
    ),
    Mutation(
        "trap 3", "test convergence on the first iteration too",
        "                if iter_ != 1:\n                    if ircode == 0:",
        "                if True:\n                    if ircode == 0:",
    ),
    Mutation(
        "trap 4", "drop the irrigation convergence criterion",
        "                            if abs(dirrig - odirrig) <= tol * (dirrig + odirrig) * 0.5:\n"
        "                                break",
        "                            if True:\n                                break",
        inert="Once irrigation fires, tvol is pinned to vvol from the first "
              "iteration on, so every later iteration computes the same vd and "
              "returns bit-identical rtemp1/rtemp2/tempvol. The second criterion "
              "can force extra iterations -- it does so in 48 of the frozen "
              "cases -- but it cannot change the answer.",
    ),
    Mutation(
        "trap 5", "drop the store-empties fix-up",
        "            if tempvol == 0.0:\n                rtemp3 = rtemp1 + rtemp2",
        "            if False:\n                rtemp3 = rtemp1 + rtemp2",
    ),
    Mutation(
        "trap 6", "<= instead of < in the macropore split",
        "                if rtemp1 < mflowmax * tstep:",
        "                if rtemp1 <= mflowmax * tstep:",
        inert="On exact equality both branches add the same amount to "
              "macropore flow and add rtemp1 - mflowmax*tstep == 0.0 to runoff. "
              "The branches are numerically identical at the boundary, so the "
              "comparison operator cannot change any result.",
    ),
    Mutation(
        "trap 7", "release the whole tail element instead of (1 - frdelay)",
        "        f.recharge = f.recharge + (1.0 - frdelay) * drainsub[irdelay]",
        "        f.recharge = f.recharge + drainsub[irdelay]",
    ),
    Mutation(
        "trap 8", "floor the volume at a true 1e-10",
        "ELEVATION_FLOOR = 1.00000001335143196e-10",
        "ELEVATION_FLOOR = 1e-10",
    ),
    Mutation(
        "trap 9", "subtract lower-store evaporation from gw_pot_evap",
        "gw_pot_evap = f.potevapn - f.evapn",
        "gw_pot_evap = f.potevapn - f.evapn - f.evapn_br",
    ),
    Mutation(
        "trap 10", "expose the lower store's crop factor instead of hardwiring 1.0",
        "cropfac_br = 1.0            # trap 10: hardwired in the source",
        "cropfac_br = 0.9",
    ),
    Mutation(
        "trap 14", "repair recharge_for_br so the lower store gets claimed water",
        "            recharge_for_br = drainsub[irdelay]",
        "            recharge_for_br = recharge_for_br + drainsub[irdelay]",
    ),
)


def _build(mutation: Mutation | None, tmp_path: Path):
    source = KERNEL_SOURCE.read_text()
    name = "lumpyrem_core_pristine"
    if mutation is not None:
        assert source.count(mutation.old) == 1, (
            f"{mutation.trap}: anchor no longer matches the kernel; the mutation "
            f"test has drifted from the code it is meant to guard"
        )
        source = source.replace(mutation.old, mutation.new)
        name = "lumpyrem_core_" + mutation.trap.replace(" ", "_")
    path = tmp_path / f"{name}.py"
    path.write_text(source)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _worst_ulp(simulate, cases, rundir) -> float:
    worst = 0.0
    for case in cases:
        want = golden_values(case.name, rundir / f"m_{case.name}.csv")
        try:
            got = run_case_through(simulate, case)
        except Exception:
            return float("inf")
        if got.shape != want.shape:
            return float("inf")
        worst = max(worst, float(ulp_diff(got, want).max()))
    return worst


@pytest.mark.parametrize("mutation", MUTATIONS, ids=[m.trap for m in MUTATIONS])
def test_golden_set_rejects_each_mutation(mutation, frozen_cases, rundir, tmp_path):
    module = _build(mutation, tmp_path)
    worst = _worst_ulp(module.simulate, frozen_cases, rundir)

    if mutation.inert:
        assert worst == 0.0, (
            f"{mutation.trap} was documented as unable to change results, but "
            f"mutating it moved the output by {worst:.0f} ULP. The reasoning "
            f"recorded with it is wrong:\n  {mutation.inert}"
        )
        return

    assert worst > 10, (
        f"{mutation.trap} ({mutation.what}) changed nothing across "
        f"{len(frozen_cases)} golden cases. The gate does not constrain this "
        f"trap, so a port could get it wrong and still pass."
    )


def test_unmutated_kernel_is_the_one_under_test(frozen_cases, rundir, tmp_path):
    """Guard the harness: a fresh build of the real kernel must still be exact."""
    module = _build(None, tmp_path)
    assert _worst_ulp(module.simulate, frozen_cases[:30], rundir) == 0.0
