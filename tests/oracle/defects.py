"""Known defects in the reference implementation.

The golden files record what LUMPREM2 actually does, which in two situations
is not what it intends: its own water-balance column does not close.  Both are
reproducible, both are pinned by dedicated cases, and both are recorded here
so that no test quietly asserts a closure the reference cannot deliver.

Phase 1 must decide, for each, whether the port reproduces the defect or
departs from it.  Reproducing keeps bit-fidelity; departing makes the model
conserve mass but breaks the gate on these cases.  Either way the decision
should be explicit rather than discovered.

IRRIGATION_BYPASSES_EMPTY_FIXUP
    rechmod2.f applies its store-empties fix-up only when ``tempvol == 0``,
    but irrigation raises ``tvol`` to ``vvol`` before ``tempvol`` is taken.
    When the negative-volume clamp fires and irrigation then tops the store
    up, the fix-up is skipped and the full pre-clamp demand is charged against
    water that was never available.  Mass is created; ``balance`` goes
    negative.  Pinned by ``irrig_bypasses_empty_fixup``.

LOWER_STORE_LOSES_CLAIMED_BUFFER_WATER
    ``iflag_br`` is initialised to 1 and is only cleared inside the branch
    that never runs, so ``recharge_for_br`` is overwritten each day instead of
    accumulating.  Water claimed by the first-pass sweep -- anything sitting
    in the recharge buffer beyond ``irdelay`` when a call begins -- is added to
    ``recharge`` but never routed into the lower store, and vanishes.  Mass is
    lost; ``balance`` goes positive by exactly that amount, once.  Pinned by
    ``two_store_buffer_claim_lost``.
"""

from __future__ import annotations

from .infile import Case


def irrigation_may_bypass_fixup(case: Case) -> bool:
    """Necessary conditions for the first defect (not sufficient: the
    negative-volume clamp still has to fire)."""
    return any(case.irrigcode) and case.irrigvolfrac > 0.0


def buffer_claim_may_be_lost(case: Case) -> float:
    """Recharge-buffer water beyond irdelay, which a two-store run will lose.

    Returns 0.0 when the defect cannot apply.  When it can, the return value
    is the exact amount the reference will fail to route into the lower store.
    """
    if not case.two_store:
        return 0.0
    irdelay = int(case.rdelay) + 1
    return float(sum(case.rbuf[irdelay:]))


def balance_may_break(case: Case) -> str | None:
    """Name the defect that can apply to this case, or None if neither can."""
    if irrigation_may_bypass_fixup(case):
        return "IRRIGATION_BYPASSES_EMPTY_FIXUP"
    if buffer_claim_may_be_lost(case) > 0.0:
        return "LOWER_STORE_LOSES_CLAIMED_BUFFER_WATER"
    return None
