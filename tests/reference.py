"""What the port is compared against, and how closely it must agree.

The port is a bit-faithful translation, but bit-faithful *to what* depends on
the maths library.  ``+ - * /`` are pinned by IEEE 754; ``exp`` and ``**`` are
not, and their last bit differs between Apple's libm, glibc and mingw-w64.  The
first CI run measured what that means (docs/conversion-plan.md, section 7):

- Where the Python process and the Fortran reference call the *same* ``exp``
  and ``pow`` -- macOS, Linux -- the port reproduces a Fortran build on the
  same machine bit for bit.  That is the claim, and it stays exact.
- Where they do not -- the committed goldens off macOS, and Windows, where
  CPython uses the UCRT and mingw gfortran its own libm -- the two agree to
  within rounding, which ULP cannot express: a residual such as ``balance``
  sits at ~1e-15 and flips sign under a last-bit change, which is 1e18 ULP and
  physically nothing.  There the check is a difference relative to the scale
  of what each column is computed from (``SCALE_FROM``).

Which regime applies is measured, not assumed from the platform name:
``probe_maths`` evaluates ``exp`` and ``**`` in Fortran and in Python over the
range of arguments the kernel uses and compares the bits.
"""

from __future__ import annotations

import math
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ulp import ulp_diff

#: docs/conversion-plan.md, Phase 1 gate: exact regime.
ULP_BUDGET = 10

#: Tolerance regime: |a - b| <= CROSS_RTOL * scale.  Provisional; see the
#: cross-platform report in CI for the measured worst case.
CROSS_RTOL = float(os.environ.get("LUMPYREM_CROSS_RTOL", "1e-9"))

#: The scale a column's disagreement is measured against: the largest |value|
#: in the case of the columns it is computed from.  Most columns are measured
#: against themselves.  A column that is a difference is measured against the
#: terms it is the difference of, because it cannot be more precise than they
#: are -- ``balance`` against the fluxes it closes, ``del_vol_*`` against the
#: volume it is the change in.  One case-wide scale was tried first and is
#: wrong: it lets a small flux hide behind a large volume, so that swapping the
#: two delays -- which re-times all recharge -- passed as rounding.
SCALE_FROM: dict[str, tuple[str, ...]] = {
    "del_vol_upper": ("vol_upper", "del_vol_upper"),
    "del_vol_lower": ("vol_lower", "del_vol_lower"),
    "del_vol_drain": ("vol_drain", "del_vol_drain"),
    "del_vol_macro": ("vol_macro", "del_vol_macro"),
    "total_rech": ("drain_upper", "macro_upper", "drain_lower", "overflow_lower",
                   "total_rech"),
    "net_recharge": ("total_rech", "gw_withdrawal", "net_recharge"),
    "gw_pot_evap": ("pot_evap_upper", "evap_upper", "gw_pot_evap"),
    "balance": ("rainfall", "irrigation", "del_vol_upper", "del_vol_lower",
                "total_rech", "runoff", "evap_upper", "del_vol_drain",
                "del_vol_macro", "evap_lower", "balance"),
    "depth-to-water": ("elevation", "depth-to-water"),
}

_PROBE = """\
program probe
  implicit none
  integer :: n, i
  double precision :: x, b, e
  read(*,*) n
  do i = 1, n
    read(*,*) x, b, e
    write(*,'(z16.16,1x,z16.16)') transfer(exp(x), 0_8), transfer(b ** e, 0_8)
  end do
end program probe
"""


def probe_inputs() -> list[tuple[float, float, float]]:
    """(x, b, e) for exp(x) and b**e, shaped like the kernel's arguments:
    exp(-gamma*vd), vd**(1/m) and friends, and vol**power at elevation."""
    rng = np.random.default_rng(20260928)
    n = 4096
    exp_args = -rng.uniform(0.0, 60.0, n)                 # exp(-gamma * vd)
    bases = np.concatenate([rng.uniform(1e-6, 1.0, n // 2),       # vd, 1 - vd**(1/m)
                            10.0 ** rng.uniform(-10, 1.5, n // 2)])  # store volumes
    exponents = np.concatenate([rng.uniform(0.05, 20.0, n // 2),  # 1/m, m, l
                                rng.uniform(-1.0, 3.0, n // 2)])  # elevation power
    return [(float(x), float(b), float(e))
            for x, b, e in zip(exp_args, bases, exponents, strict=True)]


def probe_maths(fc: str, workdir: Path) -> tuple[bool, str]:
    """Do Fortran and Python get bit-identical ``exp`` and ``**`` here?"""
    workdir.mkdir(parents=True, exist_ok=True)
    src = workdir / "probe.f90"
    src.write_text(_PROBE)
    exe = workdir / ("probe.exe" if os.name == "nt" else "probe")
    subprocess.run([fc, "-O2", "-ffp-contract=off", "-o", str(exe), str(src)],
                   check=True, capture_output=True, cwd=workdir)
    pairs = probe_inputs()
    feed = f"{len(pairs)}\n" + "".join(f"{x!r} {b!r} {e!r}\n" for x, b, e in pairs)
    out = subprocess.run([str(exe)], input=feed, capture_output=True, text=True,
                         check=True).stdout.split()
    exp_bad = pow_bad = 0
    for i, (x, b, e) in enumerate(pairs):
        f_exp = np.uint64(int(out[2 * i], 16)).view(np.float64)
        f_pow = np.uint64(int(out[2 * i + 1], 16)).view(np.float64)
        exp_bad += f_exp != math.exp(x)
        pow_bad += f_pow != b ** e
    same = bool(exp_bad == 0 and pow_bad == 0)
    detail = (f"exp and ** bit-identical on {len(pairs)} probe arguments" if same else
              f"maths libraries differ: exp on {exp_bad}/{len(pairs)}, "
              f"** on {pow_bad}/{len(pairs)} probe arguments")
    return same, detail


def _scale_index() -> list[list[int]]:
    from lumpyrem.core import COLUMNS
    return [[COLUMNS.index(c) for c in SCALE_FROM.get(name, (name,))]
            for name in COLUMNS]


def scaled_diff(got: np.ndarray, want: np.ndarray) -> np.ndarray:
    """|got - want| relative to each column's scale in this case (``SCALE_FROM``).

    A NaN on one side only is infinite; NaN against NaN is agreement.
    """
    got = np.asarray(got, dtype=np.float64)
    want = np.asarray(want, dtype=np.float64)
    d = np.abs(got - want)
    d[np.isnan(got) & np.isnan(want)] = 0.0
    d[np.isnan(d)] = np.inf
    colmax = np.nan_to_num(np.nanmax(np.abs(want), axis=0), nan=0.0)
    tiny = np.finfo(np.float64).tiny
    scale = np.array([max(float(colmax[idx].max()), tiny) for idx in _scale_index()])
    with np.errstate(over="ignore"):    # nonzero where the reference is all zero
        return d / scale


@dataclass(frozen=True)
class Agreement:
    """How closely the port must match the reference on this machine."""

    exact: bool
    reason: str

    def worst(self, got: np.ndarray, want: np.ndarray) -> tuple[float, int]:
        """The largest disagreement and the column it is in."""
        d = ulp_diff(got, want) if self.exact else scaled_diff(got, want)
        _, col = np.unravel_index(int(d.argmax()), d.shape)
        return float(d.max()), int(col)

    def budget(self, ulp_budget: int = ULP_BUDGET) -> float:
        return float(ulp_budget) if self.exact else CROSS_RTOL

    def accepts(self, got, want, ulp_budget: int = ULP_BUDGET) -> bool:
        if np.shape(got) != np.shape(want):
            return False
        return self.worst(got, want)[0] <= self.budget(ulp_budget)

    def unit(self) -> str:
        return "ULP" if self.exact else "x case scale"

    def describe(self) -> str:
        rule = (f"exact regime, <= {ULP_BUDGET} ULP" if self.exact
                else f"tolerance regime, <= {CROSS_RTOL:g} of case scale")
        return f"{rule} ({self.reason})"
