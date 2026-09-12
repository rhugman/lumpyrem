"""Exact ULP distance between float64 arrays, used to size fidelity budgets."""

from __future__ import annotations

import numpy as np

_INT64_MIN = np.int64(np.iinfo(np.int64).min)


def _ordered(x: np.ndarray) -> np.ndarray:
    """Map float64 bit patterns onto a monotonically increasing integer key."""
    i = np.asarray(x, dtype=np.float64).view(np.int64)
    return np.where(i >= 0, i, _INT64_MIN - i)


def ulp_diff(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Element-wise ULP distance.  NaN against NaN counts as zero."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    oa, ob = _ordered(a), _ordered(b)
    out = np.zeros(a.shape, dtype=np.float64)
    differing = oa != ob
    if differing.any():
        # Python ints for the differing elements only: int64 subtraction would
        # wrap silently when the two values sit far apart.
        flat = np.flatnonzero(differing.ravel())
        fa, fb = oa.ravel(), ob.ravel()
        out.ravel()[flat] = [abs(int(fa[i]) - int(fb[i])) for i in flat]
    both_nan = np.isnan(a) & np.isnan(b)
    out[both_nan] = 0.0
    return out
