"""Typed parameter objects for LUMPREM2.

Each object owns one block of the model specification and validates itself on
construction.  Validation is deliberately *stricter* than the Fortran reader:
where the reader silently clamps a value into range -- ``gamma_br`` into
[0.1, 10], the delays into ``MAXDELAY - 2`` -- this layer refuses it instead,
because a silently clamped parameter is a calibration that quietly stops
responding to the knob being turned.

Every failure raises :class:`ParameterError`, which carries the owning class,
the field, the offending value and the requirement it broke, so a caller can
act on it programmatically rather than by matching on message text.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

from .core import MAXDELAY

__all__ = [
    "ParameterError",
    "UpperStore",
    "LowerStore",
    "Solver",
    "VolumeToElevation",
    "InitialState",
    "MAX_DELAY_DAYS",
    "LOWER_STORE_CROP_FACTOR",
]

#: The longest delay the kernel can carry.  ``lumprem2.f`` clamps to this; the
#: objects below reject instead.  Phase 3 replaces the fixed buffers with ring
#: buffers sized to the delay, which removes the ceiling entirely.
MAX_DELAY_DAYS = MAXDELAY - 2

#: The lower store's crop factor is hardwired to 1.0 in ``rechmod2.f`` -- marked
#: ``! hardwired`` there (trap 10).  It is surfaced as a named constant rather
#: than a parameter so that exposing it stays a deliberate decision.
LOWER_STORE_CROP_FACTOR = 1.0


class ParameterError(ValueError):
    """A parameter failed validation.

    Attributes:
        owner: name of the parameter object the field belongs to.
        field: the offending field.
        value: the value that was supplied.
        requirement: the constraint it broke, phrased to follow "must".
    """

    def __init__(self, owner: str, field: str, value: object, requirement: str):
        self.owner = owner
        self.field = field
        self.value = value
        self.requirement = requirement
        super().__init__(f"{owner}.{field} must {requirement}; got {value!r}")


def _require(ok: bool, owner: str, name: str, value: object, requirement: str) -> None:
    if not ok:
        raise ParameterError(owner, name, value, requirement)


def _as_float(owner: str, name: str, value: object) -> float:
    try:
        out = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        raise ParameterError(owner, name, value, "be a real number") from None
    _require(isfinite(out), owner, name, value, "be finite")
    return out


def _coerce_floats(obj, names: tuple[str, ...]) -> None:
    """Replace each named field with its float, in place, on a frozen dataclass."""
    owner = type(obj).__name__
    for name in names:
        object.__setattr__(obj, name, _as_float(owner, name, getattr(obj, name)))


def _check_buffer(owner: str, name: str, values) -> tuple[float, ...]:
    try:
        out = tuple(float(v) for v in values)
    except (TypeError, ValueError):
        raise ParameterError(owner, name, values, "be a sequence of real numbers") from None
    _require(len(out) <= MAXDELAY, owner, name, len(out),
             f"hold at most MAXDELAY = {MAXDELAY} entries")
    for i, v in enumerate(out):
        _require(isfinite(v) and v >= 0.0, owner, f"{name}[{i}]", v, "be finite and >= 0")
    return out


# ---------------------------------------------------------------------------
# Stores
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class UpperStore:
    """The root-zone store, and the two delays on the flows that leave it.

    ``ks``, ``m`` and ``l`` are the van Genuchten-Mualem drainage parameters of
    equation 2.2: drainage is ``ks * vd**l * (1 - (1 - vd**(1/m))**m)**2``.
    ``l`` keeps its single-letter name because that is what the literature and
    both manuals call the pore-connectivity exponent.

    Args:
        maxvol: storage at field capacity, in units of depth.
        ks: saturated drainage rate, depth per day.
        m: van Genuchten ``m``, in (0, 1).  ``1/m`` is taken, so it cannot be 0.
        l: Mualem pore-connectivity exponent; commonly 0.5, and may be negative
            in fitted soil-water models, so it is only required to be finite.
        mflowmax: cap on macropore flow, depth per day.  Overflow above the cap
            becomes runoff.
        irrigvolfrac: irrigation tops the store up to this fraction of
            ``maxvol`` on any day the forcing asks for it.
        rdelay: days of delay on drainage leaving the store.
        mdelay: days of delay on macropore flow leaving the store.
    """

    maxvol: float
    ks: float
    m: float
    l: float
    mflowmax: float = 0.0
    irrigvolfrac: float = 0.0
    rdelay: float = 0.0
    mdelay: float = 0.0

    def __post_init__(self) -> None:
        _coerce_floats(self, ("maxvol", "ks", "m", "l", "mflowmax",
                              "irrigvolfrac", "rdelay", "mdelay"))
        o = "UpperStore"
        _require(self.maxvol > 0.0, o, "maxvol", self.maxvol, "be > 0")
        _require(self.ks >= 0.0, o, "ks", self.ks, "be >= 0")
        _require(self.m > 0.0, o, "m", self.m, "be > 0 (the drainage law takes 1/m)")
        _require(self.mflowmax >= 0.0, o, "mflowmax", self.mflowmax, "be >= 0")
        _require(0.0 <= self.irrigvolfrac <= 1.0, o, "irrigvolfrac",
                 self.irrigvolfrac, "be in [0, 1]")
        for name in ("rdelay", "mdelay"):
            v = getattr(self, name)
            _require(0.0 <= v <= MAX_DELAY_DAYS, o, name, v,
                     f"be in [0, {MAX_DELAY_DAYS}] "
                     f"(lumprem2.f clamps silently at MAXDELAY - 2)")


@dataclass(frozen=True)
class LowerStore:
    """The optional store below the root zone.

    Constructing one activates it: ``maxvol`` must be positive, because
    ``maxvol_br <= 0`` is exactly how the Fortran disables the second store.
    Pass ``lower=None`` to the model for the one-store configuration.

    Enabling it re-routes the delays (trap 12): upper-store drainage is delayed
    on its way *into* this store, and only macropore flow is delayed on its way
    out to the water table.

    Args:
        maxvol: storage at field capacity.
        gamma: extraction curvature.  ``lumprem2.f`` clamps this into
            [0.1, 10] without saying so; values outside are rejected here.
        ks: saturated drainage rate, depth per day.
        m: van Genuchten ``m``, > 0.
        l: pore-connectivity exponent.
        extravol: storage above field capacity.  The store overflows at
            ``maxvol + extravol``; relative saturation is still measured
            against ``maxvol`` alone.
    """

    maxvol: float
    gamma: float
    ks: float
    m: float
    l: float
    extravol: float = 0.0

    def __post_init__(self) -> None:
        _coerce_floats(self, ("maxvol", "gamma", "ks", "m", "l", "extravol"))
        o = "LowerStore"
        _require(self.maxvol > 0.0, o, "maxvol", self.maxvol,
                 "be > 0 (a lower store with maxvol <= 0 is a disabled store; "
                 "pass lower=None instead)")
        _require(0.1 <= self.gamma <= 10.0, o, "gamma", self.gamma,
                 "be in [0.1, 10] (lumprem2.f clamps silently)")
        _require(self.ks >= 0.0, o, "ks", self.ks, "be >= 0")
        _require(self.m > 0.0, o, "m", self.m, "be > 0 (the drainage law takes 1/m)")
        _require(self.extravol >= 0.0, o, "extravol", self.extravol, "be >= 0")

    @property
    def capacity(self) -> float:
        """Volume at which the store overflows, ``maxvol + extravol``."""
        return self.maxvol + self.extravol


# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Solver:
    """Sub-stepping and Picard-iteration settings.

    Raising ``nstep`` does not always help.  The drainage law has an unbounded
    slope as the store approaches full whenever ``m < 1``, so the iteration is
    not formally contractive there, and a saturated store with a small ``m``
    can converge *worse* at finer sub-steps -- more sub-steps spend more time
    pinned against the singularity.

    Args:
        nstep: sub-steps per day.
        mxiter: Picard iterations allowed per sub-step.
        tol: convergence tolerance, relative to ``maxvol``.
        on_nonconvergence: what to do when a sub-step exhausts ``mxiter``.
            ``"warn"`` (the default) issues a :class:`ConvergenceWarning` and
            returns the results; ``"raise"`` raises; ``"ignore"`` stays silent.
            The Fortran prints to stdout and carries on regardless, which is
            how a non-converged answer gets mistaken for a converged one.
    """

    nstep: int = 5
    mxiter: int = 100
    tol: float = 1.0e-6
    on_nonconvergence: str = "warn"

    def __post_init__(self) -> None:
        o = "Solver"
        for name in ("nstep", "mxiter"):
            v = getattr(self, name)
            _require(isinstance(v, int) and not isinstance(v, bool) and v >= 1,
                     o, name, v, "be an integer >= 1")
        _coerce_floats(self, ("tol",))
        _require(self.tol > 0.0, o, "tol", self.tol, "be > 0")
        _require(self.on_nonconvergence in ("warn", "raise", "ignore"),
                 o, "on_nonconvergence", self.on_nonconvergence,
                 "be one of 'warn', 'raise', 'ignore'")


# ---------------------------------------------------------------------------
# Volume to elevation
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class VolumeToElevation:
    """Converts a store's volume to a water-table elevation, and its clamps.

    ``elevation = offset + factor1 * v + factor2 * v ** power``, where ``v`` is
    the chosen store's volume capped at that store's ``maxvol`` and floored
    just above zero so a negative ``power`` cannot blow up.  The result is then
    clipped into ``[elevmin, elevmax]``.

    Args:
        offset, factor1, factor2, power: coefficients of the conversion.
        datum: topographic surface.  ``depth_to_water = datum - elevation``.
        bucket: which store drives the conversion, ``"upper"`` or ``"lower"``.
            ``"lower"`` requires a lower store.
        elevmin, elevmax: clips on the result; ``None`` means unclipped.
    """

    offset: float
    factor1: float
    factor2: float
    power: float
    datum: float
    bucket: str = "upper"
    elevmin: float | None = None
    elevmax: float | None = None

    def __post_init__(self) -> None:
        _coerce_floats(self, ("offset", "factor1", "factor2", "power", "datum"))
        o = "VolumeToElevation"
        _require(self.bucket in ("upper", "lower"), o, "bucket", self.bucket,
                 "be 'upper' or 'lower'")
        for name in ("elevmin", "elevmax"):
            v = getattr(self, name)
            if v is not None:
                object.__setattr__(self, name, _as_float(o, name, v))
        if self.elevmin is not None and self.elevmax is not None:
            _require(self.elevmin <= self.elevmax, o, "elevmax", self.elevmax,
                     f"be >= elevmin ({self.elevmin!r})")


# ---------------------------------------------------------------------------
# Initial conditions
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class InitialState:
    """Store volumes and delay-buffer contents at the start of the run.

    The buffers hold water already on its way out of the upper store but not
    yet released: ``drain_buffer[0]`` is the most recent day's drainage,
    ``drain_buffer[k]`` the drainage from ``k`` days before the run began.
    Anything sitting beyond the delay length is released on the first call.

    Args:
        vol: upper-store volume, in [0, ``UpperStore.maxvol``].
        vol_lower: lower-store volume, in [0, ``LowerStore.capacity``].
        drain_buffer: pending drainage, most recent first.
        macro_buffer: pending macropore flow, most recent first.
    """

    vol: float = 0.0
    vol_lower: float = 0.0
    drain_buffer: tuple[float, ...] = (0.0,)
    macro_buffer: tuple[float, ...] = (0.0,)

    def __post_init__(self) -> None:
        _coerce_floats(self, ("vol", "vol_lower"))
        o = "InitialState"
        _require(self.vol >= 0.0, o, "vol", self.vol, "be >= 0")
        _require(self.vol_lower >= 0.0, o, "vol_lower", self.vol_lower, "be >= 0")
        for name in ("drain_buffer", "macro_buffer"):
            buf = _check_buffer(o, name, getattr(self, name))
            _require(len(buf) >= 1, o, name, buf, "hold at least one entry")
            object.__setattr__(self, name, buf)
