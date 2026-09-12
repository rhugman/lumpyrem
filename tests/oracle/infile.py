"""Legacy LUMPREM2 input-file writer -- a test fixture, never an API.

The reference binary accepts nothing but its own fixed ``.in`` control file
and day-indexed forcing files, so driving it requires a writer for that
format.  This is the only legacy-format code in the project and it is
deliberately confined to ``tests/``; ``lumpyrem`` itself never imports it.

Forcing is always written densely -- one row per simulated day -- so that the
three gap-filling rules in the Fortran reader (linear interpolation for
vegetation, forward-fill for evaporation and irrigation, zero-fill for
rainfall) are no-ops and cannot contaminate a fidelity comparison.
"""

from __future__ import annotations

import gzip
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

# lumprem2.f: parameter(MAXDELAY=500); rdelay/mdelay are clamped to MAXDELAY-2.
MAXDELAY = 500

# The 26 output columns, in the order lumprem2.f writes them.
COLUMNS = (
    "vol_upper", "vol_lower", "vol_drain", "vol_macro",
    "del_vol_upper", "del_vol_lower", "del_vol_drain", "del_vol_macro",
    "rainfall", "irrigation", "drain_upper", "macro_upper",
    "drain_lower", "overflow_lower", "total_rech", "gw_withdrawal",
    "net_recharge", "runoff", "pot_evap_upper", "evap_upper",
    "pot_evap_lower", "evap_lower", "gw_pot_evap", "balance",
    "elevation", "depth-to-water",
)


def _r(x: float) -> str:
    """Format a real at full float64 precision so the round trip is exact."""
    return f"{float(x):.17g}"


@dataclass
class Case:
    """One fully specified LUMPREM2 run, with dense daily forcing."""

    name: str

    # -- earth properties, upper store
    maxvol: float
    irrigvolfrac: float
    rdelay: float
    mdelay: float
    ks: float
    m: float
    l: float
    mflowmax: float

    # -- earth properties, lower store (maxvol_br == 0 disables it)
    maxvol_br: float
    extravol_br: float
    gamma_br: float
    ks_br: float
    m_br: float
    l_br: float

    # -- volume to elevation
    offset: float
    factor1: float
    factor2: float
    power: float
    bucket: str            # "upper" or "lower"
    elevmin: float | None
    elevmax: float | None
    datum: float

    # -- initial conditions
    vol: float
    vol_br: float
    rbuf: list[float]
    mbuf: list[float]

    # -- solution parameters
    nstep: int
    mxiter: int
    tol: float

    # -- timing
    numdays: int
    outdays: list[int]

    # -- dense daily forcing, each of length numdays
    cropfac: list[float]
    gamma: list[float]
    rain: list[float]
    epot: list[float]
    irrigcode: list[int]
    gwirrigfrac: list[float]
    epot_br: list[float] = field(default_factory=list)

    # LUMPREM2 warns and carries on when the Picard loop hits MXITER, exiting
    # with whatever iterate it holds.  That is reproducible and is the sharpest
    # available test of traps 2 and 3 (the loop exits carrying the current
    # iterate; the first iteration never tests convergence), but it is only
    # ever wanted deliberately -- everywhere else it means the case is junk.
    allow_nonconvergence: bool = False

    # ------------------------------------------------------------------
    @property
    def two_store(self) -> bool:
        return self.maxvol_br > 0.0

    def validate(self) -> None:
        """Mirror the Fortran reader's own checks, so a bad case fails here."""
        n = self.numdays
        if self.maxvol <= 0:
            raise ValueError("maxvol must be > 0")
        if not 0.0 <= self.irrigvolfrac <= 1.0:
            raise ValueError("irrigvolfrac must be in [0, 1]")
        if not 0.0 <= self.vol <= self.maxvol:
            raise ValueError("vol must be in [0, maxvol]")
        if self.rdelay > MAXDELAY - 2 or self.mdelay > MAXDELAY - 2:
            raise ValueError("rdelay/mdelay exceed MAXDELAY-2; Fortran clamps silently")
        if self.rdelay < 0 or self.mdelay < 0:
            raise ValueError("rdelay/mdelay must be >= 0")
        if not self.rbuf or not self.mbuf:
            raise ValueError("NRBUF and NMBUF must both be > 0")
        if len(self.rbuf) > MAXDELAY or len(self.mbuf) > MAXDELAY:
            raise ValueError("buffer longer than MAXDELAY")
        if any(v < 0 for v in self.rbuf) or any(v < 0 for v in self.mbuf):
            raise ValueError("buffer entries must be >= 0")
        if self.nstep <= 0 or self.mxiter <= 0 or self.tol <= 0:
            raise ValueError("nstep, mxiter and tol must all be > 0")
        if self.bucket not in ("upper", "lower"):
            raise ValueError("bucket must be 'upper' or 'lower'")
        if self.bucket == "lower" and not self.two_store:
            raise ValueError("bucket='lower' requires maxvol_br > 0")
        if self.two_store:
            if self.extravol_br < 0:
                raise ValueError("extravol_br must be >= 0")
            if not 0.0 <= self.vol_br <= self.maxvol_br + self.extravol_br:
                raise ValueError("vol_br must be in [0, maxvol_br + extravol_br]")
            if not 0.1 <= self.gamma_br <= 10.0:
                raise ValueError("gamma_br outside [0.1, 10]; Fortran clamps silently")
            if len(self.epot_br) != n:
                raise ValueError("epot_br must be dense daily when the lower store is on")
            if any(v < 0 for v in self.epot_br):
                raise ValueError("epot_br must be >= 0")
        if not self.outdays or self.outdays != sorted(set(self.outdays)):
            raise ValueError("outdays must be non-empty, strictly increasing")
        if self.outdays[0] < 1:
            raise ValueError("outdays must be >= 1")
        for attr in ("cropfac", "gamma", "rain", "epot", "irrigcode", "gwirrigfrac"):
            if len(getattr(self, attr)) != n:
                raise ValueError(f"{attr} must have exactly numdays={n} entries")
        if any(v < 0 for v in self.cropfac) or any(v < 0 for v in self.gamma):
            raise ValueError("cropfac and gamma must be >= 0")
        if any(v < 0 for v in self.rain):
            raise ValueError("rainfall must be >= 0")
        if any(c not in (0, 1) for c in self.irrigcode):
            raise ValueError("irrigcode must be 0 or 1")
        if any(not 0.0 <= f <= 1.0 for f in self.gwirrigfrac):
            raise ValueError("gwirrigfrac must be in [0, 1]")

    # ------------------------------------------------------------------
    def write(self, directory: Path) -> Path:
        """Write the ``.in`` control file and its forcing files; return the .in."""
        self.validate()
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        days = range(1, self.numdays + 1)

        # Dense daily forcing: one row per day, so no gap-filling rule fires.
        (directory / "veg.dat").write_text(
            "".join(f"{d} {_r(c)} {_r(g)}\n"
                    for d, c, g in zip(days, self.cropfac, self.gamma, strict=True))
        )
        (directory / "rain.dat").write_text(
            "".join(f"{d} {_r(v)}\n" for d, v in zip(days, self.rain, strict=True))
        )
        (directory / "epot.dat").write_text(
            "".join(f"{d} {_r(v)}\n" for d, v in zip(days, self.epot, strict=True))
        )
        (directory / "irrig.dat").write_text(
            "".join(f"{d} {c} {_r(f)}\n"
                    for d, c, f in zip(days, self.irrigcode, self.gwirrigfrac, strict=True))
        )
        if self.two_store:
            (directory / "epot_br.dat").write_text(
                "".join(f"{d} {_r(v)}\n" for d, v in zip(days, self.epot_br, strict=True))
            )

        lines: list[str] = ["* earth properties"]
        if self.two_store:
            lines.append(f"{_r(self.maxvol)} {_r(self.irrigvolfrac)} "
                         f"{_r(self.maxvol_br)} {_r(self.extravol_br)} {_r(self.gamma_br)}")
        else:
            lines.append(f"{_r(self.maxvol)} {_r(self.irrigvolfrac)}")
        lines.append(f"{_r(self.rdelay)} {_r(self.mdelay)}")
        upper_soil = f"{_r(self.ks)} {_r(self.m)} {_r(self.l)} {_r(self.mflowmax)}"
        if self.two_store:
            upper_soil += f" {_r(self.ks_br)} {_r(self.m_br)} {_r(self.l_br)}"
        lines.append(upper_soil)

        lines.append("* volume to elevation")
        v2e = (f"{_r(self.offset)} {_r(self.factor1)} "
               f"{_r(self.factor2)} {_r(self.power)}")
        if self.two_store:
            v2e += f" {self.bucket}"
        if self.elevmin is not None:
            v2e += f" {_r(self.elevmin)}"
            if self.elevmax is not None:
                v2e += f" {_r(self.elevmax)}"
        lines.append(v2e)

        lines.append("* topographic surface")
        lines.append(_r(self.datum))

        lines.append("* initial conditions")
        lines.append(f"{_r(self.vol)} {_r(self.vol_br)}" if self.two_store else _r(self.vol))
        lines.append(f"{len(self.rbuf)} {len(self.mbuf)}")
        lines.append(" ".join(_r(v) for v in self.rbuf))
        lines.append(" ".join(_r(v) for v in self.mbuf))

        lines.append("* solution parameters")
        lines.append(f"{self.nstep} {self.mxiter} {_r(self.tol)}")

        lines.append("* timing information")
        lines.append(f"{self.numdays} {len(self.outdays)}")
        lines.extend(
            " ".join(str(d) for d in self.outdays[i:i + 20])
            for i in range(0, len(self.outdays), 20)
        )

        lines.append("* data filenames")
        lines += ["veg.dat", "rain.dat", "epot.dat", "irrig.dat"]
        if self.two_store:
            lines.append("epot_br.dat")

        path = directory / "case.in"
        path.write_text("\n".join(lines) + "\n")
        return path

    # ------------------------------------------------------------------
    def to_json(self) -> dict:
        return asdict(self)

    @classmethod
    def from_json(cls, data: dict) -> Case:
        return cls(**data)


def dump_cases(cases: list[Case], path: Path) -> None:
    """Freeze case specifications, gzipped.

    Dense daily forcing makes these bulky -- 5.4 MB of JSON for 200 cases, and
    1.8 MB compressed.  They are frozen rather than regenerated from a seed so
    that reproducing the golden files never depends on NumPy's RNG streams
    staying stable across versions.
    """
    blob = json.dumps([c.to_json() for c in cases], indent=1).encode()
    with gzip.GzipFile(path, "wb", compresslevel=9, mtime=0) as fh:
        fh.write(blob)


def load_cases(path: Path) -> list[Case]:
    with gzip.open(path, "rb") as fh:
        return [Case.from_json(d) for d in json.loads(fh.read().decode())]
