"""Generate valid, varied LUMPREM2 cases to drive the reference oracle.

Coverage is deliberate rather than purely random.  A block of structured cases
pins each fidelity trap from the conversion plan -- the store-empties fix-up,
the macropore split, fractional delay buffers, the irrigation convergence
criterion, the clamped volume-to-elevation conversion -- and a randomised
block fills in around them.

Four inputs (MFLOWMAX, NSTEP, MXITER, TOL) are absent from the
``lumprem_variables.rec`` echo and so cannot be checked by reading it back.
They are covered instead by ``SENSITIVITY_PAIRS``: pairs of cases identical
but for one of those values, which a test asserts actually produce different
output.  If the writer ever dropped one on the floor, that test fails.

All forcing is dense daily, so the Fortran's gap-filling rules never fire.
"""

from __future__ import annotations

import numpy as np

from .infile import Case

SEED = 20260912

# Cases that must differ from each other, keyed by the parameter under test.
# These exist because lumprem_variables.rec does not echo these four values.
SENSITIVITY_PAIRS = (
    ("sens_mflowmax_lo", "sens_mflowmax_hi", "mflowmax"),
    ("sens_nstep_lo", "sens_nstep_hi", "nstep"),
    ("sens_tol_loose", "sens_tol_tight", "tol"),
    ("sens_mxiter_lo", "sens_mxiter_hi", "mxiter"),
)


def _seasonal(rng: np.random.Generator, n: int, lo: float, hi: float) -> list[float]:
    """A smooth annual cycle with mild noise, scaled into [lo, hi]."""
    t = np.arange(n)
    phase = rng.uniform(0, 2 * np.pi)
    wave = 0.5 * (1.0 + np.sin(2 * np.pi * t / 365.0 + phase))
    wave = np.clip(wave + rng.normal(0, 0.03, n), 0.0, 1.0)
    return [float(v) for v in lo + (hi - lo) * wave]


def _storms(rng: np.random.Generator, n: int, wet_frac: float, scale: float) -> list[float]:
    """Intermittent rainfall: mostly dry days with exponential wet-day depths."""
    wet = rng.random(n) < wet_frac
    depth = rng.exponential(scale, n) * wet
    return [float(v) for v in depth]


def _blocks(rng: np.random.Generator, n: int, on_frac: float) -> list[int]:
    """Irrigation on/off in runs of days rather than day-by-day noise."""
    out, i = [], 0
    while i < n:
        run = int(rng.integers(5, 40))
        val = 1 if rng.random() < on_frac else 0
        out.extend([val] * min(run, n - i))
        i += run
    return out[:n]


def _outdays(rng: np.random.Generator, n: int) -> list[int]:
    """Output times, always ending exactly on the last simulated day."""
    step = int(rng.integers(5, 40))
    days = list(range(step, n, step)) + [n]
    return sorted(set(days))


def _base(n: int) -> dict:
    """A valid, unremarkable one-store case to vary from."""
    return dict(
        maxvol=0.5, irrigvolfrac=0.3, rdelay=2.0, mdelay=1.0,
        ks=0.9, m=0.2, l=0.55, mflowmax=0.1,
        maxvol_br=0.0, extravol_br=0.0, gamma_br=0.0,
        ks_br=0.0, m_br=0.0, l_br=0.0,
        offset=15.0, factor1=1.0, factor2=1.11, power=2.0,
        bucket="upper", elevmin=None, elevmax=None, datum=120.0,
        vol=0.15, vol_br=0.0, rbuf=[0.0] * 5, mbuf=[0.0] * 5,
        nstep=5, mxiter=500, tol=1e-5,
        numdays=n, outdays=list(range(30, n + 1, 30)) or [n],
        cropfac=[0.7] * n, gamma=[4.5] * n,
        rain=[0.01] * n, epot=[0.004] * n,
        irrigcode=[0] * n, gwirrigfrac=[0.0] * n, epot_br=[],
    )


def structured_cases() -> list[Case]:
    """Hand-built cases, one per fidelity trap and per coverage requirement."""
    rng = np.random.default_rng(SEED)
    out: list[Case] = []

    def add(name: str, n: int = 180, **over):
        spec = _base(n)
        spec.update(over)
        if spec["maxvol_br"] > 0 and not spec["epot_br"]:
            spec["epot_br"] = [0.001] * n
        out.append(Case(name=name, **spec))

    # -- stores driven to their limits (trap 5: the store-empties fix-up)
    # Starts part-full and is drawn down to nothing, which is what makes the
    # store-empties fix-up fire.  Starting empty would leave it a no-op.
    add("dry_no_rain", rain=[0.0] * 180, epot=[0.02] * 180, vol=0.3)
    add("dry_high_epot", rain=[0.0] * 180, epot=[0.5] * 180, cropfac=[1.4] * 180)
    add("start_full", vol=0.5)
    add("start_empty", vol=0.0)

    # -- overflow, macropore split and runoff (trap 6)
    add("deluge", rain=[2.0] * 180, mflowmax=0.1)
    add("no_macropore", rain=[2.0] * 180, mflowmax=0.0)
    add("all_macropore", rain=[2.0] * 180, mflowmax=1.0e6)
    add("sens_mflowmax_lo", rain=[1.0] * 180, mflowmax=0.05)
    add("sens_mflowmax_hi", rain=[1.0] * 180, mflowmax=0.50)

    # -- delay buffers (trap 7)
    add("delay_zero", rdelay=0.0, mdelay=0.0)
    add("delay_fractional", rdelay=2.5, mdelay=0.25)
    add("delay_tiny_fraction", rdelay=0.1, mdelay=0.9)
    add("delay_long", rdelay=30.0, mdelay=17.5)
    add("delay_mismatched", rdelay=0.0, mdelay=25.5)
    # buffers pre-loaded past irdelay, so the first-pass claim sweep fires
    add("buffers_prefilled", rdelay=2.0, mdelay=1.0,
        rbuf=[0.003, 0.002, 0.001, 0.004, 0.005, 0.006],
        mbuf=[0.002, 0.001, 0.003, 0.004])
    add("buffers_long_prefilled", rdelay=6.5, mdelay=4.5,
        rbuf=[0.001 * (i + 1) for i in range(12)],
        mbuf=[0.002 * (i + 1) for i in range(9)])

    # -- irrigation (trap 4: the second convergence criterion)
    add("irrig_always", irrigcode=[1] * 180, irrigvolfrac=0.6, gwirrigfrac=[1.0] * 180)
    add("irrig_never", irrigcode=[0] * 180, irrigvolfrac=0.6)
    add("irrig_blocks", irrigcode=_blocks(rng, 180, 0.5), irrigvolfrac=0.45,
        gwirrigfrac=[0.4] * 180, rain=_storms(rng, 180, 0.25, 0.01))
    add("irrig_zero_frac", irrigcode=[1] * 180, irrigvolfrac=0.0)
    add("irrig_full_frac", irrigcode=[1] * 180, irrigvolfrac=1.0)
    add("irrig_gw_split", irrigcode=[1] * 180, irrigvolfrac=0.5,
        gwirrigfrac=[0.0 if i % 2 else 1.0 for i in range(180)])

    # -- volume to elevation (trap 8: the clamp and the 1e-10 floor)
    add("elev_negative_power", power=-0.5, factor2=1.0e-4,
        elevmin=100.0, elevmax=130.0, vol=0.0)
    add("elev_clamped_both", elevmin=15.2, elevmax=15.4)
    add("elev_fractional_power", power=0.35, factor2=2.0)

    # -- the lower store (traps 9, 10, 12)
    add("two_store_basic", maxvol_br=0.4, extravol_br=0.0, gamma_br=2.0,
        ks_br=0.3, m_br=0.25, l_br=0.6, vol_br=0.1)
    # m is raised and ks_br lowered so the upper store stays solvable while
    # still delivering more water than the lower store can drain -- which is
    # what fills EXTRAVOL_BR and drives the overflow_lower column.  A stiffer
    # upper store simply never converges here, however fine the substep.
    add("two_store_extravol", maxvol_br=0.4, extravol_br=0.3, gamma_br=2.0,
        ks_br=0.02, m_br=0.25, l_br=0.6, vol_br=0.1, m=0.8, rain=[0.2] * 180)
    add("two_store_lower_elev", maxvol_br=0.4, extravol_br=0.1, gamma_br=2.0,
        ks_br=0.3, m_br=0.25, l_br=0.6, vol_br=0.2, bucket="lower")
    add("two_store_start_full", maxvol_br=0.4, extravol_br=0.2, gamma_br=5.0,
        ks_br=0.3, m_br=0.25, l_br=0.6, vol_br=0.6)
    add("two_store_dry", maxvol_br=0.4, extravol_br=0.0, gamma_br=9.9,
        ks_br=0.3, m_br=0.25, l_br=0.6, vol_br=0.0,
        rain=[0.0] * 180, epot_br=[0.05] * 180)
    add("two_store_no_epot_br", maxvol_br=0.4, extravol_br=0.0, gamma_br=0.1,
        ks_br=0.3, m_br=0.25, l_br=0.6, vol_br=0.2, epot_br=[0.0] * 180)
    add("two_store_delayed", maxvol_br=0.5, extravol_br=0.1, gamma_br=3.0,
        ks_br=0.2, m_br=0.3, l_br=0.7, vol_br=0.15, rdelay=4.5, mdelay=2.5)

    # -- two defects in the reference that break its own water balance.  Both
    #    are reproducible, and both are recorded deliberately: Phase 1 has to
    #    decide whether to reproduce them or depart from them, and either way
    #    it needs a reference for what the Fortran actually does.
    #
    # The store-empties fix-up keys on `tempvol == 0`, but irrigation raises
    # tvol to vvol *before* tempvol is taken.  So whenever the negative-volume
    # clamp fires and irrigation then tops the store up, the fix-up is skipped
    # and the full pre-clamp demand is charged against water that was never
    # there.  Here the evaporative demand per substep far exceeds the
    # irrigation target volume, so the clamp fires on every substep.
    add("irrig_bypasses_empty_fixup", n=60, maxvol=1.0, irrigvolfrac=0.02,
        irrigcode=[1] * 60, gwirrigfrac=[0.0] * 60, vol=0.02,
        ks=0.01, m=0.5, l=1.0, rain=[0.0] * 60, epot=[5.0] * 60,
        cropfac=[1.0] * 60, gamma=[4.5] * 60, outdays=[15, 30, 45, 60])

    # `iflag_br` is initialised to 1 and only ever cleared inside the branch
    # that never runs, so `recharge_for_br` is overwritten each day rather than
    # accumulating.  Buffer water claimed by the first-pass sweep -- anything
    # sitting beyond irdelay when a call starts -- is therefore added to
    # `recharge` but never routed into the lower store, and simply vanishes.
    add("two_store_buffer_claim_lost", n=60, maxvol_br=0.4, extravol_br=0.1,
        gamma_br=2.0, ks_br=0.3, m_br=0.5, l_br=0.6, vol_br=0.1, m=0.8,
        rdelay=2.0, rbuf=[0.01] * 8, outdays=[15, 30, 45, 60])

    # -- solver settings, and the call partitioning that drives rechmod
    add("sens_nstep_lo", nstep=1)
    add("sens_nstep_hi", nstep=20)
    add("sens_tol_loose", tol=1.0e-3)
    add("sens_tol_tight", tol=1.0e-9)
    add("sens_mxiter_lo", mxiter=2, tol=1.0e-9, allow_nonconvergence=True)
    add("sens_mxiter_hi", mxiter=800, tol=1.0e-9)
    add("outdays_single", outdays=[180])
    add("outdays_daily", n=60, outdays=list(range(1, 61)))
    add("outdays_irregular", outdays=[1, 2, 3, 17, 90, 91, 179, 180])

    return out


def random_cases(count: int, seed: int = SEED) -> list[Case]:
    """Randomised but valid parameter sets across both configurations."""
    rng = np.random.default_rng(seed + 1)
    out: list[Case] = []
    for i in range(count):
        n = int(rng.integers(90, 400))
        two_store = bool(rng.random() < 0.45)
        maxvol = float(rng.uniform(0.05, 20.0))
        # gamma is kept away from zero: evap() divides by
        # (1 - 2*exp(-gamma) + exp(-gamma*vd)), which is 0/0 at gamma == 0.
        # The reader clamps gamma_br into [0.1, 10] but leaves upper gamma free.
        gamma_lo = float(rng.uniform(0.3, 3.0))
        gamma_hi = gamma_lo + float(rng.uniform(0.0, 6.0))
        irrig_on = bool(rng.random() < 0.5)
        rain_mode = rng.choice(["dry", "normal", "wet", "extreme"], p=[0.1, 0.5, 0.3, 0.1])
        if rain_mode == "dry":
            rain = [0.0] * n
        elif rain_mode == "normal":
            rain = _storms(rng, n, 0.25, maxvol * 0.05)
        elif rain_mode == "wet":
            rain = _storms(rng, n, 0.6, maxvol * 0.2)
        else:
            rain = _storms(rng, n, 0.9, maxvol * 2.0)

        spec = dict(
            maxvol=maxvol,
            irrigvolfrac=float(rng.uniform(0.0, 1.0)) if irrig_on else 0.0,
            rdelay=float(rng.choice([0.0, rng.uniform(0.0, 40.0)])),
            mdelay=float(rng.choice([0.0, rng.uniform(0.0, 20.0)])),
            ks=float(rng.uniform(0.001, 15.0)),
            m=float(rng.uniform(0.05, 1.0)),
            l=float(rng.uniform(0.1, 3.0)),
            mflowmax=float(rng.choice([0.0, rng.uniform(0.0, 2.0)])),
            offset=float(rng.uniform(-20.0, 40.0)),
            factor1=float(rng.uniform(0.0, 3.0)),
            factor2=float(rng.uniform(0.0, 3.0)),
            power=float(rng.choice([rng.uniform(0.1, 3.0), rng.uniform(-1.0, -0.1)])),
            datum=float(rng.uniform(40.0, 250.0)),
            vol=float(rng.choice([0.0, maxvol, rng.uniform(0.0, maxvol)])),
            rbuf=[float(v) for v in rng.uniform(0.0, 0.01, int(rng.integers(1, 15)))],
            mbuf=[float(v) for v in rng.uniform(0.0, 0.01, int(rng.integers(1, 15)))],
            nstep=int(rng.integers(1, 13)),
            mxiter=int(rng.integers(200, 1000)),
            tol=float(10.0 ** rng.uniform(-9, -3)),
            numdays=n,
            outdays=_outdays(rng, n),
            cropfac=_seasonal(rng, n, 0.1, float(rng.uniform(0.8, 1.5))),
            gamma=_seasonal(rng, n, gamma_lo, gamma_hi),
            rain=rain,
            epot=_seasonal(rng, n, 0.0, float(rng.uniform(0.002, 0.05))),
            irrigcode=_blocks(rng, n, 0.5) if irrig_on else [0] * n,
            gwirrigfrac=_seasonal(rng, n, 0.0, 1.0),
            bucket="upper",
            elevmin=None,
            elevmax=None,
            epot_br=[],
        )
        # A negative power sends factor2 * dtemp**power towards 1e20 at the
        # 1e-10 volume floor, so those cases must carry an elevation clamp --
        # which is exactly why the clamp is there.
        if spec["power"] < 0:
            spec["elevmin"] = float(spec["offset"] - rng.uniform(1.0, 20.0))
            spec["elevmax"] = float(spec["offset"] + rng.uniform(1.0, 50.0))

        if two_store:
            maxvol_br = float(rng.uniform(0.05, 10.0))
            extravol_br = float(rng.choice([0.0, rng.uniform(0.0, 5.0)]))
            spec.update(
                maxvol_br=maxvol_br,
                extravol_br=extravol_br,
                gamma_br=float(rng.uniform(0.1, 10.0)),
                ks_br=float(rng.uniform(0.001, 5.0)),
                m_br=float(rng.uniform(0.05, 1.0)),
                l_br=float(rng.uniform(0.1, 3.0)),
                vol_br=float(rng.choice([0.0, rng.uniform(0.0, maxvol_br + extravol_br)])),
                epot_br=_seasonal(rng, n, 0.0, float(rng.uniform(0.0, 0.02))),
                bucket=str(rng.choice(["upper", "lower"])),
            )
        else:
            spec.update(maxvol_br=0.0, extravol_br=0.0, gamma_br=0.0,
                        ks_br=0.0, m_br=0.0, l_br=0.0, vol_br=0.0)

        out.append(Case(name=f"random_{i:04d}", **spec))
    return out


def generate_cases(n_random: int = 160, seed: int = SEED) -> list[Case]:
    cases = structured_cases() + random_cases(n_random, seed)
    names = [c.name for c in cases]
    assert len(set(names)) == len(names), "duplicate case names"
    for c in cases:
        c.validate()
    return cases
