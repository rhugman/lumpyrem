"""The three gap-filling rules, tested as rules rather than as file formats.

``lumprem2.f`` fills gaps in its forcing differently depending on which file a
value was written in: vegetation is linearly interpolated, potential
evaporation and irrigation are forward-filled as steps, and rainfall is
zero-filled.  Phase 0 deliberately generated dense daily forcing so that none
of this fires during the fidelity comparison, which leaves it untested there --
so it is tested here, twice over: against hand-built expectations, and against
the reference itself driven by genuinely sparse forcing files.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from api_compare import EPOCH, model_for
from lumpyrem import Fill, Forcing, ForcingError
from lumpyrem.forcing import FORTRAN_FILL, VARIABLES
from oracle.results import read_csv_output
from oracle.runner import run_case

DAY = pd.Timedelta(days=1)


def _frame(days, **columns) -> pd.DataFrame:
    """A forcing frame whose row ``k`` is simulation day ``days[k]``."""
    index = pd.DatetimeIndex([EPOCH + (d - 1) * DAY for d in days])
    return pd.DataFrame(columns, index=index)


#: Stand-ins for whatever a test is not exercising, so every run is well posed.
_FILLERS = {"rainfall": 0.0, "pot_evap": 0.0, "veg_gamma": 1.0}


def _daily(days, ndays, *, fill=None, **columns):
    frame = _frame(days, **columns)
    spare = {k: v for k, v in _FILLERS.items() if k not in columns}
    return Forcing.from_dataframe(frame, fill=fill, **spare) \
        .to_daily(EPOCH, EPOCH + (ndays - 1) * DAY)


# ---------------------------------------------------------------------------
# The rules, one at a time
# ---------------------------------------------------------------------------

def test_the_defaults_are_the_fortran_rules():
    assert FORTRAN_FILL == {
        "rainfall": Fill.ZERO,
        "pot_evap": Fill.FORWARD,
        "crop_factor": Fill.INTERPOLATE,
        "veg_gamma": Fill.INTERPOLATE,
        "irrigate": Fill.FORWARD,
        "gw_irrig_frac": Fill.FORWARD,
        "pot_evap_lower": Fill.FORWARD,
    }
    assert set(FORTRAN_FILL) == set(VARIABLES)


def test_rainfall_is_zero_filled():
    out = _daily([3, 5, 6], 8, rainfall=[4.0, 1.0, 2.0])["rainfall"]
    assert out.tolist() == [0.0, 0.0, 4.0, 0.0, 1.0, 2.0, 0.0, 0.0]


def test_evaporation_is_a_step_held_forward():
    out = _daily([1, 4, 6], 8, pot_evap=[2.0, 5.0, 1.0])["pot_evap"]
    assert out.tolist() == [2.0, 2.0, 2.0, 5.0, 5.0, 1.0, 1.0, 1.0]


def test_irrigation_is_forward_filled_as_a_step_not_interpolated():
    """An interpolated irrigation code would be a fraction, which means nothing."""
    out = _daily([1, 5], 6, irrigate=[0.0, 1.0], gw_irrig_frac=[0.0, 0.75])
    assert out["irrigate"].tolist() == [0.0, 0.0, 0.0, 0.0, 1.0, 1.0]
    assert out["gw_irrig_frac"].tolist() == [0.0, 0.0, 0.0, 0.0, 0.75, 0.75]


def test_vegetation_is_linearly_interpolated():
    out = _daily([1, 11], 11, crop_factor=[0.0, 10.0])["crop_factor"]
    assert out[0] == 0.0
    assert out[-1] == 10.0
    assert out.tolist() == pytest.approx(list(range(11)), abs=1e-12)


def test_interpolation_keeps_the_fortran_reciprocal_form():
    """`v1 + dv * ((day - d1) * (1/(d2 - d1)))`, not a division per day.

    The two differ in the last bit, and ``lumprem2.f`` computes the reciprocal
    once per segment.  Reproducing that is what lets a run driven by sparse
    forcing still match the reference exactly, so it is pinned here rather than
    left to look like an implementation detail.
    """
    out = _daily([1, 11], 11, crop_factor=[0.0, 10.0])["crop_factor"]
    naive = np.interp(np.arange(1.0, 12.0), [1.0, 11.0], [0.0, 10.0])
    assert out[3] == 0.0 + 10.0 * (3.0 * (1.0 / 10.0))
    assert out[3] != naive[3]
    assert (out != naive).any()


def test_a_single_listed_vegetation_value_is_held_flat():
    out = _daily([1], 5, crop_factor=[0.7])["crop_factor"]
    assert out.tolist() == [0.7] * 5


def test_dense_requires_a_value_on_every_day():
    # The same series is fine under the default forward rule ...
    assert _daily([1, 3], 4, pot_evap=[1.0, 2.0])["pot_evap"].tolist() == \
        [1.0, 1.0, 2.0, 2.0]
    # ... and refused outright under 'dense', which never guesses.
    with pytest.raises(ForcingError, match="requires a value on every day, and"):
        _daily([1, 3], 4, fill={"pot_evap": "dense"}, pot_evap=[1.0, 2.0])


def test_rules_can_be_swapped_per_variable():
    """The whole point of naming them: rainfall can be held forward instead."""
    zeroed = _daily([1, 4], 5, rainfall=[2.0, 6.0])["rainfall"]
    assert zeroed.tolist() == [2.0, 0.0, 0.0, 6.0, 0.0]
    stepped = _daily([1, 4], 5, fill={"rainfall": Fill.FORWARD},
                     rainfall=[2.0, 6.0])["rainfall"]
    assert stepped.tolist() == [2.0, 2.0, 2.0, 6.0, 6.0]


def test_nan_is_a_gap_not_a_value():
    out = _daily([1, 2, 3], 3, pot_evap=[1.0, np.nan, 4.0])["pot_evap"]
    assert out.tolist() == [1.0, 1.0, 4.0]


# ---------------------------------------------------------------------------
# What the rules refuse to guess
# ---------------------------------------------------------------------------

def test_forward_fill_refuses_to_invent_a_leading_value():
    with pytest.raises(ForcingError, match="nothing to carry forward"):
        _daily([3, 5], 6, pot_evap=[1.0, 2.0])


def test_interpolation_refuses_to_extrapolate():
    with pytest.raises(ForcingError, match="on or outside both ends"):
        _daily([1, 4], 8, crop_factor=[1.0, 2.0])
    with pytest.raises(ForcingError, match="on or outside both ends"):
        _daily([2, 8], 8, crop_factor=[1.0, 2.0])


def test_a_value_listed_before_the_run_still_anchors_the_fill():
    """The Fortran cannot express this -- its files are 1-based day numbers."""
    out = _daily([-4, 6], 7, pot_evap=[3.0, 7.0])["pot_evap"]
    assert out.tolist() == [3.0] * 5 + [7.0, 7.0]


def test_bad_forcing_names_the_day_and_the_variable():
    with pytest.raises(ForcingError, match="rainfall must be >= 0; day 2"):
        _daily([1, 2], 2, rainfall=[0.0, -1.0])
    with pytest.raises(ForcingError, match="veg_gamma must be > 0"):
        _daily([1, 2], 2, veg_gamma=[1.0, 0.0])
    with pytest.raises(ForcingError, match=r"gw_irrig_frac must be in \[0, 1\]; day 1"):
        _daily([1], 2, gw_irrig_frac=[1.5])


def test_an_interpolated_irrigation_code_is_caught_not_rounded():
    """Half-irrigating a store is not a thing; the rule swap has to be refused."""
    with pytest.raises(ForcingError, match="irrigate must be 0 or 1"):
        _daily([1, 5], 5, fill={"irrigate": "interpolate"}, irrigate=[0.0, 1.0])


def test_a_misspelled_column_is_rejected_not_ignored():
    with pytest.raises(ForcingError, match="unknown forcing column"):
        Forcing.from_dataframe(_frame([1], rain=[1.0]))
    with pytest.raises(ForcingError, match="not a forcing variable"):
        Forcing.from_dataframe(_frame([1], rainfall=[1.0]), fill={"rain": "zero"})


def test_sub_daily_stamps_are_refused_with_a_reason():
    frame = pd.DataFrame({"rainfall": [1.0]},
                         index=pd.DatetimeIndex(["2001-01-03 06:00"]))
    with pytest.raises(ForcingError, match="midnight"):
        Forcing.from_dataframe(frame)


# ---------------------------------------------------------------------------
# Constructors
# ---------------------------------------------------------------------------

def test_every_constructor_reaches_the_same_daily_arrays(tmp_path):
    days = pd.date_range("2001-01-03", periods=4, freq="D")
    rain = [0.0, 3.0, 0.0, 1.5]
    epot = [2.0, 2.0, 2.5, 2.5]
    frame = pd.DataFrame({"rainfall": rain, "pot_evap": epot}, index=days)

    from_frame = Forcing.from_dataframe(frame, veg_gamma=4.0)
    from_arrays = Forcing.from_arrays("2001-01-03", rainfall=rain, pot_evap=epot,
                                      veg_gamma=4.0)
    from_dict = Forcing.from_dict({"time": days, "rainfall": rain,
                                   "pot_evap": epot, "veg_gamma": 4.0})
    path = tmp_path / "climate.csv"
    frame.assign(veg_gamma=4.0).rename_axis("time").to_csv(path)
    from_csv = Forcing.from_csv(path)

    want = from_frame.to_daily()
    for other in (from_arrays, from_dict, from_csv):
        for name in VARIABLES:
            assert other.to_daily()[name].tolist() == want[name].tolist(), name
    assert from_csv.span == (days[0], days[-1])


def test_coerce_accepts_what_run_accepts(tmp_path):
    days = pd.date_range("2001-01-03", periods=3, freq="D")
    frame = pd.DataFrame({"rainfall": [1.0, 0.0, 2.0], "pot_evap": [1.0] * 3}, index=days)
    direct = Forcing.from_dataframe(frame, veg_gamma=3.0)
    assert Forcing.coerce(direct) is direct
    assert Forcing.coerce(frame, veg_gamma=3.0).to_daily()["rainfall"].tolist() == \
        direct.to_daily()["rainfall"].tolist()
    with pytest.raises(ForcingError, match="cannot read forcing from"):
        Forcing.coerce(42)


def test_a_missing_required_variable_is_named():
    frame = _frame([1, 2], rainfall=[1.0, 0.0])
    with pytest.raises(ForcingError, match=r"missing \['pot_evap', 'veg_gamma'\]"):
        Forcing.from_dataframe(frame).to_daily()


# ---------------------------------------------------------------------------
# Against the reference
# ---------------------------------------------------------------------------

def _sparse_days(n: int, stride: int, first: int = 1) -> list[int]:
    return list(range(first, n + 1, stride))


@pytest.mark.oracle
def test_fill_rules_match_the_reference_on_sparse_forcing(oracle_exe, frozen_cases, rundir,
                                                          agreement):
    """Drive the Fortran from genuinely sparse files and reproduce it exactly.

    This is the test dense forcing was designed to make impossible during the
    fidelity gates, and the only one that checks the fill rules against the
    implementation they were copied from rather than against a reading of it.

    Exactly means 0 ULP where Python and the Fortran share a maths library.
    Where they do not, the tolerance regime still catches any fill rule that
    lands on the wrong day or the wrong value, but not trap 16's last-bit
    reciprocal -- that is left to the platforms that can see it.
    """
    case = next(c for c in frozen_cases
                if not c.two_store and not c.allow_nonconvergence
                and c.numdays >= 90 and c.irrigvolfrac > 0.0)
    n = case.numdays

    # Chosen so each rule is exercised on a different stride, and so rainfall
    # starts late -- the one variable the Fortran lets begin after day 1.
    veg_days = _sparse_days(n, 17) + [n + 3]
    rain_days = _sparse_days(n, 3, first=4)
    epot_days = _sparse_days(n, 11)
    irrig_days = _sparse_days(n, 7)

    def rewrite(workdir):
        (workdir / "veg.dat").write_text("".join(
            f"{d} {case.cropfac[min(d, n) - 1]!r} {case.gamma[min(d, n) - 1]!r}\n"
            for d in veg_days))
        (workdir / "rain.dat").write_text("".join(
            f"{d} {case.rain[d - 1]!r}\n" for d in rain_days))
        (workdir / "epot.dat").write_text("".join(
            f"{d} {case.epot[d - 1]!r}\n" for d in epot_days))
        (workdir / "irrig.dat").write_text("".join(
            f"{d} {case.irrigcode[d - 1]} {case.gwirrigfrac[d - 1]!r}\n"
            for d in irrig_days))

    run = run_case(oracle_exe, case, rundir / "sparse", rewrite=rewrite)
    want = read_csv_output(rundir / "sparse" / "case.csv").values

    listed = sorted({*veg_days, *rain_days, *epot_days, *irrig_days})
    frame = pd.DataFrame(
        {name: [np.nan] * len(listed) for name in
         ("rainfall", "pot_evap", "crop_factor", "veg_gamma", "irrigate", "gw_irrig_frac")},
        index=pd.DatetimeIndex([EPOCH + (d - 1) * DAY for d in listed]),
    )
    at = {d: i for i, d in enumerate(listed)}
    for d in veg_days:
        frame.iloc[at[d], frame.columns.get_loc("crop_factor")] = case.cropfac[min(d, n) - 1]
        frame.iloc[at[d], frame.columns.get_loc("veg_gamma")] = case.gamma[min(d, n) - 1]
    for d in rain_days:
        frame.iloc[at[d], frame.columns.get_loc("rainfall")] = case.rain[d - 1]
    for d in epot_days:
        frame.iloc[at[d], frame.columns.get_loc("pot_evap")] = case.epot[d - 1]
    for d in irrig_days:
        frame.iloc[at[d], frame.columns.get_loc("irrigate")] = case.irrigcode[d - 1]
        frame.iloc[at[d], frame.columns.get_loc("gw_irrig_frac")] = case.gwirrigfrac[d - 1]

    got = model_for(case).run(Forcing.from_dataframe(frame),
                              start=EPOCH, end=EPOCH + (n - 1) * DAY,
                              times=case.outdays).values

    assert got.shape == want.shape
    worst, _ = agreement.worst(got, want)
    assert worst <= agreement.budget(ulp_budget=0), (
        f"{case.name} under sparse forcing: worst {worst:.3g} {agreement.unit()}; "
        f"{agreement.describe()}. The Python fill rules and the Fortran's have "
        f"diverged."
    )
    assert "itn limit exceeded" not in run.stdout


# ---------------------------------------------------------------------------
# Per-cell forcing
# ---------------------------------------------------------------------------

def _sparse_block(rng, ndays, ncell, *, same_days: bool):
    """Values on scattered days, NaN elsewhere -- listed on the same days in
    every cell, or on each cell's own days."""
    block = rng.uniform(0.5, 3.0, (ndays, ncell))
    if same_days:
        keep = rng.random(ndays) < 0.2
        keep[[0, -1]] = True
        block[~keep] = np.nan
    else:
        keep = rng.random((ndays, ncell)) < 0.2
        keep[[0, -1]] = True
        block[~keep] = np.nan
    return block


@pytest.mark.parametrize("same_days", [True, False], ids=["shared-days", "own-days"])
def test_grid_forcing_fills_each_cell_as_a_forcing_would_alone(same_days):
    """Each column, every rule: exactly what a plain Forcing makes of it."""
    from lumpyrem import GridForcing
    rng = np.random.default_rng(11)
    ndays, ncell = 90, 5
    days = pd.date_range(EPOCH, periods=ndays)
    blocks = {name: _sparse_block(rng, ndays, ncell, same_days=same_days)
              for name in ("rainfall", "pot_evap", "crop_factor", "veg_gamma")}
    shared = rng.uniform(0.0, 1.0, ndays)
    for rules in (None, {"rainfall": "forward", "pot_evap": "interpolate",
                         "crop_factor": "zero", "veg_gamma": "forward"}):
        grid = GridForcing.from_arrays(days, fill=rules, gw_irrig_frac=shared, **blocks)
        daily = grid.to_daily()
        assert daily["rainfall"].shape == (ndays, ncell)
        assert daily["gw_irrig_frac"].shape == (ndays,)
        for c in range(ncell):
            alone = Forcing.from_dataframe(
                pd.DataFrame({k: v[:, c] for k, v in blocks.items()} | {
                    "gw_irrig_frac": shared}, index=days), fill=rules).to_daily()
            for name in VARIABLES:
                col = daily[name] if daily[name].ndim == 1 else daily[name][:, c]
                assert np.array_equal(col, alone[name]), (rules, c, name)


def test_grid_forcing_names_the_cell_that_broke_a_rule():
    from lumpyrem import GridForcing
    rain = np.ones((10, 3))
    rain[4, 2] = -1.0
    with pytest.raises(ForcingError, match="rainfall must be >= 0; day 5, cell 2"):
        GridForcing.from_arrays(EPOCH, rainfall=rain, pot_evap=1.0, veg_gamma=2.0).to_daily()
    with pytest.raises(ForcingError, match="at least one"):
        GridForcing.from_arrays(EPOCH, rainfall=np.ones(10), pot_evap=1.0, veg_gamma=2.0)
    with pytest.raises(ForcingError, match="disagree on the number of cells"):
        GridForcing.from_arrays(EPOCH, rainfall=np.ones((10, 3)), pot_evap=np.ones((10, 2)),
                                veg_gamma=2.0)
    with pytest.raises(ForcingError, match="unknown forcing variable"):
        GridForcing.from_arrays(EPOCH, rain=np.ones((10, 3)))
