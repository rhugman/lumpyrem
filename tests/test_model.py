"""Output scheduling, non-convergence reporting, and the Results object."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

from api_compare import EPOCH, forcing_for, model_for
from lumpyrem import (
    ConvergenceWarning,
    Forcing,
    Model,
    ModelError,
    Solver,
    UpperStore,
    VolumeToElevation,
)
from lumpyrem.results import COLUMNS

DAY = pd.Timedelta(days=1)


@pytest.fixture
def simple():
    """A small, well-behaved run: 90 days from 2001-01-03, mild rainfall."""
    rng = np.random.default_rng(7)
    n = 90
    forcing = Forcing.from_arrays(
        EPOCH,
        rainfall=np.where(rng.random(n) < 0.25, rng.exponential(6.0, n), 0.0),
        pot_evap=np.full(n, 3.0),
        veg_gamma=5.0, crop_factor=0.9,
    )
    model = Model(
        upper=UpperStore(maxvol=40.0, ks=3.0, m=0.4, l=0.5, mflowmax=1.5, rdelay=2.0),
        elevation=VolumeToElevation(offset=5.0, factor1=0.2, factor2=0.1,
                                    power=0.5, datum=60.0),
        solver=Solver(nstep=4),
    )
    return model, forcing


# ---------------------------------------------------------------------------
# Output scheduling
# ---------------------------------------------------------------------------

def test_the_default_schedule_is_every_day(simple):
    model, forcing = simple
    res = model.run(forcing)
    assert res.days.tolist() == list(range(0, 91))
    assert res.times[0] == EPOCH
    assert res.times[-1] == EPOCH + 90 * DAY


def test_a_frequency_string_reports_on_that_schedule(simple):
    model, forcing = simple
    res = model.run(forcing, times="MS")
    assert list(res.times[1:]) == [pd.Timestamp("2001-02-01"),
                                   pd.Timestamp("2001-03-01"),
                                   pd.Timestamp("2001-04-01")]
    # The run begins mid-month, so the first interval is a partial month.
    assert res.days.tolist() == [0, 29, 57, 88]


def test_a_timestamp_marks_the_end_of_the_interval_it_labels(simple):
    """Row k covers the days after row k-1 up to and including row k.

    This is the stress-period convention: ``df.loc[a:b]`` reads as "everything
    that happened after a, up to b", and the totals line up with a MODFLOW
    period boundary without an off-by-one.
    """
    model, forcing = simple
    daily = model.run(forcing)
    monthly = model.run(forcing, times="MS")
    rain = daily["rainfall"]
    first = monthly["rainfall"].iloc[1]
    # 2001-02-01 labels the total over 2001-01-03 .. 2001-01-31 inclusive.
    assert first == pytest.approx(rain.loc["2001-01-04":"2001-02-01"].sum(), rel=1e-12)


def test_explicit_dates_and_day_numbers_agree(simple):
    model, forcing = simple
    by_day = model.run(forcing, times=[10, 40, 90])
    by_date = model.run(forcing, times=[EPOCH + 10 * DAY, EPOCH + 40 * DAY,
                                        EPOCH + 90 * DAY])
    assert np.array_equal(by_day.values, by_date.values)
    assert by_day.days.tolist() == [0, 10, 40, 90]


def test_output_times_outside_the_run_are_refused_by_name(simple):
    model, forcing = simple
    with pytest.raises(ModelError, match="falls outside the run"):
        model.run(forcing, times=[10, 200])
    with pytest.raises(ModelError, match="falls outside the run"):
        model.run(forcing, times=[0, 10])
    with pytest.raises(ModelError, match="falls outside the run"):
        model.run(forcing, times=["2002-01-01"])


def test_output_times_must_be_on_the_day_grid(simple):
    model, forcing = simple
    with pytest.raises(ModelError, match="not midnight"):
        model.run(forcing, times=["2001-01-10 12:00"])


def test_output_times_must_increase(simple):
    model, forcing = simple
    with pytest.raises(ModelError, match="strictly increasing"):
        model.run(forcing, times=[40, 10])
    with pytest.raises(ModelError, match="strictly increasing"):
        model.run(forcing, times=[10, 10])


def test_an_empty_or_unusable_schedule_says_which(simple):
    model, forcing = simple
    with pytest.raises(ModelError, match="times is empty"):
        model.run(forcing, times=[])
    with pytest.raises(ModelError, match="no output time inside the run"):
        model.run(forcing, times="YS")
    with pytest.raises(ModelError, match="sequence of dates or day numbers"):
        model.run(forcing, times=3.5)


def test_the_run_can_be_narrowed_without_touching_the_forcing(simple):
    model, forcing = simple
    res = model.run(forcing, start="2001-01-10", end="2001-02-10", times="7D")
    assert res.times[0] == pd.Timestamp("2001-01-10")
    assert res.times[-1] <= pd.Timestamp("2001-02-10")
    with pytest.raises(ModelError, match="ends .* before it starts"):
        model.run(forcing, start="2001-02-10", end="2001-01-10")


# ---------------------------------------------------------------------------
# Non-convergence
# ---------------------------------------------------------------------------

def test_non_convergence_warns_by_default(frozen_cases):
    """The Fortran prints to stdout and carries on; silence is the bug."""
    case = next(c for c in frozen_cases if c.allow_nonconvergence)
    model = model_for(case, on_nonconvergence="warn")
    with pytest.warns(ConvergenceWarning, match="hit mxiter"):
        res = model.run(forcing_for(case), times=case.outdays)
    assert not res.converged
    assert res.nonconverged_steps > 0


def test_non_convergence_can_be_made_fatal(frozen_cases):
    case = next(c for c in frozen_cases if c.allow_nonconvergence)
    model = model_for(case, on_nonconvergence="raise")
    with pytest.raises(ModelError, match="did not converge"):
        model.run(forcing_for(case), times=case.outdays)


def test_non_convergence_can_be_silenced_deliberately(frozen_cases):
    case = next(c for c in frozen_cases if c.allow_nonconvergence)
    model = model_for(case, on_nonconvergence="ignore")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        res = model.run(forcing_for(case), times=case.outdays)
    assert not res.converged


def test_a_converged_run_is_silent(simple):
    model, forcing = simple
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        res = model.run(forcing)
    assert res.converged and res.nonconverged_steps == 0


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------

def test_the_dataframe_is_indexed_by_time_and_keeps_the_day_number(simple):
    model, forcing = simple
    frame = model.run(forcing, times="MS").to_dataframe()
    assert frame.index.name == "time"
    assert list(frame.columns) == ["day", *COLUMNS]
    assert "depth_to_water" in frame.columns and "depth-to-water" not in frame.columns
    assert frame["day"].tolist()[0] == 0


def test_a_column_can_be_had_as_a_series_or_an_array(simple):
    model, forcing = simple
    res = model.run(forcing, times="MS")
    series = res["total_rech"]
    assert isinstance(series, pd.Series)
    assert series.index.equals(res.times)
    assert np.array_equal(series.to_numpy(), res.column("total_rech"))


def test_an_unknown_column_lists_the_real_ones(simple):
    model, forcing = simple
    res = model.run(forcing)
    with pytest.raises(KeyError, match="total_rech"):
        res.column("recharge")


def test_row_zero_is_the_initial_state(simple):
    model, forcing = simple
    res = model.run(forcing, times="MS")
    assert res.values[0, COLUMNS.index("vol_upper")] == model.initial.vol
    fluxes = [c for c in COLUMNS if c not in
              ("vol_upper", "vol_lower", "vol_drain", "vol_macro",
               "elevation", "depth_to_water")]
    assert all(res.column(c)[0] == 0.0 for c in fluxes)


def test_the_balance_closes_on_a_well_posed_run(simple):
    model, forcing = simple
    res = model.run(forcing, times="MS")
    assert res.balance_error() < 1e-12
    assert res.balance_error(relative=True) < 1e-15


def test_depth_to_water_is_the_datum_less_the_elevation(simple):
    model, forcing = simple
    res = model.run(forcing, times="MS")
    assert np.array_equal(res.column("depth_to_water"),
                          model.elevation.datum - res.column("elevation"))


def test_to_csv_round_trips(simple, tmp_path):
    model, forcing = simple
    res = model.run(forcing, times="MS")
    path = tmp_path / "out.csv"
    res.to_csv(path)
    back = pd.read_csv(path, index_col="time", parse_dates=True)
    assert list(back.columns) == ["day", *COLUMNS]
    assert back["total_rech"].to_numpy() == pytest.approx(res.column("total_rech"))


def test_plot_returns_axes_and_refuses_unknown_columns(simple):
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    model, forcing = simple
    res = model.run(forcing, times="MS")
    ax = res.plot()
    assert ax.get_ylabel() == "depth per interval"
    assert res.plot(["vol_upper"], ax=ax) is ax
    with pytest.raises(KeyError):
        res.plot(["nonsense"])


def test_repr_says_whether_it_converged(simple):
    model, forcing = simple
    text = repr(model.run(forcing, times="MS"))
    assert "3 output times" in text and "converged" in text
