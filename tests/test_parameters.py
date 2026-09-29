"""Validation: what the objects refuse, and how they say so.

The Fortran reader clamps several parameters into range without saying so --
``gamma_br`` into [0.1, 10], both delays at ``MAXDELAY - 2``.  That is worse
than an error in a calibration context, because the knob keeps turning and the
model stops responding.  These objects refuse instead, and every refusal names
the field and the value so a caller can act on it without parsing prose.
"""

from __future__ import annotations

import pytest

from lumpyrem import (
    InitialState,
    LowerStore,
    Model,
    ParameterError,
    Solver,
    UpperStore,
    VolumeToElevation,
)
from lumpyrem.parameters import LOWER_STORE_CROP_FACTOR

GOOD_UPPER = dict(maxvol=0.3, ks=0.05, m=0.4, l=0.5)


def test_a_refusal_carries_the_field_and_the_value():
    with pytest.raises(ParameterError) as exc:
        UpperStore(maxvol=-1.0, ks=0.05, m=0.4, l=0.5)
    err = exc.value
    assert err.owner == "UpperStore"
    assert err.field == "maxvol"
    assert err.value == -1.0
    assert "be > 0" in err.requirement
    assert "UpperStore.maxvol" in str(err) and "-1.0" in str(err)
    assert isinstance(err, ValueError)


@pytest.mark.parametrize("bad,field", [
    (dict(maxvol=0.0), "maxvol"),
    (dict(ks=-1.0), "ks"),
    (dict(m=0.0), "m"),
    (dict(mflowmax=-0.1), "mflowmax"),
    (dict(irrigvolfrac=1.5), "irrigvolfrac"),
    (dict(rdelay=-1.0), "rdelay"),
    (dict(mdelay=float("nan")), "mdelay"),
    (dict(l=float("inf")), "l"),
])
def test_upper_store_rejects(bad, field):
    with pytest.raises(ParameterError) as exc:
        UpperStore(**{**GOOD_UPPER, **bad})
    assert exc.value.field == field


def test_there_is_no_delay_ceiling():
    """lumprem2.f silently sets any delay over MAXDELAY - 2 = 498 to 498.  The
    buffers are sized to the delay now, so there is no ceiling to clamp to or
    refuse at; tests/test_delay_buffers.py runs such a model."""
    assert UpperStore(**GOOD_UPPER, rdelay=5000.5, mdelay=499.0).rdelay == 5000.5
    InitialState(vol=0.0, drain_buffer=(0.1,) * 2000)
    with pytest.raises(ParameterError) as exc:
        UpperStore(**GOOD_UPPER, mdelay=-0.5)
    assert exc.value.field == "mdelay"


def test_the_lower_store_gamma_range_is_refused_rather_than_clamped():
    LowerStore(maxvol=0.5, gamma=0.1, ks=0.01, m=0.5, l=0.5)
    LowerStore(maxvol=0.5, gamma=10.0, ks=0.01, m=0.5, l=0.5)
    for gamma in (0.09, 10.1):
        with pytest.raises(ParameterError, match="clamps silently") as exc:
            LowerStore(maxvol=0.5, gamma=gamma, ks=0.01, m=0.5, l=0.5)
        assert exc.value.field == "gamma"


def test_a_disabled_lower_store_is_none_not_a_zero_object():
    with pytest.raises(ParameterError, match="pass lower=None instead"):
        LowerStore(maxvol=0.0, gamma=1.0, ks=0.01, m=0.5, l=0.5)


def test_lower_store_capacity_includes_the_extra_volume():
    lower = LowerStore(maxvol=0.5, gamma=1.0, ks=0.01, m=0.5, l=0.5, extravol=0.25)
    assert lower.capacity == 0.75


def test_the_lower_store_crop_factor_is_a_constant_not_a_parameter():
    """rechmod2.f marks it `! hardwired`; exposing it stays a decision."""
    assert LOWER_STORE_CROP_FACTOR == 1.0
    with pytest.raises(TypeError):
        LowerStore(maxvol=0.5, gamma=1.0, ks=0.01, m=0.5, l=0.5, crop_factor=0.8)


@pytest.mark.parametrize("bad,field", [
    (dict(nstep=0), "nstep"),
    (dict(mxiter=0), "mxiter"),
    (dict(tol=0.0), "tol"),
    (dict(nstep=2.5), "nstep"),
    (dict(on_nonconvergence="explode"), "on_nonconvergence"),
])
def test_solver_rejects(bad, field):
    with pytest.raises(ParameterError) as exc:
        Solver(**bad)
    assert exc.value.field == field


def test_elevation_clips_must_be_ordered():
    VolumeToElevation(offset=0.0, factor1=1.0, factor2=0.0, power=1.0,
                      datum=10.0, elevmin=1.0, elevmax=2.0)
    with pytest.raises(ParameterError, match="be >= elevmin") as exc:
        VolumeToElevation(offset=0.0, factor1=1.0, factor2=0.0, power=1.0,
                          datum=10.0, elevmin=5.0, elevmax=2.0)
    assert exc.value.field == "elevmax"


def test_buffers_must_be_non_negative_and_non_empty():
    InitialState(vol=0.0, drain_buffer=(0.0, 1.0, 2.0))
    with pytest.raises(ParameterError) as exc:
        InitialState(vol=0.0, drain_buffer=(1.0, -1.0))
    assert exc.value.field == "drain_buffer[1]"
    with pytest.raises(ParameterError, match="at least one entry"):
        InitialState(vol=0.0, macro_buffer=())


def test_the_model_checks_what_no_single_object_can():
    upper = UpperStore(**GOOD_UPPER)
    Model(upper=upper, initial=InitialState(vol=0.3))
    with pytest.raises(ParameterError, match="be <= UpperStore.maxvol") as exc:
        Model(upper=upper, initial=InitialState(vol=0.4))
    assert exc.value.field == "vol"

    with pytest.raises(ParameterError, match="no lower store"):
        Model(upper=upper, initial=InitialState(vol=0.1, vol_lower=0.2))

    lower = LowerStore(maxvol=0.5, gamma=1.0, ks=0.01, m=0.5, l=0.5, extravol=0.1)
    Model(upper=upper, lower=lower, initial=InitialState(vol=0.1, vol_lower=0.6))
    with pytest.raises(ParameterError, match="maxvol \\+ extravol"):
        Model(upper=upper, lower=lower, initial=InitialState(vol=0.1, vol_lower=0.7))


def test_lower_bucket_elevation_needs_a_lower_store():
    upper = UpperStore(**GOOD_UPPER)
    elev = VolumeToElevation(offset=0.0, factor1=1.0, factor2=0.0, power=1.0,
                             datum=10.0, bucket="lower")
    with pytest.raises(ValueError, match="needs a lower store"):
        Model(upper=upper, elevation=elev)
    Model(upper=upper, lower=LowerStore(maxvol=0.5, gamma=1.0, ks=0.01, m=0.5, l=0.5),
          elevation=elev)


def test_wrong_object_types_are_caught_at_construction():
    with pytest.raises(ValueError, match="must be a UpperStore"):
        Model(upper="not a store")
    with pytest.raises(ValueError, match="must be a LowerStore or None"):
        Model(upper=UpperStore(**GOOD_UPPER), lower=object())


def test_with_replaces_one_field():
    model = Model(upper=UpperStore(**GOOD_UPPER), solver=Solver(nstep=3))
    other = model.with_(solver=Solver(nstep=9))
    assert model.solver.nstep == 3 and other.solver.nstep == 9
    assert other.upper is model.upper
