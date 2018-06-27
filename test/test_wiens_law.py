from octotribble.wiens_law import wien_displacement, inverse_wien_displacement

from hypothesis import given, strategies as st
import numpy as np


@given(st.floats(min_value=1, max_value=10000))
def test_wien_displacement(temp):
    wave = wien_displacement(temp)
    temp2 = inverse_wien_displacement(wave)
    assert np.allclose(temp, temp2)


@given(st.floats(min_value=0.1, max_value=50))
def test_inverse_wien_displacement(wave):
    temp = inverse_wien_displacement(wave)
    wave2 = wien_displacement(temp)
    assert np.allclose(wave, wave2)

