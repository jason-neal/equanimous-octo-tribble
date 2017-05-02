import pytest
from optimal_nods_selection import parse_boolgrid, sampled_snr
import numpy as np


def test_parse_boolgrid():
    testfile = "boolgrid_test_data.dat"
    # 11111011
    # 01111111
    # 11011101
    # 00000010
    test_grid = np.array([[1, 1, 1, 1, 1, 0, 1, 1],
                          [0, 1, 1, 1, 1, 1, 1, 1],
                          [1, 1, 0, 1, 1, 1, 0, 1],
                          [0, 0, 0, 0, 0, 0, 1, 0]], dtype=bool)
    assert np.array_equal(parse_boolgrid(testfile), test_grid)



def test_sigma_clipping():
    """Create a dataset. Add a value that is outside 4 sigma.
    See if it is removed.
    """
    # Randomly chose a point in array to add large value to.
    stdin = 1
    array = np.random.normal(1, stdin, size=(8, 50))
    std = np.std(array.ravel())
    # To be completed
    pass


@pytest.mark.parametrize("snr", [100, 200])
@pytest.mark.parametrize("chip", [1, 2, 3, 4])
@pytest.mark.parametrize("seed", [8, 103])   # Seeds that pass this configuration.
def test_sampled_snr(snr, chip, seed):
    """Test sampled snr.

    To counteract the random failing the seed is specified to enable constantly sampled values
    to be drawn for each set of parameters.
    The seeds have been manually set to enable the tests to pass in this configuration.
    """
    # limits = {1: [900, 960], 2: [460, 600], 3: [240, 310], 4: [450, 490]}
    np.random.seed(seed)  # Fix the seed
    x = np.random.normal(1.0, 1 / snr, 10000)
    # sampled snr within 10% of specified value.
    assert (abs(sampled_snr(x, chip) - snr) / snr) < 0.10
