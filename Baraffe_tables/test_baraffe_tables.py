import pytest
import numpy as np
from BDmass_to_flux_ratio import flux_mag_ratio, calculate_flux_ratios





# Test Mass to flux ratio
def test_flux_mag_ratio():
    """Test flux-magnitude ratio.

    A magnitude difference of 5 should be 100 (by definition).

    """
    # positive 5 magnitude difference
    assert np.allclose(flux_mag_ratio(1, 6), 100)

    # negative 5 magnitude difference
    assert np.allclose(flux_mag_ratio(7, 2), 1/100)



def test_calculate_flux_ratios():
    """Returns a dict with the 3 overlapping values."""
    star_params = {"FLUX_J": 1, "FLUX_H": 1, "FLUX_K": 2}   # Names from simbad
    companion_params = {"Mj": 6, "Mh": 11, "Mk": 7}      # Names from Baraffe

    flux_ratios = calculate_flux_ratios(star_params, companion_params)
    assert isinstance(flux_ratios, dict)

    assert np.allclose(flux_ratios["J"], 100)
    assert np.allclose(flux_ratios["H"], 10000)
    assert np.allclose(flux_ratios["K"], 100)



# Test Flux ratio to Mass
