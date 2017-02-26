import pytest
import numpy as np
from calculations import flux_mag_ratio, calculate_flux_ratios


# Test Mass to flux ratio
def test_flux_mag_ratio():
    """Test flux-magnitude ratio.

    A magnitude difference of 5 should be 100 (by definition).

    """
    # positive 5 magnitude difference
    assert np.allclose(flux_mag_ratio(1, 6), 100)

    # negative 5 magnitude difference
    assert np.allclose(flux_mag_ratio(7, 2), 1. / 100)


def test_calculate_flux_ratios():
    """Returns a dict with the 3 overlapping values."""
    star_params = {"FLUX_J": 1, "FLUX_H": 1, "FLUX_K": 2}   # Names from simbad
    companion_params = {"Mj": 6, "Mh": 11, "Mk": 7}      # Names from Baraffe

    flux_ratios = calculate_flux_ratios(star_params, companion_params)
    assert isinstance(flux_ratios, dict)

    assert np.allclose(flux_ratios["J"], 100)
    assert np.allclose(flux_ratios["H"], 10000)
    assert np.allclose(flux_ratios["K"], 100)


@pytest.mark.xfail  # Query fails when offline
def test_main():
    from BDmass_to_flux_ratio import mass_main
    # Check it returns 0 (Runs normally)
    assert mass_main("HD30501", 90, 5) is 0





# Test Flux ratio to Mass

def test_mag_conversions():
    """ Test converting from flux ratio to magnitude back to flux ratio etc.
    Tests show the conversion goes back and forward."""

    vals = dict()
    star = dict()
    vals["Mj"] = 4
    vals["Mk"] = 6
    vals["Mh"] = 12
    star["FLUX_J"] = 3
    star["FLUX_K"] = 11
    star["FLUX_H"] = 5
    ratios = calculate_flux_ratios(star, vals)
    print("Ratios from function", ratios)
    magnitudes = calculate_companion_magnitude(star, 1. / ratios["J"])
    print("magnitudes from ratios", magnitudes)
    new_ratios = calculate_flux_ratios(star, magnitudes)
    print("new_ratios from mags", new_ratios)

    assert np.allclose(new_ratios, ratios)

