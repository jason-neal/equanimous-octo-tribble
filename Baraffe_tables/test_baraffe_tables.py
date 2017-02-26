import pytest
import numpy as np
from astropy.constants import M_jup, M_sun

from calculations import flux_mag_ratio, calculate_flux_ratio, calculate_companion_magnitude
from table_search import mass_table_search, magnitude_table_search, age_table
from BDmass_to_flux_ratio import main as mass_main
from flux_ratio_to_BDmass import main as ratio_main
from db_queries import get_sweet_cat_temp


# Test Mass to flux ratio
def test_flux_mag_ratio():
    """Test flux-magnitude ratio.

    A magnitude difference of 5 should be 100 (by definition).

    """
    # positive 5 magnitude difference
    assert np.allclose(flux_mag_ratio(1, 6), 100)

    # negative 5 magnitude difference
    assert np.allclose(flux_mag_ratio(7, 2), 1. / 100)


def test_calculate_flux_ratio():
    """Returns a dict with the 3 overlapping values."""
    star_params = {"FLUX_J": 1, "FLUX_H": 1, "FLUX_K": 2}   # Names from simbad
    companion_params = {"Mj": 6, "Mh": 11, "Mk": 7}      # Names from Baraffe

    flux_ratios = calculate_flux_ratio(star_params, companion_params, bands=["J", "H", "K"])
    assert isinstance(flux_ratios, dict)

    assert np.allclose(flux_ratios["J"], 100)
    assert np.allclose(flux_ratios["H"], 10000)
    assert np.allclose(flux_ratios["K"], 100)


@pytest.mark.xfail  # Query fails when offline
def test_main():
    # Check it returns 0 (Runs normally)
    assert mass_main("HD30501", 90, 5) is 0


@pytest.mark.xfail
def test_get_sweet_cat_temp():
    """Test getting from sweet-cat."""
    # hd number in SweetCat
    a = get_sweet_cat_temp("HD107383")
    assert isinstance(a, float)
    assert a == 4830
    # hd number not in sweet
    b = get_sweet_cat_temp("HD1")
    assert b is False
    # non hd id
    with pytest.raises(NotImplementedError):
        get_sweet_cat_temp("GJ 422")  # in sweetcat but not an hd number


# Test Flux ratio to Mass
@pytest.mark.xfail  # Query fails when offline
def test_ratio_main():
    # Check it returns 0 (Runs normally)
    assert ratio_main("HD30501", 0.01, 5) is 0


@pytest.mark.parametrize("band", ["H", "J", "K"])
def test_mag_conversions(band):
    """ Test converting from flux ratio to magnitude back to flux ratio etc.

    Tests show the conversion goes back and forward.
    """
    comp_vals = {"Mj": 4, "Mk": 6, "Mh": 12}
    star_vals = {"FLUX_J": 3, "FLUX_K": 11, "FLUX_H": 5}

    ratios = calculate_flux_ratio(star_vals, comp_vals, band)
    print("Ratios from function", ratios)
    magnitudes = calculate_companion_magnitude(star_vals, 1. / ratios[band], band)
    print("magnitudes from ratios", magnitudes)
    magnitudes["M{}".format(band.lower())] = magnitudes[band]
    print("magnitudes from ratios", magnitudes)
    new_ratios = calculate_flux_ratio(star_vals, magnitudes, band)
    print("new_ratios from mags", new_ratios)

    assert np.allclose(new_ratios[band], ratios[band])


# @pytest.mark.xfail()
@pytest.mark.parametrize("mass_model", [(50, "2003"), (150, "2015")])
@pytest.mark.parametrize("age", [1, 5, 10])
@pytest.mark.parametrize("band", ["H", "J", "K"])
# @pytest.mark.parametrize("model", ["2003", "2015"])
def test_table_searches(mass_model, age, band):
    """That a given mass calculates a magnitude and that mag finds a correct mass."""
    starting_mass = mass_model[0]  # M_jup
    model = mass_model[1]
    start_sol_mass = starting_mass * (M_jup / M_sun).value
    print("Starting_mass = {0} (Mjup), {1} (Msun)".format(starting_mass, start_sol_mass))

    comp_params = mass_table_search(start_sol_mass, age=age, model=model)
    found_mag = comp_params["M{0!s}".format(band.lower())]
    print("Band", band, "Found magnitude", found_mag)
    magnitude = {band: found_mag}
    print("magnitude", magnitude)

    mag_params = magnitude_table_search(magnitude, age=age, band=band, model=model)
    print("mag_params", mag_params)

    found_mass = mag_params["M/Ms"]
    print("found_mass", found_mass, "start_mass", start_sol_mass)

    assert np.allclose(found_mass, start_sol_mass, rtol=-2)
    assert np.allclose(found_mass * M_sun / M_jup, starting_mass)  # back int M_jup


@pytest.mark.parametrize("age", [0.1, 1, 10])
@pytest.mark.parametrize("model", ["2003", "2015", pytest.mark.xfail("2016")])
def test_age_table(age, model):
    """Selects a baraffe table of certian age and model."""
    model_03, cols_03 = age_table(5, "2003")

    model_15, cols_15 = age_table(5, "2015")

    assert isinstance(model_03, np.ndarray)
    assert isinstance(cols_03, list)
    assert isinstance(model_15, np.ndarray)
    assert isinstance(cols_15, list)

    # Assert actually different
    assert np.all(cols_15 != cols_03)
    assert np.all(model_15 != model_03)
    assert len(cols_15) != len(cols_03)
    assert len(model_15) != len(model_03)

    # check common columns are present
    common_cols = ["M/Ms", "Teff", "L/Ls", "g", "Mv",
                   "Mr", "Mi", "Mj", "Mh", "Mk", "Mll", "Mm"]
    for header in common_cols:
        assert header in cols_15
        assert header in cols_03


@pytest.mark.parametrize("model", ["2016", "", 2003, 15, "word"])
def test_bad_model_age_table(model):
    """Call model with wrong values."""
    with pytest.raises(ValueError):
        age_table(5, model)


def test_calculate_comp_magnitude():
    """Test t returns values fr all bands given."""
    bands = ["H", "J", "K"]
    star_vals = {"FLUX_J": 1, "FLUX_H": 1, "FLUX_K": 2}   # Names from simbad
    magnitudes = calculate_companion_magnitude(star_vals, 0.001, bands)

    assert isinstance(magnitudes, dict)
    assert len(magnitudes) == 3
    for band in bands:
        assert band in magnitudes.keys()


@pytest.mark.xfail
def test_mag_table_search_band():
    """One one band value is allowed."""
    age = 5
    magnitudes = {"H": 1, "J": 4, "K": 5}

    with pytest.raises(ValueError):
        magnitude_table_search(magnitudes, age, band=["H", "J", "K"], model="2003")


def test_mass_table_search_03():
    """That a value from the table returns the correct row."""
    comp_params = mass_table_search(0.09, 5, model="2003")
    print(comp_params)
    assert comp_params["M/Ms"] == 0.09
    assert comp_params["Teff"] == 2622
    assert comp_params["R"] == 0.113
    assert comp_params["Mk"] == 10.04

def test_mass_table_search_15():
    comp_params_15 = mass_table_search(0.09, 5, model="2015")
    print(comp_params_15)
    assert comp_params_15["M/Ms"] == 0.09
    assert comp_params_15["Teff"] == 2644
    assert comp_params_15["R/Rs"] == 0.113
    assert comp_params_15["Mk"] == 9.91


def test_magnitude_table_search_03():
    """That a value from the table returns the correct row."""
    mag_params = magnitude_table_search({"K": 9.91}, 5, band="K", model="2003")
    print("mag_params", mag_params)
    assert mag_params["M/Ms"] == 0.09
    assert mag_params["Teff"] == 2622
    assert mag_params["R"] == 0.113
    assert mag_params["Mk"] == 10.04


def test_magnitude_table_search_15():
    mag_params_15 = magnitude_table_search({"K": 10.04}, 5, band="K", model="2015")
    print("mag_params", mag_params_15)
    assert mag_params_15["M/Ms"] == 0.09
    assert mag_params_15["Teff"] == 2644
    assert mag_params_15["R/Rs"] == 0.113
    assert mag_params_15["Mk"] == 9.91
