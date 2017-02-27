"""Test companion/flux ratio codes with baraffe tables."""
import sys
import pytest
import numpy as np
from astropy.constants import M_jup, M_sun

from calculations import flux_mag_ratio, calculate_flux_ratio, calculate_companion_magnitude
from table_search import mass_table_search, magnitude_table_search, age_table
from BDmass_to_flux_ratio import main as mass_main
from BDmass_to_flux_ratio import _parser as mass_parser
from flux_ratio_to_BDmass import main as ratio_main
from flux_ratio_to_BDmass import _parser as ratio_parser
from db_queries import get_sweet_cat_temp, get_stellar_params, get_temperature

org_sysargv = sys.argv


# Test two main functions.
@pytest.mark.xfail  # Query fails when offline
def test_main(raises=Exception):
    """Check it returns 0 (Runs normally).

    If there is no internet then an Exception is raised.
        Exception: Query failed: HTTPConnectionPool(host='simbad.u-strasbg.fr', port=80): Max retries exceeded with
        url: /simbad/sim-script (Caused by NewConnectionError('<requests.packages.urllib3.connection.HTTPConnection
        object at 0x7fe4e509b438>: Failed to establish a new connection: [Errno -3] Temporary failure in name
        resolution',)).
    """
    assert mass_main("HD30501", 90, 5) is 0
    assert mass_main("HD30501", 90, 5, area_ratio=True) is 0


@pytest.mark.xfail  # Query fails when offline
def test_ratio_main(raises=Exception):
    """Check it returns 0 (Runs normally).

    If there is no internet then an Exception is raised.
        Exception: Query failed: HTTPConnectionPool(host='simbad.u-strasbg.fr', port=80): Max retries exceeded with
        url: /simbad/sim-script (Caused by NewConnectionError('<requests.packages.urllib3.connection.HTTPConnection
        object at 0x7fe4e509b438>: Failed to establish a new connection: [Errno -3] Temporary failure in name
        resolution',)).
    """
    assert ratio_main("HD30501", 0.01, 5) is 0


def test_flux_mag_ratio():
    """Test flux-magnitude ratio.

    A magnitude difference of 5 should be 100 (by definition).

    """
    # positive 5 magnitude difference
    assert np.allclose(flux_mag_ratio(1, 6), 100)

    # negative 5 magnitude difference
    assert np.allclose(flux_mag_ratio(7, 2), 1. / 100)


def test_calculate_flux_ratio():
    """Test return of a dict with the 3 overlapping values."""
    star_params = {"FLUX_J": 1, "FLUX_H": 1, "FLUX_K": 2}   # Names from simbad
    companion_params = {"Mj": 6, "Mh": 11, "Mk": 7}      # Names from Baraffe

    flux_ratios = calculate_flux_ratio(star_params, companion_params, bands=["J", "H", "K"])
    assert isinstance(flux_ratios, dict)

    assert np.allclose(flux_ratios["J"], 100)
    assert np.allclose(flux_ratios["H"], 10000)
    assert np.allclose(flux_ratios["K"], 100)


# @pytest.mark.xfail
def test_get_sweet_cat_temp():
    """Test getting from sweet-cat."""
    # hd number in SweetCat
    a = get_sweet_cat_temp("HD107383")
    assert isinstance(a, float)
    assert a == 4830
    # hd number not in sweet
    b = get_sweet_cat_temp("HD10")
    assert b is False
    # non hd id
    with pytest.raises(NotImplementedError):
        get_sweet_cat_temp("GJ 422")  # in sweetcat but not an hd number

    # 2 tries that line in sweetcat that has no temperature
    # These are the only two that have hd numbers
    c = get_sweet_cat_temp("HD145934")  # This may change if sweetcat is updated
    assert c is False
    d = get_sweet_cat_temp("HD41004B")  # This may change if sweetcat is updated
    assert d is False


@pytest.mark.xfail(raises=Exception)
def test_get_temperature_without_params_input():
    """Test it calls for params itself, if needed.

    If there is no internet then an Exception is raised.
        Exception: Query failed: HTTPConnectionPool(host='simbad.u-strasbg.fr', port=80): Max retries exceeded with
        url: /simbad/sim-script (Caused by NewConnectionError('<requests.packages.urllib3.connection.HTTPConnection
        object at 0x7fe4e509b438>: Failed to establish a new connection: [Errno -3] Temporary failure in name
        resolution',)).
    """
    name = "HD30501"
    params = get_stellar_params(name)
    assert get_temperature(name, star_params=None) == get_temperature(name, star_params=params)


@pytest.mark.xfail(raises=Exception)
def test_get_temerature_examples():
    """Test some temperatures.

    If there is no internet then an Exception is raised.
        Exception: Query failed: HTTPConnectionPool(host='simbad.u-strasbg.fr', port=80): Max retries exceeded with
        url: /simbad/sim-script (Caused by NewConnectionError('<requests.packages.urllib3.connection.HTTPConnection
        object at 0x7fe4e509b438>: Failed to establish a new connection: [Errno -3] Temporary failure in name
        resolution',)).
    """
    star = "HD12116"   # Has zero temp in simbad (Check this is not used)
    assert get_temperature(star) != 0

    star_simbad = "HD215909"    # has a temperature in simbad
    simbad_temp = 4328
    assert get_temperature(star_simbad) == simbad_temp

    # Sweet_cat temperaute and NO Simbad temperature!
    sweet_name = "HD343246"
    sweet_temp = 5754
    assert get_temperature(sweet_name) == sweet_temp


@pytest.mark.xfail(raises=Exception)
def test_get_stellar_params():
    """Test some values from simbad database reult.

    If there is no internet then an Exception is raised.
        Exception: Query failed: HTTPConnectionPool(host='simbad.u-strasbg.fr', port=80): Max retries exceeded with
        url: /simbad/sim-script (Caused by NewConnectionError('<requests.packages.urllib3.connection.HTTPConnection
        object at 0x7fe4e509b438>: Failed to establish a new connection: [Errno -3] Temporary failure in name
        resolution',)).
    """
    name = "HD219828"
    params = get_stellar_params(name)

    # Check Name put into params
    assert params["name"] == name
    # Check some values
    assert params["Fe_H_Teff"][0] == 5842
    assert params["FLUX_B"] == 8.68
    assert params["FLUX_V"] == 8.01
    assert params["FLUX_K"] == 6.530
    assert params["PLX_VALUE"][0] == 13.83  # parallax
    assert params['Fe_H_log_g'] == 4.19     # log g
    assert params['Fe_H_Fe_H'] == 0.16      # metalicity


# Test Flux ratio to Mass
@pytest.mark.parametrize("band", ["H", "J", "K"])
def test_mag_conversions(band):
    """Test converting from flux ratio to magnitude back to flux ratio etc.

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

    assert np.allclose(found_mass, start_sol_mass)
    assert np.allclose(found_mass * M_sun / M_jup, starting_mass)  # back int M_jup


@pytest.mark.parametrize("age", [0.1, 1, 10])
@pytest.mark.parametrize("model", ["2003", "2015"])
def test_age_table(age, model):
    """Select a baraffe table of certian age and model."""
    model_table, cols = age_table(age, model=model)

    assert isinstance(model_table, dict)
    assert isinstance(cols, list)

    for key in model_table:
        assert isinstance(model_table[key], np.ndarray)

    # check common columns are present
    common_cols = ["M/Ms", "Teff", "L/Ls", "g", "Mv",
                   "Mr", "Mi", "Mj", "Mh", "Mk", "Mll", "Mm"]
    for header in common_cols:
        assert header in cols


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

    # Test error
    bands = ["H", "P"]
    with pytest.raises(ValueError):
        calculate_companion_magnitude(star_vals, 0.001, bands)


def test_mag_table_search_band():
    """One one band value is allowed."""
    age = 5
    magnitudes = {"H": 1, "J": 4, "K": 5}

    with pytest.raises(ValueError):
        magnitude_table_search(magnitudes, age, band=["H", "J", "K"], model="2003")


def test_mass_table_search_03():
    """That a value from the table returns the correct row."""
    mass = 0.09
    comp_params = mass_table_search(mass, 5, model="2003")
    print(comp_params)
    assert comp_params["M/Ms"] == 0.09
    assert comp_params["Teff"] == 2622
    assert comp_params["R"] == 0.113
    assert comp_params["Mk"] == 10.04


def test_mass_table_search_15():
    """Manual test of a 2015 table mass search."""
    mass = 0.09
    comp_params_15 = mass_table_search(mass, 5, model="2015")
    print(comp_params_15)
    assert comp_params_15["M/Ms"] == 0.09
    assert comp_params_15["Teff"] == 2644
    assert comp_params_15["R/Rs"] == 0.113
    assert comp_params_15["Mk"] == 9.91


def test_magnitude_table_search_03():
    """That a value from the table returns the correct row."""
    mag_params = magnitude_table_search({"K": 10.04}, 5, band="K", model="2003")
    print("mag_params", mag_params)
    assert mag_params["M/Ms"] == 0.09
    assert mag_params["Teff"] == 2622
    assert mag_params["R"] == 0.113
    assert mag_params["Mk"] == 10.04


def test_magnitude_table_search_15():
    """Manual test of a 2015 table magnitude search."""
    mag_params_15 = magnitude_table_search({"K": 9.91}, 5, band="K", model="2015")
    assert mag_params_15["M/Ms"] == 0.09
    assert mag_params_15["Teff"] == 2644
    assert mag_params_15["R/Rs"] == 0.113
    assert mag_params_15["Mk"] == 9.91


def test_magnitude_table_search_errors():
    """Test only takes one band as a str."""
    # Check single band as list is ok.
    with pytest.raises(ValueError):
        # More than one band not allowed
        magnitude_table_search({"K": 9.91, "H": 9.91}, 5, band=["H", "K"], model="2015")

    with pytest.raises(ValueError):
        # More than one band not allowed
        magnitude_table_search({"K": 9.91, "H": 9.91}, 5, band=("K",), model="2015")

    with pytest.raises(ValueError):
        # More then one band not allowed
        magnitude_table_search({"K": 9.91, "H": 9.91}, 5, band=["K"], model="2015")

    with pytest.raises(ValueError):
        # Bad band (Not present on cols)
        magnitude_table_search({"K": 9.91, "H": 9.91}, 5, band="Z", model="2015")

    with pytest.raises(ValueError):
        # Band not in magnitudes given
        magnitude_table_search({"K": 9.91, "H": 9.91}, 5, band="J", model="2015")


# Need sys.argv fixture with teardown
def test_BDmass_parser():
    """Test argparse function using sys.argv."""
    sys.argv = []
    test_args = "pytest HD30501 90 5 -a".split()
    for arg in test_args:
        sys.argv.append(arg)
    args = mass_parser()
    assert args.star_name == "HD30501"
    assert args.companion_mass == 90
    assert args.age == 5
    assert args.bands == ["All"]
    assert args.model == "2003"
    assert args.area_ratio is True
    # Not the correct way for teardown
    sys.argv = org_sysargv


# Need sys.argv fixture with teardown
def test_ratio_parser():
    """Test argparse function using sys.argv."""
    sys.argv = []
    test_args = "pytest HD30501 0.001 5 -b H -m 2015".split()
    for arg in test_args:
        sys.argv.append(arg)

    args = ratio_parser()
    assert args.star_name == "HD30501"
    assert args.flux_ratio == 0.001
    assert args.bands == ["H"]
    assert args.age == 5
    assert args.model == "2015"
    sys.argv = org_sysargv


def test_ratio_parser2():
    """Test argparse function using sys.argv more than one band."""
    sys.argv = []
    test_args = "pytest HD30501 0.001 5 -m 2015 -b H K".split()
    for arg in test_args:
        sys.argv.append(arg)

    args = ratio_parser()
    assert args.star_name == "HD30501"
    assert args.flux_ratio == 0.001
    assert args.bands == ["H", "K"]
    assert args.age == 5
    assert args.model == "2015"
    sys.argv = org_sysargv


# Need sys.argv fixture with teardown
def test_failing_ratio_parsers():
    """Test argparse function using sys.argv."""
    sys.argv = []
    test_args = "pytest HD30501 ratio 5 -b H K -m 2015".split()
    for arg in test_args:
        sys.argv.append(arg)

    with pytest.raises(SystemExit):
        ratio_parser()   # ratio is not a number

    sys.argv = []
    test_args = "pytest HD30501 0.001 5 -b H K -m 2015 -z 4".split()  # unknown parameter -z
    for arg in test_args:
        sys.argv.append(arg)

    with pytest.raises(SystemExit):
        ratio_parser()   # unknown parameter -z

    # Invalid choice of band
    sys.argv = []
    test_args = "pytest HD305010 0.01 5 -b Q".split()  # invlaid choice Q
    for arg in test_args:
        sys.argv.append(arg)

    with pytest.raises(SystemExit):
        ratio_parser()   # SystemExit

    sys.argv = org_sysargv


# Need sys.argv fixture with teardown
# @pytest.mark.xfail
def test_failing_mass_parsers():
    """Test argparse function using sys.argv."""
    sys.argv = []
    test_args = "pytest HD30501 mass 5 -b H K -m 2015".split()
    for arg in test_args:
        sys.argv.append(arg)

    with pytest.raises(SystemExit):
        mass_parser()   # ratio is not a number

    sys.argv = []
    test_args = "pytest HD30501 mass 5 -b H K -m 2003-z 4".split()  # unknown parameter -z
    for arg in test_args:
        sys.argv.append(arg)

    with pytest.raises(SystemExit):
        mass_parser()   # unknown parameter -z

    # Invalid choice of band
    sys.argv = []
    test_args = "pytest HD305010 100 5 -b Q".split()  # invlaid choice Q
    for arg in test_args:
        sys.argv.append(arg)

    with pytest.raises(SystemExit):
        mass_parser()   # SystemExit

    sys.argv = org_sysargv
