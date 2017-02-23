#!/usr/bin/python

""" Brown Dwarf Mass calculator:

    Uses stellar parameter databases to find host star parameters. The
    magnitude of the low mass companion from the provided flux ratio and the
    coresponding mass is looked up in the Barraffe evolutionary models.

    Inputs
    ------
    Star name: str
        Stellar idenification number. eg. HD30501
    flux_ratio: float
        Flux ratio between host and companion.
    age: float
        Stellar Age. (Closest model is used)
    """

from __future__ import division, print_function
import argparse
from db_queries import get_stellar_params
from calculations import calculate_companion_magnitude
from table_search import magnitude_table_search
# import pandas as pd
# import matplotlib.pyplot as plt


def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Determine mass of stellar ' +
                                     'companion from a flux ratio')
    parser.add_argument('star_name', help='Name of host star.', type=str)
    parser.add_argument('flux_ratio', type=float, help='Flux ratio between host ' +
                        'and companion - (F_companion/F_host)')
    parser.add_argument('age', help='Star age (Gyr)', type=float)
    parser.add_argument("-b", "--band", choices=["J", "K"], default="K",
                        type=str, help='Magnitude band for the flux ratio value')
    parser.add_argument('-m', '--model', choices=['03', '15', '2003', '2015'],
                        help='Baraffe model to use [2003, 2015]',
                        default='2003', type=str)
    args = parser.parse_args()
    return args


def main(star_name, flux_ratio, stellar_age, band="K", model="2003"):
    """Compute companion mass from flux ratio value

    Parameters
    ----------
    star_name: str
        Stellar idenification number. eg. HD30501
    flux_ratio: float
        Flux ratio for the system (F_companion/F_host).
    stellar_age: float
        Age of star/system (Gyr).
    band: str
        Wavelength band to use. (optional)
    model: int (optional)
       Year of Barraffe model to use [2003 (default), 2015].
"""

    Jup_mass = 1047.56  # Jupiters in 1 M_sol

    # Obtain Stellar parameters from astroquery
    star_params = get_stellar_params(star_name)  # returns a astroquesry result table

    # Calculate Magnitude J and K band magnitude for this flux ratio
    magnitudes = calculate_companion_magnitude(star_params, flux_ratio)
    print("Magnitude calculate for companion", magnitudes)
    # Find companion parameters that match these magnitudes
    companion_params = magnitude_table_search(magnitudes, stellar_age,
                                              band=band, model=model)

    # Print flux ratios using a generator
    print("Estimated Companion Mass from {} band Flux ratio\n".format(band.upper()))
    print("M/M_S = {0} (M_star)".format(companion_params["M/Ms"][0]) +
          " = {} (M_Jup)".format(Jup_mass * companion_params["M/Ms"][0]) +
          ", Temp = {} K".format(companion_params["Teff"][0]))



def test_mag_conversions():
    """ Test converting from flux ratio to magnitude back to flux ratio etc.
    Tests show the conversion goes back and forward."""

    vals = dict()
    star = dict()
    vals["Mj"] = 4
    vals["Mk"] = 6
    star["FLUX_J"] = 3
    star["FLUX_K"] = 4
    ratios = calculate_flux_ratios(star, vals)
    print("Ratios from function", ratios)
    magnitudes = calculate_companion_magnitude(star, 1./ratios["J"])
    print("magnitudes from ratios", magnitudes)
    new_ratios = calculate_flux_ratios(star, magnitudes)
    print("new_ratios from mags", new_ratios)
if __name__ == '__main__':
    args = vars(_parser())
    star_name = args.pop('star_name')
    flux_ratio = args.pop('flux_ratio')
    age = args.pop('age')
    opts = {k: args[k] for k in args}

    main(star_name, flux_ratio, age, **opts)
