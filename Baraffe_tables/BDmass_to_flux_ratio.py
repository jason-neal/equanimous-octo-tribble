#!/usr/bin/python
"""Brown Dwarf Flux ratio calculator.

Calculates the flux/contrast ratio between a host star and a brown dwarf of a specified mass.

This script uses the SIMBAD database to obtain the host star parameters, such as magnitude and age.
The companion/brown dwarf mass is a given input (in Mjup) and  is used to obtain the band magnitudes
of the companion from the Barraffe tables.

The magnitude difference between the host and companion are used to caluate the flux/contrast ratio.

Inputs
------
Star name: str
    Stellar idenification number. eg. HD30501
companion_mass: float
    Mass of companion in Jupiter masses
age: float
    Stellar Age. (Closest model is used)
bands: list of str
    Spectral bands to obtain ratio.
model: str
    Choose between the 2003 and 2015 Barraffe modeling.

"""
# TODO: Interpolate between tables?
from __future__ import division, print_function
import sys
import argparse
import numpy as np
from typing import List, Optional

from calculations import calculate_flux_ratio, calculate_stellar_radius
from table_search import mass_table_search
from db_queries import get_stellar_params


def _parser() -> object:
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Determine flux ratio of stellar companion')
    parser.add_argument('star_name', help='Input fits file to calibrate')
    parser.add_argument('companion_mass', help='Mass of companion (M_Jup)', type=float)
    parser.add_argument('age', help='Star age (Gyr)', type=float)
    parser.add_argument('-b', '--bands', help='Spectral Band to measure. Options=["All", "K", ""]',
                        choices=["All", "J", "H", "K"], default=["All"], nargs="+", type=str)
    parser.add_argument('-m', '--model', choices=['03', '15', '2003', '2015'],
                        help='Baraffe model to use [2003, 2015]', default='2003', type=str)
    parser.add_argument("-a", "--area_ratio", help="Calculate the area ratio.", default=False, action="store_true")
    args = parser.parse_args()
    return args


def main(star_name: str, companion_mass: float, stellar_age: float, bands: Optional[List[str]]=None,
         model: str="2003", area_ratio: bool=False) -> int:
    """Compute flux/contrast ratio between a stellar host and companion.

    Parameters
    ----------
    star_name: str
        Stellar idenification number. eg. HD30501.
    companion_mass: float
        Mass of companion in Jupiter masses.
    stellar_age: float
        Stellar Age. (Closest model is used).
    bands: list of str
        Spectral bands to obtain ratio.
    model: str (optional)
        Year of Barraffe model to use [2003 (default), 2015].
    area_ratio: bool default=False
        Perform simple radius and area comparions calculations.

    """
    if (bands is None) or ("All" in bands):
        bands = ["J", "H", "K"]

    # Obtain Stellar parameters from astroquery
    star_params = get_stellar_params(star_name)  # returns a astroquesry result table

    companion_mass_solar = companion_mass * (Mjup / Msun).value    # transform to solar mass for table search
    # Get parameters for this mass and age
    companion_params = mass_table_search(companion_mass_solar, stellar_age, model=model)
    flux_ratios = calculate_flux_ratio(star_params, companion_params, bands)

    # Print flux ratios using a generator
    print("\nFlux ratios:")
    print_generator = (("{0!s} band star/companion Flux ratio = {1:4.2f},"
                       " >>> companion/star Flux ratio = {2:0.4f}").format(key, val[0], 1. / val[0])
                      for key, val in flux_ratios.items() if key in bands)

    for print_string in print_generator:
        print(print_string)

    if area_ratio:
        # Compare to area ratio
        Rstar = calculate_stellar_radius(star_params)
        print(Rstar)
        Rcomp_Rstar = companion_params["R"] / Rstar

        print("\nRadius Calculation")
        print("Host radius         = {} R_sun".format(Rstar[0]))
        print("companion radius    = {} R_sun".format(np.round(companion_params["R"], 4)))
        print("Radius Ratio of companion/star    = {}".format(Rcomp_Rstar[0]))
        print("Area Ratio of companion/star      = {}".format(Rcomp_Rstar[0]**2))

    return 0


if __name__ == '__main__':
    args = vars(_parser())
    star_name = args.pop('star_name')
    companion_mass = args.pop('companion_mass')
    age = args.pop('age')
    opts = {k: args[k] for k in args}

    sys.exit(main(star_name, companion_mass, age, **opts))
