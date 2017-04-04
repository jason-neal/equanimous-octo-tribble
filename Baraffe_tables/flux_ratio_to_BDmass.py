#!/usr/bin/python
"""Brown Dwarf Mass calculator.

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
import sys
import argparse
from typing import List, Optional
from astropy.constants import M_sun, M_jup

from db_queries import get_stellar_params
from table_search import magnitude_table_search
from calculations import calculate_companion_magnitude


def _parser() -> object:
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Determine mass of stellar ' +
                                     'companion from a flux ratio')
    parser.add_argument('star_name', help='Name of host star.', type=str)
    parser.add_argument('flux_ratio', type=float, help='Flux ratio between host ' +
                        'and companion - (F_companion/F_host)')
    parser.add_argument('age', help='Star age (Gyr)', type=float)
    parser.add_argument("-b", "--bands", choices=["All", "J", "H", "K"], default=["K"],
                        type=str, help='Magnitude bands for the flux ratio value', nargs="+")
    parser.add_argument('-m', '--model', choices=['03', '15', '2003', '2015'],
                        help='Baraffe model to use [2003, 2015]',
                        default='2003', type=str)
    args = parser.parse_args()
    return args


def main(star_name: str, flux_ratio: float, stellar_age: float,
         bands: Optional[List[str]]=None, model: str="2003") -> int:
    """Compute companion mass from flux ratio value.

    Parameters
    ----------
    star_name: str
        Stellar idenification number. eg. HD30501
    flux_ratio: float
        Flux ratio for the system (F_companion/F_host).
    stellar_age: float
        Age of star/system (Gyr).
    bands: str
        Wavelength band to use. (optional)
    model: int (optional)
       Year of Barraffe model to use [2003 (default), 2015].

    """
    Jup_sol_mass = (M_sun / M_jup).value  # Jupiters in 1 M_sol

    if (bands is None) or ("All" in bands):
        bands = ["H", "J", "K"]

    # Obtain Stellar parameters from astroquery
    star_params = get_stellar_params(star_name)  # returns a astroquesry result table

    # Calculate Magnitude J, H and K bands magnitude for this flux ratio
    magnitudes = calculate_companion_magnitude(star_params, flux_ratio, bands=bands)
    print("Magnitude calculation for companion {}".format([magnitudes[key] for key in magnitudes]))

    for band in bands:
        # Find companion parameters that match these magnitudes
        companion_params = magnitude_table_search(magnitudes, stellar_age,
                                                  band=band, model=model)

        # Print flux ratios using a generator
        print("Estimated Companion Mass from {} band Flux ratio\n".format(band.upper()))
        print("M/M_S = {0} (M_star)".format(companion_params["M/Ms"]) +
              " = {} (M_Jup)".format(Jup_sol_mass * companion_params["M/Ms"]) +
              ", Temp = {} K".format(companion_params["Teff"]))

    return 0


if __name__ == '__main__':
    args = vars(_parser())
    star_name = args.pop('star_name')
    flux_ratio = args.pop('flux_ratio')
    age = args.pop('age')
    opts = {k: args[k] for k in args}

    sys.exit(main(star_name, flux_ratio, age, **opts))
