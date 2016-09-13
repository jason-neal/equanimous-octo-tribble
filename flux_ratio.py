#!/usr/bin/python

"""Brown Dwarf Flux ratio calculator:

Inputs :
Star name:
Mass of companion in Jupiter masses
Age of star: """

from __future__ import division, print_function
import numpy as np

def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Wavelength Calibrate Spectra')
    
    parser.add_argument('star_name', help='Input fits file to calibrate')

    parser.add_argument('companion_mass', help='Mass of companion')
    parser.add_argument('age', help='Star age')
    args = parser.parse_args()
    return args


def main(star_name, companion_mass, stellar_age):
    """Compute flux ratio of star to companion """
    pass
if __name__ == '__main__':
    args = vars(_parser())
    star_name = args.pop('star_name')
    companion_mass = args.pop('companion_mass')
    Age = args.pop('Age')
    opts = {k: args[k] for k in args}

    main(star_name, companion_mass, age, **opts)