#!/usr/bin/env python

# relative_rv.py

# Script to calculate the relative Radial Velocity between two given wavelengths

# Jason Neal
# December 2016
from __future__ import division, print_function
import argparse


def relative_rv(wav_1, wav_2):
    """Calculate the radial velocity difference between two wavelength values."""
    c = 299792.458   # km / s
    difference = wav_2 - wav_1
    relative = difference / wav_1
    return relative * c


def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Calculate RV doppler shift between two wavelengths')
    parser.add_argument('original_wavelength', help='The original wavelength', type=float)
    parser.add_argument('new_wavelength', help='The new wavelength', type=float)
    parser.add_argument('-m', '--mode', choices=['normal', 'reverse', 'iterable'], default='normal',
                        help='Output mode', type=str)   # reverse, both, (iteratorable input)
    # TODO: Add these extra modes
    parser.add_argument('-u', '--unit', help='Distance unit of measurement [Mm, km, m, cm, mm]',
                        choices=["Mm", "km", "m", "cm", "mm"], default="km")
    # parser.add_argument('-r', '--round', help='Rounding on ouput, e.g. ["2dp", "3sf"]')
    args = parser.parse_args()
    return args


def main(original_wavelength, new_wavelength, mode=False, unit="km"):  # unit, rounding
    """Obtian RV offset between wavelengths and print to screen.

    Prints rv between the wavelength values.
    """
    if mode == "reverse":
        # Switch direction of doppler shift
        original_wavelength, new_wavelength = new_wavelength, original_wavelength

    rv = relative_rv(original_wavelength, new_wavelength)

    if unit == "mm":
        rv *= 1e6
    elif unit == "cm":
        rv *= 1e5
    elif unit == "m":
        rv *= 1e3
    elif unit == "km":
        rv *= 1e0
    elif unit == "Mm":
        rv *= 1e-3

    print("\n" + "-" * 20)
    print("Original   -->   New   |  RV ({0}/s)".format(unit))
    print("-" * 20)
    print("{0}  -->  {1}  |  {2}".format(original_wavelength, new_wavelength, rv))
    print("-" * 20)


if __name__ == "__main__":
    args = vars(_parser())
    wave1 = args.pop('original_wavelength')
    wave2 = args.pop('new_wavelength')
    opts = {k: args[k] for k in args}

    main(wave1, wave2, **opts)
