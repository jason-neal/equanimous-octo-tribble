"""cont_norm.py

Based on http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaGuiDoc/continuumFinder.html.
"""
import sys
import argparse

from astropy.io import fits
from PyAstronomy import pyaGui as pg
import numpy as np
import matplotlib.pylab as plt


def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Interactively normalize fits spectra.')

    parser.add_argument("fitsname", type=str, help="Specturm to continuum normalize.")

    parser.add_argument("--suffix", default="inorm", help="Indicator to add to filename before '.fits'.")
    parser.add_argument("-e", "--flux_errors", default=False, action="store_true",
                        help="Calculate normalization fucntion errors on the flux.")
    args = parser.parse_args()
    return args


def main():
    args = _parser()

    data = fits.getdata(args.fitsname)
    hdr = fits.getheader(args.fitsname)
    print(data)
    print(type(data))
    print(data.shape)
    if isinstance(data, fits.fitsrec.FITS_rec):
        wave = data.field(0)
        flux = data.field(1)
    else:
        wave = np.arange(len(data))
        flux = data

    cf = pg.ContinuumInteractive(wave, flux)

    cf.plot([wave[0], wave[-1]], [1, 1], 'k--')

    # Opens the GUI and starts the interactive session.
    c = cf.findContinuum()

    if args.flux_errors:
        flux_error(wave, flux, c)
    # plt.title("Normalized data")
    # plt.plot(wave, c["normalizedData"], 'b.--')
    # plt.show()
    norm_flux = c["normalizedData"]

    fitssave = args.fitsname.replace(".fits", ".{!s}.fits".format(args.suffix))
    if len(data) == 2:
        fits.writeto(fitssave, norm_flux, header=hdr)
    else:
        fits.writeto(fitssave, norm_flux, header=hdr)


def flux_error(wave, flux, c):
    """Take the points in c and perform all the different spline normlazations.
    Calcualate the change in flux based on the difference methods.
    """
    # TODO: This
    raise NotImplementedError


if __name__ == "__main__":
    sys.exit(main())
