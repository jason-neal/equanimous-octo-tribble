#!/usr/bin/env python

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import argparse
from astropy.io import fits
# Testing for a snr of my spectra

def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Telluric Removal')
    parser.add_argument('fname', help='Input fits file')
    #parser.add_argument('-x', '--export', default=False,
    #                    help='Export result to fits file True/False')
    parser.add_argument('-b', '--binsize', default=5,
                        help='Snr binsize')
    #parser.add_argument('-k', '--kind', default="linear",
    #                    help='Interpolation order, linear, quadratic or cubic')
    #parser.add_argument('-m', '--method', default="scipy",
    #                    help='Interpolation method numpy or scipy')
    args = parser.parse_args()
    return args

def snr(spectra, ddof=0):
	snr = np.mean(spectra) / np.std(spectra, ddof=ddof)
	return snr

def snrscan(Spec, width):
    snr_list = []
    for i in range(len(Spec)-width):
        #snr = np.mean(Spec[i:i+width]) / np.std(Spec[i:i+width], ddof=0)
	SNR = snr(Spec[i:i+width])
        snr_list.append(SNR)

    return snr_list

def click_snr(wl, Spec):
    """ Calculate snr in a specific range given by clicks on a plot """
    plt.figure()
    plt.plot(wl, Spec)
    plt.show(block=False)
    # points from click

    # temp values untill implement above
    point2 = np.max(wl) 
    point1 = np.min(wl) 
    map2 = wl < point2
    map1 = wl > point1
    wl_slice = wl[map1*map2]
    Spec_slice = Spec[map1*map2]
 
    # Calculate SNR on the slice
    SNR = snr(Spec_slice)

    return SNR


def main(fname, binsize=5):
    binsize = int(binsize)
    # load the fits file 
    data = fits.getdata(fname)
    hdr = fits.getheader(fname)
    try:
        test = data["Extracted_OPT"]
    except:
        test = data["Extracted_DRACS"]
#test = [1,2,3,4,2,1,2,1,1,2,3,2,1,4,2,1,2,3,4,1,2,3,4,1,2,4,2,1,2,3,2,2.5,2,2,2,1,1,2,3,4]
#tests = [t / 100 for t in test]
    x = range(len(test))

    xbin = [xi + binsize/2 for xi in x]  # center snr values
    xbin = xbin[:-binsize]
    
    snrlist = snrscan(test, binsize)

    plt.figure
    plt.plot(x, test, label="data")
    plt.plot(xbin, snrlist, label="snr sze"+ str(binsize))
    plt.legend()
    plt.show()



if __name__ == "__main__":
    args = vars(_parser())
    fname = args.pop('fname')
    opts = {k: args[k] for k in args}

    main(fname, **opts)
