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

def snrscan(Spec, width):
    snr_list = []
    for i in range(len(Spec)-width):
        snr = 1 / np.std(Spec[i:i+width])
        snr_list.append(snr)

    return snr_list




def main(fname, binsize=5):

#load the fits file 
    data = fits.getdata(fname)
    hdr = fits.getheader(fname)
    try:
        test = data["Extracted_OPT"]
    else:
        test = data["Extracted_DRACS"]
#test = [1,2,3,4,2,1,2,1,1,2,3,2,1,4,2,1,2,3,4,1,2,3,4,1,2,4,2,1,2,3,2,2.5,2,2,2,1,1,2,3,4]
#tests = [t / 100 for t in test]
    x = range(len(test))
    x5 = [xi + 2.5 for xi in x]
    x10 = [xi + 5 for xi in x]
    #print(snrscan(tests,10))
    #print(tests)
    print(len(x),len(x10),len(x5))

    snrlist10 = snrscan(test,10)
    snrlist5 = snrscan(test,5)
    plt.figure
    plt.plot(x, test, label="data")
    plt.plot(x10[:-10], snrlist10, label="snr10")
    plt.plot(x5[:-5], snrlist5, label="snr5")
    plt.legend()
    plt.show()



if __name__ == "__main__":
    args = vars(_parser())
    fname = args.pop('fname')
    opts = {k: args[k] for k in args}

    main(fname, **opts)
