#!/usr/bin/env python

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import argparse
from astropy.io import fits

from numpy.lib import stride_tricks
import pandas as pd
import seaborn as sns
# Testing for a snr of my spectra

import SpectralTools as s_tools
def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Telluric Removal')
    parser.add_argument('fname', help='Input fits file')
    #parser.add_argument('-x', '--export', default=False,
    #                    help='Export result to fits file True/False')
    parser.add_argument('-i', '--interactive_mode',  action='store_true',
                        help='Interactive mode lets you click the range '
                        + 'to calculate the signal to noise of.')
    parser.add_argument('-b', '--binsize', default=5,
                        help='SNR binsize')
    parser.add_argument('-r', '--bin_range', type=int, nargs=3,
                        help='Values for np.range of bin values [start, stop, step]')
    parser.add_argument('-l', '--limits', type=int, nargs=2,
                        help='Bounds for the wavelength range to '
                            +'calculate snr for in format [min, max]')
    parser.add_argument('-p', '--plot', action='store_true',
                        help="Display plots")
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
    plt.show(block=True)
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

def simple_plot(wav, flux, label="Spectra", xscale="nm"):
    plt.plot(wav, flux, label=label)
    plt.ylabel("FLux")
    if xscale== "nm":
        plt.xlabel("Wavelenght (nm)")
    else:
        plt.xlabel("Pixel number")
    plt.show()

def strided_snr(data, frame_length, hop_length=1):
    num_frames = 1 + (len(data)-frame_length)/hop_length
    row_stride = data.itemsize * hop_length # *hopesize
    col_stride = data.itemsize
    data_strided = stride_tricks.as_strided(data, shape=(num_frames, frame_length), strides=(row_stride, col_stride))
    
    #print("length of data_strided",len(data_strided))
    snr = 1/np.std( data_strided, axis=1)
    #print("frame_length", frame_length)
    #print("num_frames", num_frames)
    #print("len(snr)",len(snr))
    #print(snr)
    
    #zeropad to make uniform length of spectra
    missing_size = len(data)-len(snr)
    #print("missing size", missing_size)
    before = missing_size//2
    end = missing_size//2
    if missing_size%2 is not 0:
        print("missing size is not even !!!! You need to use odd values with even steps")
    padded_snr = np.pad(snr, (before, end),"constant")
    #print("padded length",len(padded_snr))
    #print(padded_snr)
    return padded_snr

def main(fname, binsize=5, interactive_mode=False, bin_range=None, limits=None, plot=False):
    # load the fits file 
    data = fits.getdata(fname)
    hdr = fits.getheader(fname)

    try:
        # Gasgano reduction
        I = data["Extracted_OPT"]
        I = np.array(I, dtype="float64")
    except:
        # Dracs reduction
        try:
            I = data["Extracted_DRACS"]
            print(type(I))
            I = np.array(I, dtype="float64")
            print(type(I))
        except:
            I = data
    try:
        wl = data["Wavelength"]
        xscale = "nm"
    except:
        # no wavelength values so use pixel positions
        wl = range(len(I))
        xscale = "pixel"
        
    if interactive_mode:
        # do interactve mode stuff
        pass

    elif limits: 
        limited_wl, limited_flux = s_tools.fast_wav_selector(wl, I, limits[0], limits[1])
        limit_snr = snr(limited_flux)
        print("SNR of the spectra between the given limits {0}-{1} is {2}".format(limits[0], limits[1], limit_snr))
        if plot:
            simple_plot(wl, I, label="Spectra", xscale=xscale)

    else:
        # Default to scan mode

        binsize = int(binsize)
       
        snr_list = snrscan(I, binsize)
        print("snr list", snr_list)
        # Try using stride on np.array
        plt.hist(snr_list)
        plt.show()
        # striding
        if bin_range:

            bins = np.arange(bin_range[0], bin_range[1], bin_range[2])  # TODO ability for less than 3 args?
            hopper = 1
            print("itemsize", I.itemsize, "dtype", I.dtype)
            store_list = []
            for i, b in enumerate(bins):
                store_list.append(strided_snr(I, b, hop_length=hopper))
            
            # pandas dataframe (not nessesary)
            df_list = pd.DataFrame(store_list, index=bins, columns=np.round(wl,2))
            # ploting heatmap
            plt.subplot(221)
            plt.plot(wl, I, label="spectra")
            plt.ylabel("Flux")
            plt.xlabel("Wavelength")
            plt.title("Spectra of {0}".format(fname))
            plt.subplot(212)
            sns.set()
            cmap = sns.diverging_palette(220, 10, as_cmap=True)
            #ax = sns.heatmap(store_list, cmap=cmap, xticklabels=200, vmax=300, vmin=10)
            ax = sns.heatmap(df_list, cmap=cmap, xticklabels=200, vmax=300, vmin=0)
            #ax = sns.heatmap(df_list)
            #plt.xticks(np.arange(int(np.min(wl)), int(np.max(wl)+1), 1.0))
            ax.set(ylabel="Binsize", xlabel="Wavelength")
            plt.title("SNR")
            
            plt.show()

       # fancy colour plot if snr for different bin sizes verse central wavelenght value
       #
       #  b   |
       #  i   |
       #  n   |
       #      |    SNR Ratio
       #  s   |
       #  i   |
       #  z   |_________________
       #  e         Position
       #
       #

#test = [1,2,3,4,2,1,2,1,1,2,3,2,1,4,2,1,2,3,4,1,2,3,4,1,2,4,2,1,2,3,2,2.5,2,2,2,1,1,2,3,4]
#tests = [t / 100 for t in test]

# old code
    #x = range(len(test))

    #xbin = [xi + binsize/2 for xi in x]  # center snr values
    #xbin = xbin[:-binsize]
    
    #snrlist = snrscan(test, binsize)

    #plt.figure
    #plt.plot(x, test, label="data")
    #plt.plot(xbin, snrlist, label="snr sze"+ str(binsize))
    #plt.legend()
    #plt.show()

def fast_wav_selector(wav, flux, wav_min, wav_max, verbose=False):
    """ Faster Wavelenght selector
    
    If passed lists it will return lists.
    If passed np arrays it will return arrays
    
    Fastest is using np.ndarrays
    fast_wav_selector ~1000-2000 * quicker than wav_selector
    """

    if isinstance(wav, list): # if passed lists
        wav_sel = [wav_val for wav_val in wav if (wav_min < wav_val < wav_max)]
        flux_sel = [flux_val for wav_val, flux_val in zip(wav, flux) if (wav_min < wav_val < wav_max)]
    elif isinstance(wav, np.ndarray):
        # Super Fast masking with numpy
        mask = (wav > wav_min) & (wav < wav_max)
        if verbose:
            print("mask=", mask)
            print("len(mask)", len(mask))
            print("wav", wav)
            print("flux", flux)
        wav_sel = wav[mask]
        flux_sel = flux[mask]
    else:
          raise TypeError("Unsupported input wav type")
    return [wav_sel, flux_sel]


if __name__ == "__main__":
    args = vars(_parser())
    fname = args.pop('fname')
    opts = {k: args[k] for k in args}

    main(fname, **opts)
