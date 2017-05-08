#!/home/jneal/anaconda3/bin/python
"""Select optimal/non-optimal spectra in regards to specification
   Median and norm Combine spectra after
"""
from __future__ import division, print_function
import os
import sys
import logging
import argparse
import numpy as np
from tqdm import tqdm
from typing import Any
from astropy.io import fits
import matplotlib.pyplot as plt

from Get_filenames import get_filenames

# TODO: argparse file with suitbale nods. (optimal nod mask)
# TODO: optimal nod mask)

nod_mask_file = "optimal_nod_masks.txt"


# Try parse 4*8 bools.
def parse_boolgrid(filename: str, nod: int=8, chip: int=4) -> Any:
    """Parse file with 4*8 bool values."""
    line_num = 0
    boolgrid = np.empty((chip, nod), dtype=bool)
    with open(filename) as f:
        for line in f:
            if line.startswith("#"):
                pass
            else:
                line = line.strip()

                assert len(line) == nod, "line length is wrong={}".format(len(line))
                for nod_num, val in enumerate(line):
                    boolgrid[line_num, nod_num] = bool(int(val))
                line_num += 1
    return boolgrid
# DRACS Output quicklook/status.


# Plot all 8 reduced spectra
# Plot all 8 normalized reduced spectra
# plot mean combined and median combined spectra.
def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    # parser = GooeyParser(description='Remove : from data files')
    parser = argparse.ArgumentParser(description='Combines Nods using ')
    parser.add_argument('listspectra', help='List of spectra to combine.', default=False)
    parser.add_argument('-o', "--optimal-nods", help="Optimal nod bool matix file.")
    parser.add_argument("-s", "--spectralcoords", default=False, action="store_true",
                        help="Turn spectra into spectral corrdinates first before adding. Default=False")
    parser.add_argument("-n", "--nod_num", help="Number of nods in the nod cycle, default=8", default=8, type=int)
    parser.add_argument("-c", "--combination", help="Nod combination method, default=all means do all three.",
                        default="all", choices=["all", "optimal", "non-opt", "mix"])
    parser.add_argument("--snr", help="Show snr of continuum.", action="store_true")
    args = parser.parse_args()
    return args


def main(**kwargs):
    """Main function."""
    pass
            if kwargs["snr"]:
                # Analysis signal to noise in a part of the continuim of each spectra.
                # Normlazied result.
                print("For chip {}".format(chip_num))
                print("opt_norm_mean snr = {}".format(sampled_snr(opt_norm_mean, chip_num)))
                print("opt_norm_median snr = {}".format(sampled_snr(opt_norm_median, chip_num)))
                print("opt_norm_sum snr = {}".format(sampled_snr(opt_norm_sum, chip_num)))
                print("nonopt_norm_mean snr = {}".format(sampled_snr(nonopt_norm_mean, chip_num)))
                print("nonopt_norm_median snr = {}".format(sampled_snr(nonopt_median, chip_num)))
                print("nonopt_norm_sum snr = {}".format(sampled_snr(nonopt_norm_sum, chip_num)))
                print("mix_norm_mean snr = {}".format(sampled_snr(mix_norm_mean, chip_num)))
                print("mix_norm_median snr = {}".format(sampled_snr(mix_norm_median, chip_num)))
                print("mix_norm_sum snr = {}".format(sampled_snr(mix_norm_sum, chip_num)))


    return 0


def sampled_snr(spectrum, chip):
    """Sample SNR with Predefined continuim locations per chip."""
    limits = {1: [900, 960], 2: [460, 600], 3: [240, 310], 4: [450, 490]}
    section = spectrum[slice(limits[chip][0], limits[chip][1])]
    return np.mean(section) / np.std(section)


# TODO Add a level of Iteration!
# np.nan the bad pixels the loop untill doesnt change maybe?
# or just do 2-3 loops

def sigma_detect(nods, plot=True):
    """Detect the local pixels that are outside 4sigma from all nods.

    Local means 2 pixels either side.
    """
    if isinstance(nods, list):
        raise TypeError("Input an nod*pixel array please.")
    sig_clip = 4  # Sigma clipping Value.

    nods = np.asarray(nods)
    if nods.shape[0] > 8:
        raise ValueError("Too many nods (>8), check dimensions of input. ([nod, pixel])")

    print(nods.shape)
    bad_pixel_count = 0
    bad_pixel_record = []
    for pixel in range(nods.shape[1]):
        if (pixel < 2):
            near_pixels = nods[:, :5]      # First 5 pixels to do the 2 end pixels..
            grid_index = pixel
        elif pixel > (nods.shape[1] - 3):
            near_pixels = nods[:, -5:]     # Last 5 pixels to do the last 2 end pixels..
            grid_index = pixel - nods.shape[1]
        else:
            near_pixels = nods[:, slice(pixel - 2, pixel + 3)]
            grid_index = 2

        # ravel pixels near this picel output
        ravel_pixels = near_pixels.ravel()
        median = np.median(ravel_pixels)
        std = np.std(ravel_pixels)

        # The values of this pixel for all nods, taken from near_pixels. Should be same as nods[:, pixel]
        this_pixel = near_pixels[:, grid_index]
        assert np.all(this_pixel == nods[:, pixel])

        sig_over = this_pixel > (median + sig_clip * std)
        sig_below = this_pixel < (median - sig_clip * std)
        sig = sig_over | sig_below
        if np.any(sig):
            bad_pixel_count += 1
            print("A sigma clip value was found")
            print("pixel num", pixel)
            print("median", median, "std", std)
            bad_nod = sig.nonzero()
            if len(bad_nod) > 1:
                raise ValueError("More then one bad pixel this time")
            print(np.asarray(sig_over, dtype=int) - np.asarray(sig_below, dtype=int))
            # print("location of bad pixel above should be in position {}".format(bad_nod))
            # print(near_pixels)

            print("bad nod", bad_nod)
            bad_pixel_record += [[bad_nod[0][0], pixel, this_pixel[bad_nod[0][0]]]]
        # print(near_pixels)
    print("bad_pixel_count", bad_pixel_count)
    print("bad_pixel_record", bad_pixel_record)


    if plot:
        for i, nod in enumerate(nods):
            plt.plot(nod, label="nod {}".format(i))
        bad_pixel_x = [x[1] for x in bad_pixel_record]
        bad_pixel_y = [x[2] for x in bad_pixel_record]

        plt.plot(bad_pixel_x, bad_pixel_y, "x", label=">4 sigma")
        plt.legend()
        plt.title("Identifying bad pixels")
        plt.xlabel("pixel")
        plt.ylabel("norm flux")
        plt.show()

    return [pixel[0:2] for pixel in bad_pixel_record]





if __name__ == "__main__":
    args = vars(_parser())

    opts = {k: args[k] for k in args}
    sys.exit(main(**opts))

    #
