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

    # Raising errors on non implemented features.
    if kwargs["spectralcoords"]:
        raise NotImplementedError("Spectral coordinate not available.")

    if kwargs["nod_num"] != 8:
        raise NotImplementedError("A nod number that is not 8 has not been implemented.")

    # List of combination methods to use.
    if kwargs["combination"] == "all":
        comb_methods = ["optimal", "non-opt", "mix"]
    else:
        comb_methods = kwargs["combination"]

    # Get optimal nod grid
    if kwargs["optimal_nods"]:
        try:
            nod_mask = parse_boolgrid(kwargs["optimal_nods"])
        except Exception as e:
            print("Issue with the parse_boolgrid")
            raise e
    else:
        if "mix" in comb_methods:
            raise ValueError("No optimal_nod file supplied for 'mix' combination.")

    # Now need to load in the spectra. Cycle over combination options []

    # def file structure
    dir_path = os.getcwd()
    intermediate_path = dir_path + "/Intermediate_steps/"
    combined_path = dir_path + "/Combined_Nods/"
    image_path = dir_path + "/images/"
    observation_name = os.path.split(dir_path)[-1]


    for chip_num in tqdm(range(1, 5)):
        # print("Starting Chip # {}".format(chip_num))
        combined_name = get_filenames(combined_path, 'CRIRE*norm.sum.fits', "*_{0}.*".format(chip_num))
        nod_names = get_filenames(intermediate_path, 'CRIRE*.ms.fits', "*_{0}.*".format(chip_num))
        norm_names = get_filenames(intermediate_path, 'CRIRE*.ms.norm.fits', "*_{0}.*".format(chip_num))

        # print(norm_names)

        combined_data = fits.getdata(combined_path + combined_name[0])
    # for combo in comb_methods:
        if combined_data.shape != (3, 1, 1024):

            raise ValueError("Fits files do no have the right shape, [3,1,1024]")
        else:
            optimal_nods = [fits.getdata(name)[0, 0] for name in nod_names]
            optimal_norm_nods = [fits.getdata(name)[0, 0] for name in norm_names]
            nonoptimal_nods = [fits.getdata(name)[1, 0] for name in nod_names]
            nonoptimal_norm_nods = [fits.getdata(name)[1, 0] for name in norm_names]

            # print("this chips nod bools =", nod_mask[chip_num - 1])
            if "mix" in comb_methods:
                # True values map to index zero. whereas False maps to index 1.
                band_index = [int(not x) for x in nod_mask[chip_num - 1]]
                mix_nods = [fits.getdata(name)[indx, 0] for name, indx in zip(nod_names, band_index)]
                mix_norm_nods = [fits.getdata(name)[indx, 0] for name, indx in zip(norm_names, band_index)]
            else:
                mix_nods, mix_norm_nods = [], []


            opt_median = np.median(optimal_nods, axis=0)
            opt_mean = np.mean(optimal_nods, axis=0)
            opt_norm_median = np.median(optimal_norm_nods, axis=0)
            opt_norm_mean = np.mean(optimal_norm_nods, axis=0)
            opt_norm_sum = np.sum(optimal_norm_nods, axis=0)
            nonopt_median = np.median(nonoptimal_nods, axis=0)
            nonopt_mean = np.mean(nonoptimal_nods, axis=0)
            nonopt_norm_median = np.median(nonoptimal_norm_nods, axis=0)
            nonopt_norm_mean = np.mean(nonoptimal_norm_nods, axis=0)
            nonopt_norm_sum = np.sum(nonoptimal_norm_nods, axis=0)
            mix_median = np.median(mix_nods, axis=0)
            mix_mean = np.mean(mix_nods, axis=0)
            mix_norm_median = np.median(mix_norm_nods, axis=0)
            mix_norm_mean = np.mean(mix_norm_nods, axis=0)
            mix_norm_sum = np.sum(mix_norm_nods, axis=0)

            fig = plt.figure()
            # ax1 = plt.subplot(211)
            # plt.plot(opt_median, label="opt_median")
            # plt.plot(opt_mean, label="opt_mean")
            # plt.plot(nonopt_median, label="nonopt_median")
            # plt.plot(nonopt_mean, label="nonopt_mean")
            # plt.plot(mix_median, label="mix_median")
            # plt.plot(mix_mean, label="mix_mean")
            # plt.legend()


            #ax2 = plt.subplot(212)
            plt.plot(opt_norm_median, label="opt_norm_median")
            plt.plot(opt_norm_mean, label="opt_norm_mean")
            plt.plot(nonopt_norm_median, label="nonopt_norm_median")
            plt.plot(nonopt_norm_mean, label="nonopt_norm_mean")
            plt.plot(mix_norm_median, label="mix_norm_median")
            plt.plot(mix_norm_mean, label="mix_norm_mean")

            plt.legend()
            plt.show()
            # print result

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
        # plot Results

        # save results
            # try sigma clipping on norm data

            # turn to array here ( can do before jsut need to work out how/where what to change)
            mix_norm_nods_arr = np.asarray(mix_norm_nods)

            bad_pixel_record = sigma_detect(mix_norm_nods_arr, plot="True")

            fix_mix_nods = interp_badpixels(mix_norm_nods_arr, bad_pixel_record)
            if len(bad_pixel_record) > 0:
                assert np.any(fix_mix_nods != mix_norm_nods_arr)
            plt.plot(mix_norm_nods_arr.T, ".")
            plt.plot(fix_mix_nods.T, label="fixed")
            plt.legend()
            plt.show()


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
    if nods.shape[0] > 8:
        raise ValueError("Too many nods (>8), check dimensions of input. ([nod, pixel])")
    # aviod mutation
    old_nods = np.empty_like(nods)
    old_nods[:] = np.nan   # set to nans to start iteration
    new_nods = np.empty_like(nods)
    new_nods[:] = nods

    bad_pixel_count = 0
    bad_pixel_record = []

    iteration = 0
    # Iterate untill no more bad pixels are replaced by nans during an iteration. or less than 5.
    while ((iteration < 5) & np.any(np.isnan(old_nods) != np.isnan(new_nods))):
        old_nods[:] = new_nods
        for pixel in range(new_nods.shape[1]):
            if (pixel < 2):
                near_pixels = new_nods[:, :5]      # First 5 pixels to do the 2 end pixels..
                grid_index = pixel
            elif pixel > (new_nods.shape[1] - 3):
                near_pixels = new_nods[:, -5:]     # Last 5 pixels to do the last 2 end pixels..
                grid_index = pixel - new_nods.shape[1]
            else:
                near_pixels = new_nods[:, slice(pixel - 2, pixel + 3)]
                grid_index = 2

            # ravel pixels near this picel output
            ravel_pixels = near_pixels.ravel()
            median = np.nanmedian(ravel_pixels)   # ignore nan values
            std = np.nanstd(ravel_pixels)         # ignore nan values

            # The values of this pixel for all nods, taken from near_pixels. Should be same as new_nods[:, pixel]
            this_pixel = near_pixels[:, grid_index]
            assert np.all((this_pixel == new_nods[:, pixel]) | np.isnan(this_pixel))

            # Find values outside the sig_clip level.
            sig = (this_pixel > (median + sig_clip * std)) | (this_pixel < (median - sig_clip * std))

            if np.any(sig):
                bad_nod = sig.nonzero()[0]
                if len(bad_nod) > 1:
                    logging.Warning("More then one nod has a bad pixel value here, in pixel #{}".format(pixel))
                for val in bad_nod:
                    bad_pixel_count += 1
                    bad_pixel_record += [(val, pixel, this_pixel[val])]
        for record in bad_pixel_record:
            new_nods[record[0], record[1]] = np.nan
        iteration += 1

    print("# Pixels outside {0}sigma = {1}".format(sig_clip, bad_pixel_count))
    # print("bad_pixel_record", bad_pixel_record)

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



def interp_badpixels(nods, bad_pixels):
    """Linearly interpolate over nearby pixels.

    If it is at the end then just replace with the next pixel value.

    Parameters:
    nods: array
        A nod*pixel array of the pixel flux values.
    bad_pixels: list of lists of ints
        Index position [nod, pixel] of bad pixels to replace.

    Returns
    -------
    nods: array
        Array of nods with bad pixels interpolated over.
    """
    if isinstance(nods, list):
        raise TypeError("Input an nod*pixel array please.")

    # Warn about consecutive bad_pixels
    warn_consec_badpixels(bad_pixels, stop=False)
    output = np.empty_like(nods)
    output[:] = nods            # aviods mutation

    for pixel in bad_pixels:
        if pixel[1] == 0:
            replacement = nods[pixel[0], pixel[1] + 1]   # pixel number 2
        elif pixel[1] == (nods.shape[1] - 1):
            replacement = nods[pixel[0], pixel[1] - 1]    # pixel number -2
        else:
            if (left_consec_search(pixel, bad_pixels) + right_consec_search(pixel, bad_pixels)) != 0:
                replacement = -1
                logging.warning("Need to finish up here")
            else: # No consecutives
                # Interpolation when have both side pixels.
                x = [0, 1, 2]
                xp = [0, 2]
                fp = [nods[pixel[0], pixel[1] - 1], nods[pixel[0], pixel[1] + 1]]
                y = np.interp(x, xp, fp)
                replacement = y[1]
                # print(x, xp, fp, y)

        output[pixel[0], pixel[1]] = replacement

    return output


def left_consec_search(pixel, bad_pixels):
    """Count number of consecutive bad pixels to the left of this pixel."""
    prev_pixel = [pixel[0], pixel[1] - 1]
    if (prev_pixel in bad_pixels) or (tuple(prev_pixel) in bad_pixels):
        # print("prev_pixel in recursion", prev_pixel)
        return left_consec_search(prev_pixel, bad_pixels) + 1
    else:
        return 0


def right_consec_search(pixel, bad_pixels):
    """Count number of consecutive bad pixels to the right of this."""
    next_pixel = [pixel[0], pixel[1] + 1]
    if (next_pixel in bad_pixels) or (tuple(next_pixel) in bad_pixels):
        # print("next_pixel in recursion", next_pixel)
        return right_consec_search(next_pixel, bad_pixels) + 1
    else:
        return 0


def consec_badpixels(bad_pixels):
    """Check for consecutive badpixels in the same nod.

    Consecutive in in axis=1 for the same axis=0 value.

    parameters
    ----------
    bad_pixels: list of list of ints
        List of index locations of bad pixels [nod, pixel].

    returns
    -------
    is_consec: bool
        True if a consecutve bad pixel found."""

    for pixel in bad_pixels:
        left_pix = [pixel[0], pixel[1] - 1]
        right_pix = [pixel[0], pixel[1] + 1]
        if (left_pix in bad_pixels) or (right_pix in bad_pixels):
            return True

    return False


def warn_consec_badpixels(bad_pixels, stop=True):
    """Raise erro on consecutive badpixels in the same nod.

    parameters
    ----------
    bad_pixels: list of list of ints
        List of index locations of bad pixels [nod, pixel]."""

    if consec_badpixels(bad_pixels):
        if stop:
            raise ValueError("Consective bad pixels were found. Need to deal with these.")
        else:
            logging.warning("Consective bad pixels in a nod were found.")
    return None


if __name__ == "__main__":
    args = vars(_parser())

    opts = {k: args[k] for k in args}
    sys.exit(main(**opts))

    #
    # import numpy as np
    # from tqdm import tqdm
    # from astropy.io import fits
    # import matplotlib.pyplot as plt
    # from Get_filenames import get_filenames
    #
    # # observation_name = "HD162020-1"
    #
    # # To run
    # dir_path = os.getcwd()
    #
    # intermediate_path = dir_path + "/Intermediate_steps/"
    # combined_path = dir_path + "/Combined_Nods/"
    # image_path = dir_path + "/images/"
    # observation_name = os.path.split(dir_path)[-1]
    #
    # for chip_num in tqdm(range(1, 5)):
    #     print("Starting Chip # {}".format(chip_num))
    #     combined_name = get_filenames(combined_path, 'CRIRE*norm.sum.fits', "*_{0}.*".format(chip_num))
    #     nod_names = get_filenames(intermediate_path, 'CRIRE*.ms.fits', "*_{0}.*".format(chip_num))
    #     norm_names = get_filenames(intermediate_path, 'CRIRE*.ms.norm.fits', "*_{0}.*".format(chip_num))
    #
    #     combined_data = fits.getdata(combined_path + combined_name[0])
    #     print("length of combined_data =", len(combined_data))
    #     if len(combined_data) == 3:
    #         optimal_median = []
    #         optimal_mean = []
    #         non_optimal_median = []
    #         non_optimal_mean = []
    #
    #         for indx in range(3):
    #             print("index of extras ", indx)
    #             nod_data = [fits.getdata(name)[indx, 0] for name in nod_names]
    #
    #             norm_data = [fits.getdata(name)[indx, 0] for name in norm_names]
    #             median_nod = np.median(norm_data, axis=0)    # Median combine normalzied spectra
    #             # print("nod_data", nod_data)
    #             # print("norm_data", norm_data)
    #             # print("median_nod", median_nod)
    #
    #             # Plot Results
    #             fig = plt.figure(figsize=(10, 10))
    #             plt.suptitle("{0}, Chip-{1}".format(observation_name, chip_num), fontsize=16)
    #             ax1 = plt.subplot(411)
    #             for i, data in enumerate(nod_data):
    #                 # print("i for plot 1", i)
    #                 # print("data", data)
    #                 ax1.plot(data, label=i + 1)
    #             plt.ylabel("Intensity")
    #             plt.title("Extracted Nod Spectra")
    #             plt.xlim([0, 1024])
    #             # start, end = ax1.get_ylim()
    #             # ax1.yaxis.set_ticks(np.arange(start, end, 2000))
    #             ax1.legend()
    #
    #             print("starting plot 2")
    #             ax2 = plt.subplot(412)
    #             for data in norm_data:
    #                 ax2.plot(data)
    #             plt.ylabel("Normalized\nIntensity")
    #             plt.title("Normalized Nod Spectra")
    #             plt.xlim([0, 1024])
    #             start, end = ax2.get_ylim()
    #             #ax2.yaxis.set_ticks(np.arange(start, end, 0.1))
    #
    #             print("starting plot 3")
    #             ax3 = plt.subplot(413)
    #             ax3.plot(combined_data[indx][0], label="Nod Mean")
    #             ax3.plot(median_nod, '--r', label="Nod Median")
    #             plt.xlabel("Pixel Position", fontsize=14)
    #             plt.ylabel("Normalized\nIntensity")
    #             plt.title("Combined Nod Spectra")
    #             ax3.legend(loc=0)
    #             plt.xlim([0, 1024])
    #             start, end = ax3.get_ylim()
    #             ax3.yaxis.set_ticks(np.arange(start, end, 0.1))
    #             plt.tight_layout(pad=2.5, w_pad=0.0, h_pad=0.5)
    #             # plt.show()
    #
    #             print("starting plot 4")
    #             ax4 = plt.subplot(414)
    #             ax4.plot(combined_data[indx][0] - median_nod, label="Mean-median")
    #             plt.xlabel("Pixel Position", fontsize=14)
    #             plt.ylabel(r"\Delta norm I")
    #             plt.title("Mean-median difference.")
    #             ax4.legend(loc=0)
    #             plt.xlim([0, 1024])
    #             start, end = ax4.get_ylim()
    #             ax4.yaxis.set_ticks(np.arange(start, end, 0.01))
    #             plt.tight_layout(pad=2.5, w_pad=0.0, h_pad=0.5)
    #             # plt.show()
    #             print("Saving ...")
    #             # Save figure
    #            # fig.savefig(image_path + "quicklook_{0}_{1}_reduction_band{2}.pdf".format(observation_name, chip_num,
    #             #            indx + 1))
    #             fig.savefig(image_path + "quicklook_{0}_{1}_reduction_band{2}.png".format(observation_name, chip_num,
    #                         indx + 1))
    #             plt.close(fig)
    #
    #
    #             if indx == 0:
    #                 optimal_median = median_nod
    #                 optimal_mean = combined_data[indx][0]
    #             elif indx == 1:
    #                 non_optimal_median = median_nod
    #                 non_optimal_mean = combined_data[indx][0]
    #
    #         # After looping though orders plot difference between mean and median spectra
    #         fig = plt.figure(figsize=(10, 8))
    #         ax1 = plt.subplot(111)
    #         plt.plot(optimal_mean - non_optimal_mean, label="mean")
    #         plt.plot(optimal_median - non_optimal_median, "--", label="medain")
    #         plt.title("Differences between optimal - non-optimal combined spectra.")
    #         plt.ylabel("\Delta norm I")
    #         plt.legend(loc=0)
    #         fig.savefig(image_path + "combine_diff_{0}_{1}_reduction_band{2}.png".format(observation_name, chip_num,
    #                     indx + 1))
    #         plt.close(fig)
    #
    #     else:
    #         nod_data = [fits.getdata(name) for name in nod_names]
    #         norm_data = [fits.getdata(name) for name in norm_names]
    #         median_nod = np.median(norm_data, axis=0)    # Median combine normalzied spectra
    #
    #         # Plot Reuslts
    #         fig = plt.figure()
    #         plt.suptitle("{0}, Chip-{1}".format(observation_name, chip_num), fontsize=16)
    #         ax1 = plt.subplot(311)
    #         for i, data in enumerate(nod_data):
    #             ax1.plot(data, label=i + 1)
    #         plt.ylabel("Intensity")
    #         plt.title("Extracted Nod Spectra")
    #         plt.xlim([0, 1024])
    #         # start, end = ax1.get_ylim()
    #         # ax1.yaxis.set_ticks(np.arange(start, end, 2000))
    #         # ax1.legend()
    #
    #         ax2 = plt.subplot(312)
    #         for data in norm_data:
    #             ax2.plot(data)
    #         plt.ylabel("Normalized\nIntensity")
    #         plt.title("Normalized Nod Spectra")
    #         plt.xlim([0, 1024])
    #         start, end = ax2.get_ylim()
    #         ax2.yaxis.set_ticks(np.arange(start, end, 0.1))
    #
    #         ax3 = plt.subplot(313)
    #         ax3.plot(combined_data, label="Nod Mean")
    #         ax3.plot(median_nod, '--r', label="Nod Median")
    #         plt.xlabel("Pixel Position", fontsize=14)
    #         plt.ylabel("Normalized\nIntensity")
    #         plt.title("Combined Nod Spectra")
    #         plt.legend(loc=0)
    #         plt.xlim([0, 1024])
    #         start, end = ax3.get_ylim()
    #         ax3.yaxis.set_ticks(np.arange(start, end, 0.1))
    #         plt.tight_layout(pad=2.5, w_pad=0.0, h_pad=0.5)
    #         # plt.show()
    #
    #         # Save figure
    #         fig.savefig(image_path + "quicklook_{0}_{1}_reduction.pdf".format(observation_name, chip_num))
    #         fig.savefig(image_path + "quicklook_{0}_{1}_reduction.png".format(observation_name, chip_num))
    #         plt.close(fig)