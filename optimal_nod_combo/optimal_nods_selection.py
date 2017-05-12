#!/home/jneal/anaconda3/bin/python

"""Select optimal/non-optimal spectra in regards to specification.

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

import bp_replacement as bp
from Get_filenames import get_filenames

# TODO: argparse file with suitbale nods. (optimal nod mask)
# TODO: optimal nod mask)

nod_mask_file = "optimal_nod_masks.txt"


# Try parse 4*8 bools.
def parse_boolgrid(filename, nod=8, chip=4):
    # type: (str, int, int) -> np.ndarray
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
    # type: () -> argparse.Namespace
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
    # type: (...) -> bool
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

            # ax2 = plt.subplot(212)
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

            # turn to array here ( can do before just need to work out how/where what to change)
            fix_opt_nods = clean_nods(optimal_norm_nods)
            fix_nonopt_nods = clean_nods(nonoptimal_norm_nods)
            fix_mix_nods = clean_nods(mix_norm_nods)

            # bad_pixel_record = bp.sigma_detect(mix_norm_nods_arr, plot=True)
            #
            # bp.warn_consec_badpixels(bad_pixel_record, stop=True)
            #
            # fix_mix_nods = bp.interp_badpixels(mix_norm_nods_arr, bad_pixel_record)
            # if len(bad_pixel_record) > 0:
            #     assert np.any(fix_mix_nods != mix_norm_nods_arr)

            plt.plot(fix_opt_nods.T, ".", label="opt")
            plt.plot(fix_mix_nods.T, label="fixed")
            plt.title("After bad pixel correction")
            plt.legend()
            plt.show()

    return 0


def clean_nods(nods):
    # type: (Union[np.ndarray, List[List[float]]]) -> np.ndarray
    """Clean bad pixels from the nods."""
    nod_array = np.asarray(nods)
    bad_pixels = bp.sigma_detect(nod_array, plot=True)

    bp.warn_consec_badpixels(bad_pixels, stop=True)

    fixed_nods = bp.interp_badpixels(nod_array, bad_pixels)
    print("Number of bad pixels = {}".format(len(bad_pixels)))
    if len(bad_pixels) > 0:
        assert np.any(fixed_nods != nod_array)

    return fixed_nods


def sampled_snr(spectrum, chip):
    # type: (Any, int) -> np.float64
    """Sample SNR with Predefined continuim locations per chip."""
    limits = {1: [900, 960], 2: [460, 600], 3: [240, 310], 4: [450, 490]}
    section = spectrum[slice(limits[chip][0], limits[chip][1])]
    return np.mean(section) / np.std(section)


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
