# DRACS Output quicklook/status.

# Plot all 8 reduced spectra
# Plot all 8 normalized reduced spectra
# plot mean combined and median combined spectra.

from __future__ import division, print_function
import os
import sys
import argparse
import numpy as np
from tqdm import tqdm
from astropy.io import fits
import matplotlib.pyplot as plt
from Get_filenames import get_filenames

# observation_name = "HD162020-1"

# path = "C:/Users/Jason/Dropbox/PhD/CriresReduction/{0}/".format(observation_name)

# intermediate_path =  path + "Intermediate_steps/"
# combined_path =  path + "Combined_Nods/"
# image_path = path + "quicklooks/"


def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Make Dracs quicklook plots.')
    parser.add_argument("-p", "--showplots", help="Show the plots.", default=False, action="store_true")
    parser.add_argument("-c", "--chip", help="Specify chip to plot. Default is all.", choices=["1", "2", "3", "4"], default=None)
    parser.add_argument("-b", "--band", help="Specify band to plot. Default is all.", choices=["1", "2", "3"], default=None)
    args = parser.parse_args()
    return args


def main(chip=None, band=None, showplots=False):

    if chip is None:
        chip = range(1, 5)
    else:
        chip = [int(chip)]

    if band is None:
        band = range(3)
    else:
        chip = [int(band)]
    # To run
    dir_path = os.getcwd()

    intermediate_path = dir_path + "/Intermediate_steps/"
    combined_path = dir_path + "/Combined_Nods/"
    image_path = dir_path + "/images/"
    observation_name = os.path.split(dir_path)[-1]

    for chip_num in tqdm(chip):
        print("Starting Chip # {}".format(chip_num))
        combined_name = get_filenames(combined_path, 'CRIRE*norm.sum.fits', "*_{0}.*".format(chip_num))
        nod_names = get_filenames(intermediate_path, 'CRIRE*.ms.fits', "*_{0}.*".format(chip_num))
        norm_names = get_filenames(intermediate_path, 'CRIRE*.ms.norm.fits', "*_{0}.*".format(chip_num))

        combined_data = fits.getdata(combined_path + combined_name[0])
        print("length of combined_data =", len(combined_data))
        if len(combined_data) == 3:
            optimal_median = []
            optimal_mean = []
            non_optimal_median = []
            non_optimal_mean = []

            for indx in band:
                print("index of extras ", indx)
                nod_data = [fits.getdata(name)[indx, 0] for name in nod_names]

                norm_data = [fits.getdata(name)[indx, 0] for name in norm_names]
                median_nod = np.median(norm_data, axis=0)    # Median combine normalzied spectra

                mean_nod = combined_data[indx][0]   # From norm.sum.fits file.

                mean_median_diff = mean_nod - median_nod   # For plot 4
                diff_mask = np.abs(mean_median_diff) > 0.025
                mean_median_diff[diff_mask] = np.nan

                mean_mask = (mean_nod > 1.15) | (mean_nod < 0.0)
                mean_nod[mean_mask] = np.nan
                median_mask = (median_nod > 1.15) | (median_nod < 0.0)
                median_nod[median_mask] = np.nan

                # Plot Results
                fig = plt.figure(figsize=(10, 10))
                plt.suptitle("{0}, Chip-{1}".format(observation_name, chip_num), fontsize=16)
                ax1 = plt.subplot(411)
                for i, data in enumerate(nod_data):
                    # add some limits to data displayed, help the scaling.
                    data_mask = (data > 4 * np.median(data)) | (data < 0.0)
                    data[data_mask] = np.nan
                    ax1.plot(data, label=i + 1)
                plt.ylabel("Intensity")
                plt.title("Extracted Nod Spectra")
                plt.xlim([0, 1024])
                # start, end = ax1.get_ylim()
                # ax1.yaxis.set_ticks(np.arange(start, end, 2000))
                plt.tight_layout(pad=2.5, w_pad=0.0, h_pad=0.5)
                ax1.legend()
                # del nod_data

                print("starting plot 2")
                ax2 = plt.subplot(412)
                for data in norm_data:
                    data_mask = (data > 4 * 1.2) | (data < 0.0)
                    data[data_mask] = np.nan
                    ax2.plot(data)
                plt.ylabel("Normalized\nIntensity")
                plt.title("Normalized Nod Spectra")
                plt.xlim([0, 1024])
                start, end = ax2.get_ylim()
                # ax2.yaxis.set_ticks(np.arange(start, end, 0.1))

                print("starting plot 3")
                ax3 = plt.subplot(413)
                ax3.plot(mean_nod, label="Nod Mean")
                ax3.plot(median_nod, '--r', label="Nod Median")
                # plt.xlabel("Pixel Position", fontsize=14)
                plt.ylabel("Normalized\nIntensity")
                plt.title("Combined Nod Spectra")
                ax3.legend(loc=0)
                plt.xlim([0, 1024])
                # start, end = ax3.get_ylim()
                # ax3.yaxis.set_ticks(np.arange(start, end, 0.1))
                plt.tight_layout(pad=2.5, w_pad=0.0, h_pad=0.5)
                # plt.show()

                print("starting plot 4")
                mean_median_diff = mean_nod - median_nod
                mean_median_diff[np.abs(mean_median_diff) > 0.02] = np.nan   # mask out the large differences.
                ax4 = plt.subplot(414)
                ax4.plot(mean_median_diff, label="Mean-median")
                plt.xlabel("Pixel Position", fontsize=10)
                plt.ylabel("Flux diff")
                plt.title("Mean-median difference.")
                ax4.legend(loc=0)
                plt.xlim([0, 1024])
                start, end = ax4.get_ylim()
                # ax4.yaxis.set_ticks(np.arange(start, end, 0.01))
                plt.tight_layout(pad=2.5, w_pad=0.0, h_pad=0.5)
                if showplots:
                    plt.show()
                print("Saving ...")
                # Save figure
                # fig.savefig(image_path + "quicklook_{0}_{1}_reduction_band{2}.pdf".format(observation_name, chip_num,
                #            indx + 1))
                fig.savefig(image_path + "quicklook_{0}_{1}_reduction_band{2}.png".format(observation_name, chip_num,
                            indx + 1))
                plt.close(fig)

                if indx == 0:
                    optimal_median = median_nod
                    optimal_mean = mean_nod
                    print(np.max(optimal_mean))
                    print(np.max(optimal_median))
                elif indx == 1:
                    non_optimal_median = median_nod
                    non_optimal_mean = mean_nod
                    print(np.max(non_optimal_mean))
                    print(np.max(non_optimal_median))

            diff_of_mean = optimal_mean - non_optimal_mean
            diff_of_median = optimal_median - non_optimal_median


            print("diff_of_mean", diff_of_mean)

            # After looping though orders plot difference between mean and median spectra
            fig = plt.figure(figsize=(10, 8))
            plt.suptitle("{0}, Chip-{1}".format(observation_name, chip_num), fontsize=16)
            ax1 = plt.subplot(111)
            plt.plot(diff_of_mean, label="mean")
            plt.plot(diff_of_median, "--", label="medain")
            plt.ylim([-0.02, 0.02])    # Limit to +- 2%
            plt.title("Differences between optimal - non-optimal combined spectra.")
            plt.ylabel("Flux diff")
            plt.legend(loc=0)
            fig.savefig(image_path + "combine_diff_{0}_{1}_reduction_opt_minus_nonopt.png".format(observation_name,
                                                                                                  chip_num))
            if showplots:
                plt.show()
            plt.close(fig)

        else:
            nod_data = [fits.getdata(name) for name in nod_names]
            norm_data = [fits.getdata(name) for name in norm_names]
            median_nod = np.median(norm_data, axis=0)    # Median combine normalzied spectra

            # Plot Reuslts
            fig = plt.figure()
            plt.suptitle("{0}, Chip-{1}".format(observation_name, chip_num), fontsize=16)
            ax1 = plt.subplot(311)
            for i, data in enumerate(nod_data):
                ax1.plot(data, label=i + 1)
            plt.ylabel("Intensity")
            plt.title("Extracted Nod Spectra")
            plt.xlim([0, 1024])
            # start, end = ax1.get_ylim()
            # ax1.yaxis.set_ticks(np.arange(start, end, 2000))
            # ax1.legend()

            ax2 = plt.subplot(312)
            for data in norm_data:
                ax2.plot(data)
            plt.ylabel("Normalized\nIntensity")
            plt.title("Normalized Nod Spectra")
            plt.xlim([0, 1024])
            start, end = ax2.get_ylim()
            ax2.yaxis.set_ticks(np.arange(start, end, 0.1))

            ax3 = plt.subplot(313)
            ax3.plot(combined_data, label="Nod Mean")
            ax3.plot(median_nod, '--r', label="Nod Median")
            plt.xlabel("Pixel Position", fontsize=14)
            plt.ylabel("Normalized\nIntensity")
            plt.title("Combined Nod Spectra")
            plt.legend(loc=0)
            plt.xlim([0, 1024])
            start, end = ax3.get_ylim()
            ax3.yaxis.set_ticks(np.arange(start, end, 0.1))
            plt.tight_layout(pad=2.5, w_pad=0.0, h_pad=0.5)
            if showplots:
                plt.show()

            # Save figure
            fig.savefig(image_path + "quicklook_{0}_{1}_reduction.pdf".format(observation_name, chip_num))
            fig.savefig(image_path + "quicklook_{0}_{1}_reduction.png".format(observation_name, chip_num))
            plt.close(fig)

    return 0


if __name__ == '__main__':
    args = vars(_parser())
    opts = {k: args[k] for k in args}
    sys.exit(main(**opts))
