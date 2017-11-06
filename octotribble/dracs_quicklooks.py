# DRACS Output quicklook/status.

# Plot all 8 reduced spectra
# Plot all 8 normalized reduced spectra
# plot mean combined and median combined spectra.

from __future__ import division, print_function

import os

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

from Get_filenames import get_filenames

# observation_name = "HD162020-1"

# path = "C:/Users/Jason/Dropbox/PhD/CriresReduction/{0}/".format(observation_name)

# intermediate_path =  path + "Intermediate_steps/"
# combined_path =  path + "Combined_Nods/"
# image_path = path + "quicklooks/"

# To run
dir_path = os.getcwd()

intermediate_path = dir_path + "/Intermediate_steps/"
combined_path = dir_path + "/Combined_Nods/"
image_path = dir_path + "/images/"
observation_name = os.path.split(dir_path)[-1]

for chip_num in range(1, 5):
    combined_name = get_filenames(combined_path, 'CRIRE*.sum.fits', "*_{0}.*".format(chip_num))
    nod_names = get_filenames(intermediate_path, 'CRIRE*.ms.fits', "*_{0}.*".format(chip_num))
    norm_names = get_filenames(intermediate_path, 'CRIRE*.ms.norm.fits', "*_{0}.*".format(chip_num))

    combined_data = fits.getdata(combined_path + combined_name[0])
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
    # plt.show()

    # Save figure
    fig.savefig(image_path + "quicklook_{0}_{1}_reduction.pdf".format(observation_name, chip_num))
    fig.savefig(image_path + "quicklook_{0}_{1}_reduction.png".format(observation_name, chip_num))
    plt.close(fig)
