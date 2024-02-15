#!/usr/bin/env python
""" Plot the Iraf sum and the mixed combination spectra together and the differences"""

from __future__ import division, print_function

import argparse
import os
import sys

import matplotlib.pyplot as plt
from astropy.io import fits
from tqdm import tqdm

from octotribble.Get_filenames import get_filenames


def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description="Make Dracs quicklook plots.")
    parser.add_argument(
        "-p", "--showplots", help="Show the plots.", default=False, action="store_true"
    )
    parser.add_argument(
        "-c",
        "--chip",
        help="Specify chip to plot. Default is all.",
        choices=["1", "2", "3", "4"],
        default=None,
    )
    return parser.parse_args()


def main(chip=None, showplots=False):
    chips = range(1, 5) if chip is None else [int(chip)]
    dir_path = os.getcwd()

    intermediate_path = os.path.join(dir_path, "Intermediate_steps", "")
    combined_path = os.path.join(dir_path, "Combined_Nods", "")
    image_path = os.path.join(dir_path, "images", "")
    print(intermediate_path)
    observation_name = os.path.split(dir_path)[-1]

    for chip_num in tqdm(chips):
        print("Starting Chip # {}".format(chip_num))
        irafsum_name = get_filenames(
            combined_path, "CRIRE*norm.sum.fits", "*_{0}.*".format(chip_num)
        )
        mix_name = get_filenames(
            combined_path, "CRIRE*.mixavg.fits", "*_{0}.*".format(chip_num)
        )
        print(mix_name)
        optavg_name = get_filenames(
            combined_path, "CRIRE*.ms.norm.optavg.fits", "*_{0}.*".format(chip_num)
        )
        nonoptavg_name = get_filenames(
            combined_path, "CRIRE*.ms.norm.nonoptavg.fits", "*_{0}.*".format(chip_num)
        )

        irafsum_data = fits.getdata(combined_path + irafsum_name[0])
        iraf_opt = irafsum_data[0][0]
        iraf_nonopt = irafsum_data[1][0]
        mixavg_data = fits.getdata(combined_path + mix_name[0])
        optavg_data = fits.getdata(combined_path + optavg_name[0])
        nonoptavg_data = fits.getdata(combined_path + nonoptavg_name[0])

        fig = plt.figure()
        plt.plot(iraf_opt + 0, label="iraf optimal")
        plt.plot(iraf_nonopt + 0.04, label="iraf nonopt")
        plt.plot(optavg_data + 0.02, label="optavg")
        plt.legend()
        plt.xlabel("Pixel")
        plt.ylabel("Flux")
        plt.savefig(
            image_path + "opt-nonop-mix-{0}_{1}.png".format(observation_name, chip_num)
        )
        if showplots:
            plt.show()
        plt.close(fig)

        # mix_and_optimal
        fig = plt.figure()
        plt.plot(iraf_opt, label="iraf optimal")
        plt.plot(mixavg_data, "--", label="mixed")
        plt.xlabel("Pixel")
        plt.ylabel("Flux")
        plt.savefig(
            image_path
            + "mix_and_optimal-{0}_{1}.png".format(observation_name, chip_num)
        )
        plt.legend()
        if showplots:
            plt.show()
        plt.close(fig)

        # "Python - Iraf optimal combination"
        fig = plt.figure()
        plt.plot(optavg_data - iraf_opt)
        plt.xlabel("Pixel")
        plt.ylabel("Flux difference")
        plt.title("Python - Iraf optimal combination")
        plt.legend()
        plt.savefig(
            image_path
            + "optavg-optimalsum-{0}_{1}.png".format(observation_name, chip_num)
        )
        if showplots:
            plt.show()
        plt.close(fig)

        # Python - Iraf non-optimal combination
        fig = plt.figure()
        plt.plot(nonoptavg_data - iraf_nonopt)
        plt.xlabel("Pixel")
        plt.ylabel("Flux difference")
        plt.title("Python - Iraf non-optimal combination")
        plt.legend()
        plt.savefig(
            image_path
            + "nonoptavg-nonoptimalsum-{0}_{1}.png".format(observation_name, chip_num)
        )
        if showplots:
            plt.show()
        plt.close(fig)

        fig = plt.figure()
        plt.plot(mixavg_data - iraf_opt)
        plt.title("Mixed - Optimal")
        plt.xlabel("Pixel")
        plt.ylabel("Flux difference")
        plt.legend()
        plt.savefig(
            image_path + "mix-optimal-{0}_{1}.png".format(observation_name, chip_num)
        )
        if showplots:
            plt.show()
        plt.close(fig)

    return 0


if __name__ == "__main__":
    args = vars(_parser())
    opts = {k: args[k] for k in args}
    sys.exit(main(**opts))
