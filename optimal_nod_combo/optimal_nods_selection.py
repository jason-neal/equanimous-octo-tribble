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
    parser.add_argument("-u", "--unnorm", help="Combine the unnormalized nods.", action="store_true")
    parser.add_argument("--snr", help="Show snr of continuum.", action="store_true")
    parser.add_argument("-p", "--plot", help="Show the plots.", action="store_true")
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
    print("Combiunation mehtod", comb_methods)

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
        else:
            nod_mask = None

    # def file structure
    dir_path = os.getcwd()
    intermediate_path = os.path.join(dir_path, "Intermediate_steps", "")
    combined_path = os.path.join(dir_path, "Combined_Nods", "")
    image_path = os.path.join(dir_path, "images", "")
    observation_name = os.path.split(dir_path)[-1]

    for chip_num in tqdm(range(1, 5)):
        combined_name = get_filenames(combined_path, 'CRIRE*norm.sum.fits', "*_{0}.*".format(chip_num))

        if kwargs["snr"]:
            snr_nod_names = get_filenames(intermediate_path, 'CRIRE*.ms.fits', "*_{0}.*".format(chip_num))
            snr_norm_names = get_filenames(intermediate_path, 'CRIRE*.ms.norm.fits', "*_{0}.*".format(chip_num))
            snr_calculations(snr_nod_names, snr_norm_names, nod_mask, chip_num)

        if kwargs["unnorm"]:
            nod_names = get_filenames(intermediate_path, 'CRIRE*.ms.fits', "*_{0}.*".format(chip_num))
        else:
            nod_names = get_filenames(intermediate_path, 'CRIRE*.ms.norm.fits', "*_{0}.*".format(chip_num))

        combined_data = fits.getdata(os.path.join(combined_path, combined_name[0]))
        for combo in comb_methods:
            print("This combo = {}".format(combo))
            if combined_data.shape != (3, 1, 1024):
                raise ValueError("Fits files do no have the right shape, [3, 1, 1024]")

            if combo == "optimal":  # "optimal", "non-opt", "mix"
                nods = [fits.getdata(name)[0, 0] for name in nod_names]
                nod_combo_name = "Optimal"
            elif combo == "non-opt":
                nods = [fits.getdata(name)[1, 0] for name in nod_names]
                nod_combo_name = "Non-optimal"
            elif combo == "mix":
                # True values map to index zero. whereas False maps to index 1.
                band_index = [int(not x) for x in nod_mask[chip_num - 1]]
                nods = [fits.getdata(name)[indx, 0] for name, indx in zip(nod_names, band_index)]
                nod_combo_name = "Mixed"

            # Replace bad pixels in normalized spectra
            pbfix_nods, bad_pixels = clean_nods(nods)

            mean_nods = np.mean(nods, axis=0)

            mean_pbfix_nods = np.mean(pbfix_nods, axis=0)

            # Plot Results
            if kwargs["plot"]:
                plt.figure()
                plt.subplot(211)
                plt.plot(mean_nods, label="{}".format(nod_combo_name))
                plt.plot(mean_pbfix_nods, label="Fixed {}".format(nod_combo_name))
                plt.ylabel("Flux")
                plt.legend()
                if kwargs["unnorm"]:
                    plt.title("Combined Spectra")
                else:
                    plt.title("Combined Normalized Spectra.")

                # Add difference plot
                plt.subplot(212)
                if kwargs["unnorm"]:
                    original_combine = combined_data[1, 0, :]
                else:
                    original_combine = combined_data[0, 0, :]
                plt.plot(mean_pbfix_nods - original_combine, label="{} - Optimal".format(nod_combo_name))
                plt.title("Difference from optimal combination.")
                plt.ylabel("Flux")
                plt.xlabel("Pixel")

                plt.show()
            else:
                pass

        # save results
        if combo == "optimal":  # "optimal", "non-opt", "mix"
            Output_name = "..... .norm.optavg.fits"
        elif combo == "non-opt":
            Output_name = "..... .norm.nonoptavg.fits"
        elif combo == "mix":
            Output_name = "..... .norm.mixavg.fits"

        # Need to do the save to fits file stuff... See old codes.



    return 0


def nod_calcs(nods):
    # type: (np.ndarray) -> List[np.ndarray, np.ndarray, np.ndarray]
    """Returns calculation across the nods."""
    nod_mean = np.mean(nods, axis=0)
    nod_median = np.median(nods, axis=0)
    nod_sum = np.sum(nods, axis=0)
    return nod_mean, nod_median, nod_sum


def clean_nods(nods):
    # type: (Union[np.ndarray, List[List[float]]]) -> np.ndarray
    """Clean bad pixels from the nods."""
    nod_array = np.asarray(nods)
    bad_pixels = bp.sigma_detect(nod_array, plot=False)

    bp.warn_consec_badpixels(bad_pixels, stop=False)

    fixed_nods = bp.interp_badpixels(nod_array, bad_pixels)
    print("Number of bad pixels = {}".format(len(bad_pixels)))
    if len(bad_pixels) > 0:
        assert np.any(fixed_nods != nod_array)

    return fixed_nods, bad_pixels


def snr_calculations(nod_names, norm_names, nod_mask, chip_num):
    # Analysis signal to noise in a part of the continuim of each spectra.
    optimal_nods = [fits.getdata(name)[0, 0] for name in nod_names]
    optimal_norm_nods = [fits.getdata(name)[0, 0] for name in norm_names]
    nonoptimal_nods = [fits.getdata(name)[1, 0] for name in nod_names]
    nonoptimal_norm_nods = [fits.getdata(name)[1, 0] for name in norm_names]

    if nod_mask is not None:
        band_index = [int(not x) for x in nod_mask[chip_num - 1]]
        mix_nods = [fits.getdata(name)[indx, 0] for name, indx in zip(nod_names, band_index)]
        mix_norm_nods = [fits.getdata(name)[indx, 0] for name, indx in zip(norm_names, band_index)]
    else:
        mix_nods = []
        mix_norm_nods = []

    d = {"optimal_nods": optimal_nods, "optimal_norm_nods": optimal_norm_nods,
         "nonoptimal_nods": nonoptimal_nods, "nonoptimal_norm_nods": nonoptimal_norm_nods,
         "mix_nods": mix_nods, "mix_norm_nods": mix_norm_nods}
    print("For chip {}".format(chip_num))
    for key, this_nods in d.items():
        fixed_nods, bad_pix = clean_nods(this_nods)
        avg, med, sum_ = nod_calcs(this_nods)
        fix_avg, fix_med, fix_sum_ = nod_calcs(fixed_nods)

        print("{} mean combined   = {}".format(key, sampled_snr(avg, chip_num)))
        print("{} mean combined   = {} (> 4 sigma removed)".format(key, sampled_snr(fix_avg, chip_num)))
        print("{} median combined = {}".format(key, sampled_snr(med, chip_num)))
        print("{} median combined = {}(> 4 sigma removed)".format(key, sampled_snr(fix_med, chip_num)))
        print("{} sum combined    = {}".format(key, sampled_snr(sum_, chip_num)))
        print("{} sum combined    = {} (> 4 sigma removed)".format(key, sampled_snr(fix_sum_, chip_num)))
        print("{} bad pixels num   = {}".format(key, len(bad_pix)))


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
