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


    if __name__ == "__main__":
        args = vars(_parser())

        opts = {k: args[k] for k in args}
        sys.exit(main(**opts))
