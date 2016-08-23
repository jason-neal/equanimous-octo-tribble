#!/usr/bin/env python3
# -*- coding: utf8 -*-

""" 
Rename tapas spectra based on the some of the header information in it.

"""

import argparse
from astropy.io import fits
from time import strftime, strptime
import os
import subprocess
import Obtain_Telluric as obt


def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Xml generator for tapas submission')
    parser.add_argument("fname", help='Input file name', type=str)
    parser.add_argument("-x", "--extract", help='Extract the file from gz first', action="store_true")
    parser.add_argument("-p", "--prefix", help="Prefix to add", type=str)
    parser.add_argument("-k", "--keep_original", help="Keep the original file", action="store_true")
    parser.add_argument("-s", "--species", help="Specify constituents if known e.g. h20, noh20", type=str)

    args = parser.parse_args()
    return args


def main(fname, extract=False, prefix="", keep_original=False, species=False):
    
    split_fname = fname.split("/")
    if len(split_fname) == 1:
        # No path 
        filename = fname
        path = ""
    else:
        filename = fname.split("/")[-1]
        path_folders = fname.split("/")[:-1]
        path = "/".join(path_folders) + "/"
    
    ext = filename.split(".")[-1]

    if keep_original:
        old_fname = fname
        if ext == "gz":
            fname = path + "tmp_file." + filename.split(".")[-2] + ".gz"
        else:
            fname = path + "tmp_file." + ext   
        subprocess.call(["cp", old_fname, fname])
    
    if extract and ext == ".gz":
        subprocess.call(["gunzip ", fname], shell=True)
        print(filename)
        filename = filename[:-3]
        print(filename)
        ext = filename.split(".")[-1]
    elif ext == "gz":
        print("The file ends in .gz. Add -x to call to extract it first")
        #raise
    data, hdr = obt.load_telluric(path, filename)

    date = strftime("%Y-%m-%dT%H-%M-%S", strptime(hdr["DATE-OBS"], "%Y/%m/%d %H:%M:%S"))
    req_id = int(float(hdr["req_id"]))  # Id of request 
    berv = hdr["BARYDONE"]
    
    try:
        sampling = int(float(hdr["SAMPRATI"]))
    except:
        sampling = -1
   
    try:
        respower = int(float(hdr["RESPOWER"]))
    except:
        respower = -1

    if species:
        new_name = "{}_ReqId_{}_R-{}_sratio-{}_barydone-{}_species-{}".format(date, req_id, respower, sampling, berv, species)
    elif sampling is -1 and respower is -1:
        new_name = "{}_ReqId_{}_No_Ifunction_barydone-{}".format(date, req_id, berv)
    else:
        new_name = "{}_ReqId_{}_R-{}_sratio-{}_barydone-{}".format(date, req_id, respower, sampling, berv)
  
    if prefix:
        new_name = path + "tapas_" + prefix + "_" + new_name + "." + ext
    else:
        new_name = path + "tapas_" + new_name + "." + ext
    print("Name for file = ", new_name)

    subprocess.call(["mv", fname, new_name])
    print("Renamed {} to {}".format(fname, new_name))


if __name__ == "__main__":
    args = vars(_parser())
    fname = args.pop('fname')
    opts = {k: args[k] for k in args}
    
    print(opts)

    main(fname, **opts)
