#!/usr/bin/env python
# -*- coding: utf8 -*-

"""Combine Normalized DRACS Specta

Script that take one imput parameter which is the path to the normobs directory

Needs functions/modules to get the filenames, export to a fits tablefile, change fits hdr 
"""

import time
import argparse
import os
import fnmatch
import matplotlib.pyplot as plt
from astropy.io import fits
import IOmodule
import numpy as np
import scipy as sp 


def get_filenames(path, regexp, regexp2=False):
    """ regexp must be a regular expression as a string
            eg '*.ms.*', '*_2.*', '*.ms.norm.fits*'
        resexp2 is if want to match two expressions such as 
            '*_1*' and '*.ms.fits*'
    """
    os.chdir(path)
    filelist = []
    for file in os.listdir('.'):
        if regexp2:  # Match two regular expresions
            if fnmatch.fnmatch(file, regexp) and fnmatch.fnmatch(file, regexp2):
                #print file
                filelist.append(file)
        else:
            if fnmatch.fnmatch(file, regexp):
                #print file
                filelist.append(file)
    filelist.sort()
    return filelist

def SumNods(Spectra, Headers, Pos="All", Norm="None"):
    """ Add together the nod postitions of Crires spectra"""
    """ Currently this implements baised on expected nod ordering, 
    If want to be more sure the correct nods are added the headers
    should be checked for its nod position"""
    if Pos == "All":
        NodNum = len(Spectra)
    else:
        NodNum = len(Spectra)/2.0
    #print("number of nods ", NodNum)
    #if Headers[i]["HIERARCH ESO SEQ NODPOS"] == Pos:
    NodSum = np.zeros_like(Spectra[0])
    if Pos.upper() == "ALL":
        # Sum all
        #NodNum = 8
        for Spec in Spectra:
            NodSum += np.array(Spec)

    elif Pos.upper() == "A":
        # Sum A nods
        #NodNum = 4
        #for i in [0,3,4,7]:
           # NodSum += np.array(Spectra[i])
        for Spec, hdr in zip(Spectra, Headers):
            if hdr["HIERARCH ESO SEQ NODPOS"] == "A":
                NodSum += np.array(Spec)

    elif Pos.upper() == "B":   
        # Sum B nods
        #NodNum = 4
        #for i in [1,2,5,6]:
            #NodSum += np.array(Spectra[i])
        for Spec, hdr in zip(Spectra, Headers):
            if hdr["HIERARCH ESO SEQ NODPOS"] == "B":
                NodSum += np.array(Spec)    
    if Norm.upper() == "MEDIAN":
      NodSum /= np.median(NodSum)
    elif Norm.upper() == "MEAN":
      NodSum /= np.mean(NodSum)  
    elif Norm.upper() == "DIVIDE":
      NodSum /= NodNum  
    return NodSum 


def ExportToFits(Outputfile, Norm_All, NodA, NodB, hdr, hdrkeys, hdrvals):
    """ Write Combined DRACS CRIRES NOD Spectra to a fits table file"""
    col1 = fits.Column(name="Combined", format="E", array=Norm_All) # colums of data
    col2 = fits.Column(name="Nod A", format="E", array=NodA)
    col3 = fits.Column(name="Nod B", format="E", array=NodB)
    cols = fits.ColDefs([col1, col2, col3])
    tbhdu = fits.BinTableHDU.from_columns(cols) # binary tbale hdu
    prihdr = append_hdr(hdr, hdrkeys, hdrvals)
    prihdu = fits.PrimaryHDU(header=prihdr)
    thdulist = fits.HDUList([prihdu, tbhdu])
    #print("Writing to fits file")
    thdulist.writeto(Outputfile,output_verify="silentfix")   # Fixing errors to work properly
    return None


def append_hdr(hdr, keys, values ,item=0):
    ''' Apend/change parameters to fits hdr, 
    can take list or tuple as input of keywords 
    and values to change in the header 
    Defaults at changing the header in the 0th item 
    unless the number the index is givien,
    If a key is not found it adds it to the header'''
    # open fits file
    #hdulist = fits.open(output)
    #hdr = hdulist[item].header
    #print repr(hdr[0:10])
    #assert type(keys) == type(values), 'keys and values do not match'
    if type(keys) == str:           # To handle single value
        hdr[keys] = values
    else:
        assert len(keys) == len(values), 'Not the same number of keys as values' 
        for i in range(len(keys)):
            hdr[keys[i]] = values[i]
            print repr(hdr[-2:10])
    return hdr

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Combine CRIRES Nod positions \
                                     Normized values')
    parser.add_argument('inputpath', help='Path to normobs directory')
    #parser.add_argument('-s', '--sun', help='Plot with spectra of the Sun (1)',
    #                     default=False)
    args = parser.parse_args()

    #print args
    path = args.inputpath
    # find norm values in this directory:

    chips = range(4)
    for chip in chips:
    #while True:
        #chip = 1
        org_vals = get_filenames(path, "CRIRE*.ms.fits", "*_" + str(chip + 1) + "*")
        norm_vals = get_filenames(path, "CRIRE*.ms.norm.fits", "*_" + str(chip + 1) + "*")

        I_dracs = []
        I_dracs_hdrs = []
        I_norm = []
        I_norm_hdrs = []
        for name in org_vals:
            ThisFile = path + name
            I_dracs.append(fits.getdata(ThisFile,0))
            I_dracs_hdrs.append(fits.getheader(ThisFile,0))
        
        dracs_All = SumNods(I_dracs, I_dracs_hdrs, Pos="All", Norm="Divide")
        dracs_A = SumNods(I_dracs, I_dracs_hdrs, Pos="A", Norm="Divide")
        dracs_B = SumNods(I_dracs, I_dracs_hdrs, Pos="B", Norm="Divide")
         
        #print(norm_vals)
        for name in norm_vals:
                ThisFile = path + name
                ThisNorm = fits.open(ThisFile)
                Last_normhdr = ThisNorm[0].header
                I_norm_hdrs.append(Last_normhdr)
                I_norm.append(ThisNorm[0].data)
                ThisNorm.close()

        print("Inorm",I_norm)

        NormalizeMethod = "Divide"
        Norm_All = SumNods(I_norm, I_norm_hdrs, Pos="All", Norm=NormalizeMethod)
        Norm_A = SumNods(I_norm, I_norm_hdrs, Pos="A", Norm=NormalizeMethod)
        Norm_B = SumNods(I_norm, I_norm_hdrs, Pos="B", Norm=NormalizeMethod)

#        plt.plot(dracs_All, label="dracs All")
#        plt.plot(dracs_A, label="dracs A")
#        plt.plot(dracs_B, label="dracs B")
#        plt.legend()
#        plt.show()

#        plt.plot(Norm_All, label="All")
#        plt.plot(Norm_A, label="A")
#        plt.plot(Norm_B, label="B")
#        plt.legend()
#        plt.show()

        # write ouput to fits file
#        testhdr = fits.Header()
#        testhdr['TESTVal'] = 'Edwin Hubble'
        T_Now = str(time.gmtime()[0:6])
#        testhdr['Date'] = (T_Now, ' Time fits was last changed')        

        outputfile = path + norm_vals[-1][0:-4] + "comb.fits"   
        print("Output file name", outputfile)
        hdrkeys = ["Description"]
        hdrvals = ["Combine DRACS Nomalized CRIRES Nod Spectra"]
        hdrkeys = ['COMBINEDATE', "COMBINEMETHOD", "COMBNORMALIZE"]
        hdrvals = [(T_Now, "Time Nods were combined"), ("Addition", "How nods were added"), \
                  (NormalizeMethod, "Method for normalizing combined spectra (divide/mean/median)")]
    
        for i, name in zip(range(len(org_vals)), org_vals):
            hdrkeys.append("CRIRES NOD FILE " + str(i + 1))
            hdrvals.append((name, "Input filename"))
        ExportToFits(outputfile, Norm_All, Norm_A, Norm_B, Last_normhdr, hdrkeys, hdrvals)
        print("Wrote to fits Succesfully")

        #break


    ## Try load in the combined nod spectra fits files now














