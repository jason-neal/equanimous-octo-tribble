#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Signal to Noise calculation of spectra
Jason Neal 2015 """

import numpy as np
import matplotlib.pyplot as plt
#import Obtain_Telluric as obt

def load_telluric(tapas_path, filename):
    ext = filename[-4:] 
    file_ = tapas_path + filename
    if ext == "ipac":
        with open(file_) as f:
            col1 = []
            col2 = []
            for line in f:
                firstchar = line[0]
            #print("first char =", firstchar)
                if line[0] == "\\" or line[0] == "|":
                    pass #print("Match a header line")
                else:
                    line.strip()
                    val1, val2 = line.split()
                    col1.append(val1)
                    col2.append(val2)
        tell = np.array([col1,col2], dtype="float64")
    elif ext == "fits":
        i_tell = (fits.getdata(file_,0))
        col1 = i_tell["wavelength"]
        col2 = i_tell["transmittance"]
        #print("i_tell", i_tell)
        #print("type(i_tell)", type(i_tell))
        tell = np.array([col1,col2], dtype="float64")
    else:
        print(" Could not load file", filename," with extention", ext)
        return None
    #print(tell)
    return tell  

def snr(spectra):
	snr = np.mean(spectra) / np.std(spectra, ddof=0)
	return snr

def snr_map(spectra, mindt=3, maxdt=10):
    """2D-plot of snr over the length of spectra"""

    maxdt = np.min([len(spectra), maxdt])
    dt = range(mindt+1, maxdt,1)
    #print(dt)
    store = []
    
    #for j in range(0,(len(spectra)-10)):  # shift position
    snrdt = []
    j=0

    ### TODO #### 
    #2d surface plot in python
#				_____________    
#				|	     	|
#				|	    	|
#  Delta lambda |	 snr  	|
#				|	     	| 
#				|___________|
#					lambda	
    for t in dt:
    	s = spectra[j:t+j]
    	x = snr(s)
    	snrdt.append(x)
    	#store.append(snrdt)
    	#plt.plot(dt, snrdt)
    print(store)
    #print(np.size(store))
    #plt.figure()

    plt.plot(dt, snrdt)
    plt.xlabel("Number of pixels to include in SNR")
    plt.ylabel("SNR value")
    plt.show()
    return dt



#if __name__ is "__main__":

test = load_telluric("/home/jneal/Phd/data/Tapas/", "tapas_2012-07-13T08:05:00_2680.ipac")
#test[0], test[1]
#test = [1,2,1,2,1,2,1,2,3,4,2,1,1,2,3,4,2,1,3,5,2,1,2,3,4,2,1,3,5,2,1,2,3,4,2,1,3,5,2]
#test = [t/10. + 10 for t in test]
thesnr = snr(test[1])
print("snr =", thesnr)

snr_map(test[1], mindt=3, maxdt=200)
