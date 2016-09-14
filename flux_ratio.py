#!/usr/bin/python

"""Brown Dwarf Flux ratio calculator:

Inputs :
Star name:
Mass of companion in Jupiter masses
Age of star: """

from __future__ import division, print_function
import argparse
import numpy as np
from astroquery.simbad import Simbad
import pandas as pd
from PyAstronomy import pyasl
import matplotlib.pyplot as plt

def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Determine flux ratio of stellar companion')
    parser.add_argument('star_name', help='Input fits file to calibrate')
    parser.add_argument('companion_mass', help='Mass of companion (M_Jup)', type=float)
    parser.add_argument('age', help='Star age (Gyr)',type=float)
    args = parser.parse_args()
    return args


def main(star_name, companion_mass, stellar_age):
    """Compute flux ratio of star to companion """
    
    # Obtain Stellar parameters from astroquery
    star_params = get_stellar_params(star_name)    # returns a astroquesry result table
    
    # Get parameters for this mass and age 
    companion = get_brown_dwarf_information(companion_mass, stellar_age)
  
    Flux_ratios = calculate_flux_ratios(star_params, companion)

    # Compare to area ratio
    Rstar = calculate_stellar_radius(star_params)
    Rcomp_Rstar = companion["R"] / Rstar

    # Print flux ratios using a generator
    [print("{} band companion/star Flux ratio = {} ".format(key, val)) for key, val in Flux_ratios.items()]
  
    print("Radius Ratio of companion/star    = {} ".format(Rcomp_Rstar))
    print("Area Ratio of companion/star      = {} ".format(Rcomp_Rstar**2))
    
    get_sweet_cat_data()

def calculate_flux_ratios(star_params, companion_params):
    """ Flux ratios for the different bands """
    f = 2.512
    Flux_ratios = dict()
    Flux_ratios["J"] = f ** (companion_params["Mj"]-star_params["FLUX_J"])
    #Flux_ratios["H"] = f ** (companion_params["Mh"]-star_params["Hmag"])
    Flux_ratios["K"] = f ** (companion_params["Mk"]-star_params["FLUX_K"])
    return Flux_ratios


def get_stellar_params(star_name):
    """ Astro query search """

    #return Magnitudes, parralax, Temp
    customSimbad = Simbad()
    # Can add more fluxes here if need to extend to more flux ranges. although K is the limit for simbad.
    # if want higher need to search for Wise band in VISIER probably.
    customSimbad.add_votable_fields('parallax', 'sp', 'fluxdata(B)', 'fluxdata(V)', 'fluxdata(J)', 'fluxdata(K)')
   
    result_table =  customSimbad.query_object(star_name)

    #print("Table colums", result_table.colnames)
    
    return result_table


def calculate_stellar_radius(star_params):
    """ Based on R/Rs = (Ts/T)^2(L/Ls)^(1/2) equation"""
    BminusV = star_params["FLUX_B"] - star_params["FLUX_V"]
    print(BminusV, "b-v")

    #if "Fe_H_Teff" in star_params.keys(): # need table and interpolate to this B-V
    #    teff_star = star_params["Fe_H_Teff"]
    #    print(teff_star)
    #    print(teff_star[0])
    #    print("temperature fromsimbad =", teff_star)
    #else:
    #    teff_star = 0
    if True:
    #if teff_star == 0:
        #Interpolate from B-V
        bminusvs=[-0.31, -0.24, -0.20, -0.12, 0.0, 0.15, 0.29, 0.42, 0.58, 0.69, 0.85, 1.16, 1.42, 1.61] 
        temps = [34000, 23000, 18500, 13000, 9500, 8500, 7300, 6600, 5900, 5600, 5100, 4200, 3700, 3000]
        #teff_star = (4200-5100)/(1.16-0.85) * (BminusV-0.85) + 5100    # Linear interpolation
        teff_star = np.interp(BminusV, bminusvs, temps)
        print("Temperature of star from b-v", teff_star)

    Ts_T = 5800. / teff_star              # Temperature ratio
    Dm =   4.83 - star_params["FLUX_V"]   # Differnce of aboslute magnitude
    L_Ls = 2.51 ** Dm                     # Luminosity ratio
    R_Rs = (Ts_T)**2 * np.sqrt(L_Ls)      # Raidus of Star in Solar Radii

    return R_Rs   # R star in solar radii
    #BD_R = BD_Radius / R_Rs          # Radius_bd / Radius_star
    #BD_area_ratio =  BD_R**2



def get_brown_dwarf_information(companion_mass, age):
    """ baraffe 2003 table search
    Need the tables in a file somewhere"""

    mass_solar = companion_mass / 1047.56   # covert to solar mass
    BD_parameters = dict()
    
    # Find closest age model
    modelages = ["0.001", "0.005", "0.010", "0.050", "0.100", "0.120", "0.500", "1.000", "5.000", "10.000"]
    model_age = min(modelages, key=lambda x:abs(float(x)-age)) # Closest one
    model_id = "p".join(str(model_age).split("."))

    model_name = "./data/Baraffe2003/BaraffeCOND2003-" +  model_id + "Gyr.dat"
    model_data = np.loadtxt(model_name, skiprows=18, unpack=False)
    model_data = model_data.T

    cols = ["M/Ms", "Teff", "L/Ls", "g", "R", "Mv", "Mr", "Mi", "Mj", "Mh", "Mk", "Mll", "Mm"]
    for col, data in zip(cols, model_data):
        # Interpolate columns to mass of companion
        BD_parameters[col] = np.interp(mass_solar, model_data[0], data)

    return  BD_parameters  # as a dictionary


def get_sweet_cat_data():
    sc = pyasl.SWEETCat()
    data = sc.data
    print(data.head())
    print(data.keys())
    #data.plot('teff', 'vmag', kind='scatter')
    #plt.show()

    # We can also create another column for our table, e.g. luminosity so we can make a HR diagram
    # Note that this works even though we have missing values
    data['lum'] = (data.teff/5777)**4 * data.mass
    #data.plot('teff', 'lum', kind='scatter', xlim=(8500, 2500))

    # If we only are interested in a shorter table, with stars with Teff below 5000K, we can do so
    data_cool = data[data.teff < 5000]

    # See that it worked (print rows and columns)
    print(data.shape)
    print(data_cool.shape)

    # And the plot
    #plt.figure()
    #data.plot('teff', 'lum', kind='scatter', xlim=(8500, 2500))
    #plt.plot(data_cool.teff, data_cool.lum, 'oy')
    #plt.show()
    
    # Assuming given as hd******
    number = star_name[2:]
    this_entry = data[data.hd == number]
    print("this entry = ", this_entry)
   # Need to check for empty frames if the star is not in sweet-cat



if __name__ == '__main__':
    args = vars(_parser())
    star_name = args.pop('star_name')
    companion_mass = args.pop('companion_mass')
    age = args.pop('age')
    opts = {k: args[k] for k in args}

    main(star_name, companion_mass, age, **opts)