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
    star_params = get_stellar_params(star_name)
    
    # Get parameters for this mass and age 
    companion = get_brown_dwarf_information(companion_mass, stellar_age)
  
    Flux_ratios = calculate_flux_ratios(star_params, companion)

    # Print flux ratios using a generator
    print("{} band companion/star Flux ratio = {} ".format(key, val) for key, val in Flux_ratios.unpac())

    # Compare to area ratio
    Rstar = calculate_stellar_radius(star_params)
    Rcomp_Rstar = companion["Radius"] / Rstar
    print("Radius Ratio of companion/star    = {} ".format(Rcomp_Rstar))
    print("Area Ratio of companion/star      = {} ".format(Rcomp_Rstar**2))


def calculate_flux_ratios(star_params, companion_params):
    """ Flux ratios for the different bands """
    f = 2.512
    Flux_ratios = dict()
    Flux_ratios["J"] = f ** (companion_params["Mj"]-star_params["Jmag"])
    Flux_ratios["H"] = f ** (companion_params["Mh"]-star_params["Hmag"])
    Flux_ratios["K"] = f ** (companion_params["Mk"]-star_params["Kmag"])
    return Flux_ratios

def get_stellar_params(star_name):
    """ Astro query search """

    #return Magnitudes, parralax, Temp
    Simbad
    columns = ["HD", "SpType", "B-V", "Vmag"]
    result_table = Simbad.query_object(star_name, columns=columns)
    print("Table colums", result_table.columns)
    
    


def calculate_stellar_radius(star_params):
    """ Based on R/Rs = (Ts/T)^2(L/Ls)^(1/2) equation"""
    BminusV = star_params["Bmag"] - star_params["Vmag"]

    if Teff in star_params.keys(): # need table and interpolate to this B-V
        teff_star
    else: #
        #Interpolate from B-V
        teff_star = (4200-5100)/(1.16-0.85) * (BminusV-0.85) + 5100    # Linear interpolation

    Ts_T = 5800. / teff_star    # Temperature ratio

    Dm =   4.83 - star_params["AbsVmag"] # Differnce of aboslute magnitude
    L_Ls = 2.51 ** Dm                    # Luminosity ratio

    R_Rs = (Ts_T)**2*(L_Ls)**0.5    # Raidus of Star in Solar Radii

    return R_Rs   # R star in solar radii
#BD_R = BD_Radius / R_Rs          # Radius_bd / Radius_star
#BD_area_ratio =  BD_R**2


    pass


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




    

if __name__ == '__main__':
    args = vars(_parser())
    star_name = args.pop('star_name')
    companion_mass = args.pop('companion_mass')
    Age = args.pop('Age')
    opts = {k: args[k] for k in args}

    main(star_name, companion_mass, age, **opts)