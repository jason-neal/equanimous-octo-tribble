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
# import matplotlib.pyplot as plt

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
    Rstar = calculate_stellar_radius(star_name, star_params)
    Rcomp_Rstar = companion["R"] / Rstar

    # Print flux ratios using a generator
    print("Magnitude Calculation\n")
    [print("{0} band star/companion Flux ratio = {1} >>> companion/star Flux ratio {2}".format(key, val[0], 1./val[0])) for key, val in Flux_ratios.items()]

    print("\nRadius Calculation")
    print("Star radius      = {} R_sun".format(Rstar[0]))
    print("Planet radius    = {} R_sun".format(np.round(companion["R"], 4)))
    print("Radius Ratio of companion/star    = {} ".format(Rcomp_Rstar[0]))
    print("Area Ratio of companion/star      = {} ".format(Rcomp_Rstar[0]**2))


##############################################################################
# Access Databases
##############################################################################


def get_stellar_params(star_name):
    """ Astro query search """

    #return Magnitudes, parralax, Temp
    customSimbad = Simbad()
    # Can add more fluxes here if need to extend to more flux ranges. although K is the limit for simbad.
    # if want higher need to search for Wise band in VISIER probably.
    customSimbad.add_votable_fields('parallax', 'sp', 'fluxdata(B)', 'fluxdata(V)', 'fluxdata(J)', 'fluxdata(K)', 'fe_h')

    result_table =  customSimbad.query_object(star_name)

    #print("Table colums", result_table.colnames)

    return result_table


def get_brown_dwarf_information(companion_mass, age):
    """ baraffe 2003 table search
    Need the tables in a file somewhere"""

    mass_solar = companion_mass / 1047.56   # covert to solar mass
    BD_parameters = dict()

    # Find closest age model
    modelages = ["0.001", "0.005", "0.010", "0.050", "0.100", "0.120", "0.500", "1.000", "5.000", "10.000"]
    model_age = min(modelages, key=lambda x:abs(float(x)-age)) # Closest one
    model_id = "p".join(str(model_age).split("."))

    model_name = "./Baraffe2003/BaraffeCOND2003-" +  model_id + "Gyr.dat"
    model_data = np.loadtxt(model_name, skiprows=18, unpack=False)
    model_data = model_data.T

    cols = ["M/Ms", "Teff", "L/Ls", "g", "R", "Mv", "Mr", "Mi", "Mj", "Mh", "Mk", "Mll", "Mm"]
    for col, data in zip(cols, model_data):
        # Interpolate columns to mass of companion
        BD_parameters[col] = np.interp(mass_solar, model_data[0], data)

    return  BD_parameters  # as a dictionary


def get_sweet_cat_temp(star_name):
    sc = pyasl.SWEETCat()
    data = sc.data
    #print(data.head())
    #print(data.keys())

    # Assuming given as hd******
    hd_number = star_name[2:]
    #print("hd number ", hd_number)
    if hd_number in sc.data.hd.values:
        hd_entry = data[data.hd == hd_number]

        if hd_entry.empty:
            return False
        else:
        # Sweet-cat has temperature of zero
        #if it does not have a temperature value for the star
            return hd_entry.iloc[0]["teff"]
    else:
        print("This star not in SWEET-Cat")
        return False









##############################################################################
# Calculations
##############################################################################
def calculate_flux_ratios(star_params, companion_params):
    """ Flux ratios for the different bands """
    f = 2.512
    Flux_ratios = dict()
    Flux_ratios["J"] = f ** (companion_params["Mj"]-star_params["FLUX_J"])
    #Flux_ratios["H"] = f ** (companion_params["Mh"]-star_params["Hmag"])
    Flux_ratios["K"] = f ** (companion_params["Mk"]-star_params["FLUX_K"])
    return Flux_ratios

def get_temperature(star_name, star_params):
    """ Try get temperature of star multiple ways

    1st - Try Fe_H_Teff param from Simbad.
    2nd - Try SweetCat but the star might not be there
    3rd - Calculate from B-V and interpolation"""
    good_temp = False
    # This is not the best way to do this but work atm
    if "Fe_H_Teff" in star_params.keys(): # need table and interpolate to this B-V
        #print("star_params['Fe_H_Teff'] =", star_params["Fe_H_Teff"])
        teff = star_params["Fe_H_Teff"][0]
        if teff == 0 or teff == [0]:
            # No teff given by Simbad
            print("Simbad Temperature was zero")
            teff = None
        else:
            good_temp = True
            print("Temperature obtained from Fe_H_Teff", good_temp)
            return teff

    if not good_temp:
        try:
            teff = get_sweet_cat_temp(star_name)

            if teff == 0 or np.isnan(teff):  # temp from sweet-cat
                print("No SWEET-Cat temperature, teff was", teff)
                teff = None
            else:
                print("SWEET-Cat teff =", teff)
                good_temp = True
                return teff
        except:
            print("Failed to get Sweetcat_temp")
            good_temp = False
            raise

    if not good_temp:
        """Then use the B-V method as last resort"""
        BminusV = star_params["FLUX_B"] - star_params["FLUX_V"]
        #print(BminusV, "b-v")
        #Interpolate from B-V
        bminusvs = np.array([-0.31, -0.24, -0.20, -0.12, 0.0, 0.15, 0.29, 0.42, 0.58, 0.69, 0.85, 1.16, 1.42, 1.61])
        temps = np.array([34000, 23000, 18500, 13000, 9500, 8500, 7300, 6600, 5900, 5600, 5100, 4200, 3700, 3000])
        #teff_star = (4200-5100)/(1.16-0.85) * (BminusV-0.85) + 5100    # Linear interpolation
        teff = np.interp(BminusV, bminusvs, temps)[0]
        print("Temperature of star was calculated from b-v = {} K".format(teff))
    return teff


def calculate_stellar_radius(star_name, star_params):
    """ Based on R/Rs = (Ts/T)^2(L/Ls)^(1/2) equation"""

    teff_star = get_temperature(star_name, star_params)

    Ts_T = 5800. / teff_star              # Temperature ratio
    Dm =   4.83 - star_params["FLUX_V"]   # Differnce of aboslute magnitude
    L_Ls = 2.51 ** Dm                     # Luminosity ratio
    R_Rs = (Ts_T)**2 * np.sqrt(L_Ls)      # Raidus of Star in Solar Radii

    return R_Rs   # R star in solar radii


if __name__ == '__main__':
    args = vars(_parser())
    star_name = args.pop('star_name')
    companion_mass = args.pop('companion_mass')
    age = args.pop('age')
    opts = {k: args[k] for k in args}

    main(star_name, companion_mass, age, **opts)
