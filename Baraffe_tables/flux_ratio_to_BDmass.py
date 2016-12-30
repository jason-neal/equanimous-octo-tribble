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
    parser.add_argument('star_name', help='Name of host star.', type=str)
    parser.add_argument('flux_ratio', help='Flux ratio between host and companion - (F_comp/F_host)', type=float)
    parser.add_argument('age', help='Star age (Gyr)',type=float)
    parser.add_argument("-b", "--band", help='Magnitude band of the flux ratio', choices=["J", "K"], default="K", type=str)
    args = parser.parse_args()
    return args


def main(star_name, flux_ratio, stellar_age, band="K"):
    """Compute flux ratio of star to companion """
    Jup_mass = 1047.56   # covert to Solar to Jupiter mass
    # test_mag_conversions()
    # Obtain Stellar parameters from astroquery
    star_params = get_stellar_params(star_name)    # returns a astroquesry result table

    # Calculate Magnitude J and K band magnitude for this flux ratio
    magnitudes = calculate_companion_magnitude(star_params, flux_ratio)
    print("Magnitude calculate for companion", magnitudes)
    # Find companion parameters that match these magnitudes
    companion_params = get_BD_from_flux_ratio(magnitudes, stellar_age, band=band)

    # Print flux ratios using a generator
    print("Estimated Companion Mass from {} band Flux ratio\n".format(band.upper()))
    print("M/M_S = {0} (M_star) = {1} (M_Jup), Temp = {2} K".format(companion_params["M/Ms"][0], Jup_mass * companion_params["M/Ms"][0], companion_params["Teff"][0]))


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

    return BD_parameters  # as a dictionary


def get_BD_from_flux_ratio(magnitudes, age, band="K"):
    """ Baraffe 2003 table search for given magnitude.
    """
    # mass_solar = companion_mass / 1047.56   # covert to solar mass
    BD_parameters = dict()

    # Find closest age model
    modelages = ["0.001", "0.005", "0.010", "0.050", "0.100", "0.120", "0.500", "1.000", "5.000", "10.000"]
    model_age = min(modelages, key=lambda x:abs(float(x)-age)) # Closest one
    model_id = "p".join(str(model_age).split("."))

    model_name = "./Baraffe2003/BaraffeCOND2003-" +  model_id + "Gyr.dat"
    model_data = np.loadtxt(model_name, skiprows=18, unpack=False)
    model_data = model_data.T

    cols = ["M/Ms", "Teff", "L/Ls", "g", "R", "Mv", "Mr", "Mi", "Mj", "Mh", "Mk", "Mll", "Mm"]

    if "M"+band.lower() not in cols:
        raise ValueError("Band {} is not in Baraffe tables.".format(band))
    if band not in magnitudes.keys():
        raise ValueError("The band {} given is not in the given magnitudes".format(band))

    band_index = [i for i, c in enumerate(cols) if "M"+band.lower() == c][0]
    # print("band_index", band_index)

    print("Band mag", magnitudes[band])
    print(model_data)
    print("model data", model_data[band_index])
    print("all data")
    for col, data in zip(cols, model_data):
        # Interpolate columns to magnitude of companion
        print("\nThis col =", col)
        print("Band mag", magnitudes[band])
        #print(model_data)
        print("model data", model_data[band_index])
        print("this col data ", data)

        BD_parameters[col] = np.interp(magnitudes[band], model_data[band_index], data)
        print("This col interp value =", BD_parameters[col] )

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

def calculate_companion_magnitude(star_params, flux_ratio, band="K"):
    """ Calculte companion magnitude from flux ratio

    Using the equation m - n = -2.5 * log_10(F_m / F_n)

    Parameters
    ----------
    star_params: dict
        Parameters for the host star.
    flux_ratio: float
        Flux ratio for the system (F_comp/F_host).
    band: str
        Band to use. default = "K"

    Returns
    -------
    magnitudes: dict
        Magnitudes of the companion in J and K bands.

    """
    magnitudes = dict()
    band = band.upper()
    if band in ["J", "K"]:
        magnitudes[band] = star_params["FLUX_{}".format(band)] - 2.5 * np.log10(flux_ratio)
    else:
        return ValueError("Band is not available atm.")
    return magnitudes


def test_mag_conversions():
    """ Test converting from flux ratio to magnitude back to flux ratio etc.
    Tests show the conversion goes back and forward."""

    vals=dict()
    star=dict()
    vals["Mj"] = 4
    vals["Mk"] = 6
    star["FLUX_J"] = 3
    star["FLUX_K"] = 4
    ratios = calculate_flux_ratios(star, vals)
    print("Ratios from function", ratios)
    magnitudes = calculate_companion_magnitude(star, 1./ratios["J"])
    print("magnitudes from ratios", magnitudes)
    new_ratios = calculate_flux_ratios(star, magnitudes)
    print("new_ratios from mags", new_ratios)


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
    flux_ratio = args.pop('flux_ratio')
    age = args.pop('age')
    opts = {k: args[k] for k in args}

    main(star_name, flux_ratio, age, **opts)