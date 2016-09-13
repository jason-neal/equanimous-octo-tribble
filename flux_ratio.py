#!/usr/bin/python

"""Brown Dwarf Flux ratio calculator:

Inputs :
Star name:
Mass of companion in Jupiter masses
Age of star: """

from __future__ import division, print_function
import numpy as np

def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Wavelength Calibrate Spectra')
    
    parser.add_argument('star_name', help='Input fits file to calibrate')

    parser.add_argument('companion_mass', help='Mass of companion')
    parser.add_argument('age', help='Star age')
    args = parser.parse_args()
    return args


def main(star_name, companion_mass, stellar_age):
    """Compute flux ratio of star to companion """


def calculate_flux_ratios(star_params, companion_params):
    """ Flux ratios for the different bands """
    f = 2.512
    Flux_ratios = dict()
    Flux_ratios["J"] = f ** (companion_params["Jmag"]-star_params["Jmag"])
    Flux_ratios["H"] = f ** (companion_params["Hmag"]-star_params["Hmag"])
    Flux_ratios["K"] = f ** (companion_params["Kmag"]-star_params["Kmag"])
    return Flux_ratios

def get_stellar_params(star_name):
    """ Astro query search """

    #return Magnitudes, parralax, Temp
    pass


def calculate_stellar_radius(star_params):
    """ Based on R/Rs = (Ts/T)^2(L/Ls)^(1/2) equation"""
    BminusV = star_params["Bmag"] - star_params["Vmag"]

    if "Teff in star_params.keys(): # need table and interpolate to this B-V
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
if __name__ == '__main__':
    args = vars(_parser())
    star_name = args.pop('star_name')
    companion_mass = args.pop('companion_mass')
    Age = args.pop('Age')
    opts = {k: args[k] for k in args}

    main(star_name, companion_mass, age, **opts)