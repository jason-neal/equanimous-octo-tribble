
##############################################################################
# Access Databases
##############################################################################
import numpy as np
from astroquery.simbad import Simbad
from PyAstronomy import pyasl


def get_stellar_params(star_name):
    """Astroquery SIMBAD search for stellar parameters.

    Parameters
    ----------
    star_name: str
        Stellar name to get parameters for.

    Returns
    -------
    result_table: votable, dict-like
    """

    # return Magnitudes, parralax, Temp
    customSimbad = Simbad()
    # Can add more fluxes here if need to extend flux ranges. Although K is the simbad limit.
    # if want higher need to search for Wise band in VISIER probably.
    customSimbad.add_votable_fields('parallax', 'sp', 'fluxdata(B)',
                                    'fluxdata(V)', 'fluxdata(J)', 'fluxdata(H)', 'fluxdata(K)',
                                    'fe_h')

    result_table = customSimbad.query_object(star_name)

    # Add star name to parameters
    result_table["name"] = star_name

    return result_table



def get_sweet_cat_temp(star_name):
    """Obtain spectroscopic temperature from SWEET-Cat.

    Parameters
    ----------
    star_name: str
        Star identifier. HD# only accepted currently."""
    sc = pyasl.SWEETCat()
    data = sc.data

    if star_name[0:2].lower != "hd":
        # only accept HD numbers atm
        raise NotImplementedError

    # Assuming given as hd******
    hd_number = star_name[2:]
    # print("hd number ", hd_number)
    if hd_number in sc.data.hd.values:
        hd_entry = data[data.hd == hd_number]

        if hd_entry.empty:
            return False
        elif hd_entry.iloc[0]["teff"] != 0:
            # Temp = 0 when doesn't exist
            return hd_entry.iloc[0]["teff"]
        else:
            return False
    else:
        print("{!s} was not in SWEET-Cat.".format(star_name))
        return False


def get_temperature(star_name, star_params=None):
    """Find temperature of the star multiple ways.

    1st - Try Fe_H_Teff param from Simbad.
    2nd - Try SweetCat but the star might not be there (only planet hosts).
    3rd - Calculate from B-V and interpolation.

    """
    if star_params is None:
        star_params = get_stellar_params(star_name)

    good_temp = False
    # This is not the best way to do this but work atm
    if "Fe_H_Teff" in star_params.keys():  # need table and interpolate to this B-V
        # print("star_params['Fe_H_Teff'] =", star_params["Fe_H_Teff"])
        teff = star_params["Fe_H_Teff"][0]
        if teff == 0 or teff == [0]:
            # No teff given by Simbad
            print("Simbad Temperature was zero.")
            teff = None
        else:
            good_temp = True
            print("Temperature obtained from Fe_H_Teff = {0:5.0f} K".format(good_temp))
            return teff

    if not good_temp:
        try:
            teff = get_sweet_cat_temp(star_name)

            if teff == 0 or np.isnan(teff):  # temp from sweet-cat
                print("No SWEET-Cat temperature, teff was {0} K".format(teff))
                teff = None
            else:
                print("SWEET-Cat teff = {0:.0f} K".format(teff))
                good_temp = True
                return teff
        except:
            print("Failed to get temperature from SweetCat.")
            good_temp = False
            # raise

    if not good_temp:
        print("Using the B-V method as last resort.")
        teff = bv_temp(star_params)
        print("Temperature of star was calculated from b-v ~= {} K".format(teff))
    return teff


def bv_temp(star_params):
    """Stellar Temperature from B-V magnitudes.

    Parameters
    ----------
    star_params: dict-like
        Simbad votable with FLUX_B and FLUX_V parameters.

    Returns
    -------
    temp: float
       Temperature value in Kelvin.
    """
    BminusV = star_params["FLUX_B"] - star_params["FLUX_V"]

    # Interpolate from B-V
    bminusvs = np.array([-0.31, -0.24, -0.20, -0.12, 0.0, 0.15, 0.29,
                         0.42, 0.58, 0.69, 0.85, 1.16, 1.42, 1.61])
    temps = np.array([34000, 23000, 18500, 13000, 9500, 8500, 7300,
                      6600, 5900, 5600, 5100, 4200, 3700, 3000])

    return np.interp(BminusV, bminusvs, temps)[0]


def calculate_stellar_radius(star_params):
    """Based on R/Rs = (Ts/T)^2(L/Ls)^(1/2) equation.

    Parameters
    ----------
    star_params: votable, dict
        Table of Stellar parameters.

    Returns
    -------
    R_Rs: float
        Esitmated Stellar Radius in solar radii.

    """
    star_name = star_params['name'][0]
    teff_star = get_temperature(star_name, star_params)

    Ts_T = 5800. / teff_star              # Temperature ratio
    Dm = 4.83 - star_params["FLUX_V"]     # Differnce of aboslute magnitude
    L_Ls = 2.51 ** Dm                     # Luminosity ratio
    R_Rs = (Ts_T)**2 * np.sqrt(L_Ls)      # Raidus of Star in Solar Radii

    return R_Rs   # Radius of star in solar radii
