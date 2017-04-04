"""Access Databases."""
import numpy as np
from PyAstronomy import pyasl
from astroquery.simbad import Simbad
from typing import Any, Union, Optional


def get_stellar_params(star_name: str) -> Any:
    """"Astroquery SIMBAD search for stellar parameters.

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


def get_sweet_cat_temp(star_name: str) -> Union[bool, float, int]:
    """Obtain spectroscopic temperature from SWEET-Cat.

    Parameters
    ----------
    star_name: str
        Star identifier. HD number only accepted currently.

    """
    sc = pyasl.SWEETCat()
    data = sc.data

    if star_name[0:2].lower() != "hd":
        # only accept HD numbers atm
        raise NotImplementedError

    # Assuming given as hd******
    hd_number = star_name[2:]
    # print("hd number ", hd_number)
    if hd_number in sc.data.hd.values:
        hd_entry = data[data.hd == hd_number]

        if hd_entry.empty:
            return False
        elif (hd_entry.iloc[0]["teff"] != 0) and (not np.isnan(hd_entry.iloc[0]["teff"])):
            # Temp = 0 when doesn't exist
            return hd_entry.iloc[0]["teff"]
        else:
            return False
    else:
        print("{!s} was not in SWEET-Cat.".format(star_name))
        return False


def get_temperature(star_name: str, star_params: Optional[Any]=None) -> float:
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
        teff = get_sweet_cat_temp(star_name)

        if (teff is False) or (teff == 0) or (np.isnan(teff)):  # temp from sweet-cat
            print("No SWEET-Cat temperature, teff was {0} K".format(teff))
            teff = None
        else:
            print("SWEET-Cat teff = {0:.0f} K".format(teff))
            good_temp = True
            return teff

    if not good_temp:
        print("Using the B-V method as last resort.")
        teff = calculate_bv_temp(star_params["FLUX_B"], star_params["FLUX_V"])
        print("Temperature of star was calculated from b-v ~= {} K".format(teff))
    return teff


def calculate_bv_temp(b_mag: float, v_mag: float) -> float:
    """Calcualte Stellar Temperature from B-V magnitudes.

    Parameters
    ----------
    b_mag: float
        Stellar B magnitude.
    v_mag: float
        Stellar V magnitude.
    Returns
    -------
    temp: float
       Temperature value in Kelvin.

    """
    b_v = b_mag - v_mag

    # Interpolate from B-V
    bminusvs = np.array([-0.31, -0.24, -0.20, -0.12, 0.0, 0.15, 0.29,
                         0.42, 0.58, 0.69, 0.85, 1.16, 1.42, 1.61])
    temps = np.array([34000, 23000, 18500, 13000, 9500, 8500, 7300,
                      6600, 5900, 5600, 5100, 4200, 3700, 3000])

    return np.interp(b_v, bminusvs, temps)[0]
