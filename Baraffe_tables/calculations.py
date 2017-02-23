
# ############################################################################
# Calculations for flux ratios.
# ############################################################################
import numpy as np
from db_queries import get_temperature


def flux_mag_ratio(mag1, mag2):
    """Calcualte the flux ratio between two magnitudes.

    Using the equation f1/f2 = 10**(-0.4 * (mag1 - mag2)).
    A common approximation is the equation f1/f2 \approx 2.512**(m2 - m1), but is not here.

    Parameters
    ----------
    mag1 float
        Magnitude of first object.
    mag2: float
        Magnitude of second object.

    Returns
    -------
    flux_ratio: float
        flux/contrast ratio between the two magnitudes.

    """
    # return 2.512**(mag2 - mag1)  # Approximation
    return 10**(-0.4 * (mag1 - mag2))


def calculate_flux_ratios(star_params, companion_params):
    """Flux ratios for the different bands."""
    Flux_ratios = dict()
    Flux_ratios["J"] = flux_mag_ratio(star_params["FLUX_J"], companion_params["Mj"])
    Flux_ratios["K"] = flux_mag_ratio(star_params["FLUX_K"], companion_params["Mk"])
    Flux_ratios["H"] = flux_mag_ratio(star_params["FLUX_H"], companion_params["Mh"])
    return Flux_ratios

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


def calculate_companion_magnitude(star_params, flux_ratio, band="K"):
    """ Calculte companion magnitude from flux ratio

    Using the equation m - n = -2.5 * log_10(F_m / F_n)

    Parameters
    ----------
    star_params: dict
        Parameters for the host star.
    flux_ratio: float
        Flux ratio for the system (F_companion/F_host).
    band: str
        Band to use. default = "K"

    Returns
    -------
    magnitudes: dict
        Magnitudes of the companion in J, H, and K bands.

    """
    magnitudes = dict()
    band = band.upper()
    if band in ["J", "H", "K"]:
        magnitudes[band] = star_params["FLUX_{}".format(band)] - 2.5 * np.log10(flux_ratio)
    else:
        return ValueError("Band is not available atm.")
    return magnitudes
