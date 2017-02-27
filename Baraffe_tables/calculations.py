"""Calculations for flux ratios."""
import numpy as np
from db_queries import get_temperature
from typing import Dict, List, Any


def flux_mag_ratio(mag1: float, mag2: float) -> float:
    """Calculate the flux ratio between two magnitudes.

    Using the equation f1/f2 = 10**(-0.4 * (mag1 - mag2)).
    A common approximation is the equation f1/f2 approx 2.512**(m2 - m1), but is not here.

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
    if isinstance(mag1, float) and isinstance(mag2, float):
        # return 2.512**(mag2 - mag1)  # Approximation
        return 10**(-0.4 * (mag1 - mag2))
    else:
        print(type(mag1), type(mag2))
        raise ValueError("Magnitudes given were not floats. Mag1 = {}, Mag2 = {}".format(mag1, mag2))


def calculate_flux_ratio(star_params: Any, companion_params: Dict[str, float], bands: List[str]) -> Dict[str, float]:
    """Flux ratios for the different bands in bands.

    Parameters
    ----------
    star_params: dict
       Stellar parameters with stellar magnitude values like "FLUX_K".
    companion_params: dict
        Companion parameters with magnitude values like "Mk".
    bands: list of str
        Bands to return ratios for.

    Returns
    -------
    flux_ratios: dict
       Flux ratios for each given band.

    """
    flux_ratios = dict()
    for band in bands:
        band = band.upper()
        flux_ratios[band] = flux_mag_ratio(float(star_params["FLUX_{0!s}".format(band)][0]),
                                           companion_params["M{}".format(band.lower())])

    return flux_ratios


def calculate_stellar_radius(star_params: Any) -> float:
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
    Dm = 4.83 - star_params["FLUX_V"][0]     # Differnce of aboslute magnitude
    L_Ls = 2.51 ** Dm                     # Luminosity ratio
    R_Rs = (Ts_T)**2 * np.sqrt(L_Ls)      # Raidus of Star in Solar Radii

    return R_Rs   # Radius of star in solar radii


def calculate_companion_magnitude(star_params: Any, flux_ratio: float, bands: List[str]=["K"]) -> Dict[str, float]:
    """Calculate companion magnitude from flux ratio.

    Using the equation m - n = -2.5 * log_10(F_m / F_n).

    Parameters
    ----------
    star_params: dict
        Parameters for the host star.
    flux_ratio: float
        Flux ratio for the system (F_companion/F_host).
    band: str
        Bands to use. default = "K"

    Returns
    -------
    magnitudes: dict
        Magnitudes for the companion in the J, H, and K bands.

    Note
    ----
    This is possibly not quite the correct implemenation as we are
    only using a single flux_ratio value.

    """
    magnitudes = dict()

    for band in bands:
        band = band.upper()
        if band in ["J", "H", "K"]:
            magnitudes[band] = star_params["FLUX_{0!s}".format(band)] - 2.5 * np.log10(flux_ratio)
        else:
            raise ValueError("Magnitude band {0!s} was not in parameter dictionaries.".format(band))
        print("Band in calc magnitude ", band, "mags", magnitudes)

    return magnitudes
