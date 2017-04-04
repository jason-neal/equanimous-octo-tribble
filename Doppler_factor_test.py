# Test Doppler shift

# Astropy specutils uses the doopler factor to doppler shift.
# I want to test the diference between doppler factor and non-relativistic version.

from __future__ import division, print_function
import numpy as np
from astropy import constants
from astropy import units as u

c = 299792458 * u.meter / u.second

wavelength = np.arange(400, 5000) * u.nanometer

rv = 500000 * u.meter / u.second


def doppler_factor_shift(wave, velocity):
    """Doppler shift using doppler factor."""
    beta = velocity / constants.c
    doppler_factor = ((1 + beta) / (1 - beta)) ** 0.5

    y = wave * doppler_factor
    return y


def doppler_shift(wave, velocity):
    """Doppler shift using non-realtivistically."""
    beta = velocity / constants.c
    y = wave * (1 + beta)
    return y


shift_factor = doppler_factor_shift(wavelength, rv)
shift = doppler_shift(wavelength, rv)

print(shift - shift_factor)
