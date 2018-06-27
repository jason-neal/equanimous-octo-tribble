# Calculate wiens displacement law

# cli usage for temp = 5000 K 
# $ python wien_law.py 5000

import sys
import argparse


def wien_displacement(temperature):
    """Calcualte wien displacement law
    
    Parameters
    ----------
    temperature: float
        Temperature in K

    Returns
    -------
    wavelength: float
        Peak wavelength in microns

    Notes
    -----
    Formally, Wien's displacement law states that the
    spectral radiance of black body radiation per unit
    wavelength, peaks at the wavelength λmax given by
        
        lambda = b/T
    
    where T is the absolute temperature in kelvins. b
    is a constant of proportionality called Wien's
    displacement constant, equal to
    2.8977729(17)×10−3 m⋅K[1]
    """
    b = 2897.729  # $\mu m$.K
    return b / temperature


def inverse_wien_displacement(peak_wavelength):
    """Calcualte inverse wien displacement law
    
    Parameters
    ----------
    peak_wavelength: float
        Wavelenth of peak in micon

    Returns
    -------
    temperature: float
        Temperature of black body

    Notes
    -----
    Formally, Wien's displacement law states that the
    spectral radiance of black body radiation per unit
    wavelength, peaks at the wavelength λmax given by
        
        lambda = b/T

    where T is the absolute temperature in kelvins. b
    is a constant of proportionality called Wien's
    displacement constant, equal to
    2.8977729(17)×10−3 m⋅K[1]
    """
    b = 2897.729  # $\mu m$.K

    return b / peak_wavelength


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Caclulate Wien displacement law")
    ap.add_argument("value", help="Input temperature/wavelength")
    ap.add_argument("-i","--inverse", help="Inverse displacement. Calculate temperature from wavelength", action="store_true")
    args = vars(ap.parse_args())
    
    if args["inverse"]:
        wave = float(args["value"])
        temp = inverse_wien_displacement(wave)
        print("Temperature corresponding to a peak wavelength of {} micron is {} K".format(wave, temp))
        
    else:
        temp = float(args["value"])
        wave = wien_displacement(temp)
        print("Peak wavelength for {}K is {} micron".format(temp, wave))
    
