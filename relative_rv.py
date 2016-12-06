#

# relative_rv.py

# Script to calculate the relative Radial Velocity between two given wavelengths

# Jason Neal
# December 2016

def relative_rv(wav_1, wav_2):
""" Calculate the radial velocity difference between two wavelength values"""
    c = 299792.458   #  km / s
    difference = wav_2 - wav_1
    relative = difference / wav_1
    return relative * c


def parser():
    pass


def main(wave_1, wave_2, both=False):

    print(" The relative RV is ...")
    if both:
        print(" The reverse relative RV is ...")


if __name__ == "__main__":
    pass
