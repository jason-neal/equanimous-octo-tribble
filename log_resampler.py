# Logorithmic linear resampling.

import numpy as np

def log_resampler(wavelength, flux):
    """ Log linear resampling """
    # Take log of wavelength
    log_wav = np.log10(wavelength)

    # Minimum difference in log wavelength
    min_log = np.min(abs(log_wav[1:]-log_wav[:-1]))
    print("min log", min_log)
    print("log_wav[0]", log_wav[0])
    print("log_wav[-1]", log_wav[-1])

    # Create new wavelength vector
    new_wave = np.arange(log_wav[0], log_wav[-1], min_log)

    print("new wave", new_wave)

    new_flux = np.interp(new_wave, log_wav, flux)

    return new_wave, new_flux


def test_resampler():
    """ Test reampler visually"""

    import matplotlib.pyplot as plt
    from astropy.io import fits
    from Get_filenames import get_filenames

    detector = 1
    folder = "/home/jneal/Phd/data/Crires/BDs-DRACS/HD30501-1/Combined_Nods/"
    file_names = get_filenames(folder, "C*_{0}.nod.*".format(detector), "*wavecal.fits")

    print("file names", file_names)

    data = fits.getdata(file_names[0])

    resampled_data = log_resampler(data["Wavelength"], data["Extracted_dracs"])

    plt.plot(data["Wavelength"], data["Extracted_dracs"], label="org")
    plt.plot(10 ** resampled_data[0], resampled_data[1], "o-", label="resampled")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    test_resampler()
