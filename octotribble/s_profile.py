import matplotlib.pyplot as plt
import numpy as np

from RVCalculations import Calculate_WLshift


def s_profile(wave, wav_0, depth, width, k):
    """Generate an S-profile from appendix of ferluga et al. """

    s_amp = 2 * depth * np.exp(-np.pi * depth**2 * ((wave - wav_0)**2 + (k/2)**2) / width**2)
    s_inner = (np.pi * depth**2 * (wave - wav_0) * k) / width**2
    return s_amp * np.sinh(s_inner)


# I1, W1 and wav_1 are peak parameters.

wav_center = 2110.5
depth = 0.5
width = 1

rv1 = Calculate_WLshift(wav_center, 500)
rv2 = Calculate_WLshift(wav_center, 2000)
rv3 = Calculate_WLshift(wav_center, 5000)
rv4 = Calculate_WLshift(wav_center, 10000)

lamda = np.linspace(2090, 2130, 50)

s1 = s_profile(lamda, wav_center, depth, width, rv1)
s2 = s_profile(lamda, wav_center, depth, width, rv2)
s3 = s_profile(lamda, wav_center, depth, width, rv3)
s4 = s_profile(lamda, wav_center, depth, width, rv4)

plt.plot(lamda, s1, label="0.5 km/s")
plt.plot(lamda, s2, label="2 km/s")
plt.plot(lamda, s3, label="5 km/s")
plt.plot(lamda, s4, label="10 km/s")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Amplitude")
plt.legend()
plt.show()
