#!/usr/bin/env python

"""RVCalculations.

Calulate the Radial Velocity shift that would be
present with a given wavelength shift in the NIR range
we have.
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


def Calculate_RV(Lambda, deltalambda):
    c = 299792458.  # m/s
    V = [c * deltalambda / lmda for lmda in Lambda]
    return V


def Calculate_WLshift(wav, v):
    """Lambda = wavelength.

    v = velocity in m/s.
    """
    c = 299792458.  # m/s
    if isinstance(v, (int, float)):
        return (v / c) * wav
    else:
        deltalambda = [(v / c) * lmda for lmda in wav]
        return deltalambda


def Calc_RV(Lambda, deltalambda):
    # Fixed wavelength
    c = 299792458.  # m/s
    V = [c * (delta / Lambda) for delta in deltalambda]
    return V


def Calc_WLshift(Lambda, vels):
    # Fixed wavelength
    c = 299792458.  # m/s
    V = [(vel / c) * Lambda for vel in vels]
    return V
# NIR range we have
# 2110 nm to 2170 nm


nir_wls = [806, 900, 1020, 1220, 1630, 2190, 3450]  # NIR wavelenths in nm
nir_band = ["I", "Z", "Y", "J", "H", "K", "L"]
vels = [0.1, 0.5, 1, 5, 10, 20, 100]  # m/s velocity
deltalambdas = [0.001, 0.01, 0.05, 0.1, 0.12]  # nm
print("nir_wls", nir_wls)
print("vels", vels)
print("deltalambdas", deltalambdas)
print(" length nir_wls", len(nir_wls))
print(" length vels", len(vels))
print(" length deltalambdas", len(deltalambdas))

plt.figure()
for wl, band in zip(nir_wls, nir_band):
    delta = Calc_WLshift(wl, vels)
    plt.plot(vels, delta, label=(band + " Band/" + str(wl) + " (nm)"))

plt.legend(loc="best")
plt.title("NIR Wavelenghts Doppler effect")
plt.ylabel("Delta lambda (nm)")
plt.xlabel("Velocity (m/s)")

plt.figure()
for wl, band in zip(nir_wls, nir_band):
    RV = Calc_RV(wl, deltalambdas)
    plt.plot(deltalambdas, RV, label=(band + " Band/" + str(wl) + " (nm)"))

plt.legend(loc="best")
plt.title("NIR Wavelenghts Doppler effect")
plt.xlabel("Delta lambda (nm)")
plt.ylabel("Velocity (m/s)")

plt.show()

# ###### OLD Way using range of Lambda #######
# plt.figure()
# for vel in vels:
#     print("This vel", vel)
#     delta = Calculate_WLshift(nir_wls, vel)
#     print("This delta", delta)
#     print(" length delta", len(delta))

#     plt.plot(nir_wls, delta, label="V= "+str(vel)+" (m/s)")

# plt.legend()
# plt.title("Wavelength Differences caused by different Velocities")
# plt.ylabel("Delta lambda (nm)")
# plt.xlabel("NIR Wavelengths (nm)")

# plt.figure()
# for deltawl in deltalambdas:
#     RV = Calculate_RV(nir_wls, deltawl)
#     plt.plot(nir_wls, RV, label="deltaWL= "+str(deltawl)+" (nm)")

# plt.legend()
# plt.title("RV signal caused by different Wavelength shifts")
# plt.ylabel("RV (m/s)")
# plt.xlabel("NIR Wavelengths (nm)")

# plt.show()
