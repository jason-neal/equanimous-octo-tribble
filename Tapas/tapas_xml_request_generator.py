#!/usr/bin/env python3
# -*- coding: utf8 -*-
"""Take CRIRES fits file and generate a tapas xml request to submit on.

"""
from __future__ import print_function
import os
import argparse
import numpy as np
from astropy.time import Time, TimeDelta
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import Angle


def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Xml generator for tapas submission')
    parser.add_argument("fname", help='Input file name', type=str)
    parser.add_argument("-o", "--output_file", help='Output file', type=str, default="tapas_request.xml")
    parser.add_argument("-l", "--listspectra", action="store_true",
                        help="Was filename a DRACS list of nod spectra for the observation (without fits extensions)")
    parser.add_argument("-r", "--resolvpower", default=False,
                        help="Specify Instrument resolution power, defined as λ/FWHM for convolution")
    parser.add_argument("-s", "--sampling", type=int, default=10,
                        help=("Sampling ratio - This is the number of points per FHWM interval on which the"
                              "convolved transmission will be sampled."))
    parser.add_argument("-f", "--instrument_function", help="Instrument function - gaussian or none",
                        type=str, default="gaussian", choices=["none", "gaussian"])
    parser.add_argument("-u", "--unit", help="Spectra Unit", choices=["air", "vacuum", "wavenumber"], type=str,
                        default="vacuum")
    parser.add_argument("-b", "--berv", help="Have BERV RV correction applied to the Tapas spectra",
                        action="store_true")
    parser.add_argument("-t", "--tapas_format", help="Tapas file format", type=str, default="ascii",
                        choices=["ascii", "fits", "netcdf", "vo"])
    parser.add_argument("-c", "--constituents", help="Atmospheric constituents for spectra", type=str, default="all",
                        choices=["all", "ray", "h2o", "o2", "o3", "co2", "ch4", "n2o", "not_h2o"])
    parser.add_argument("-i", "--request_id", type=int, default=0,
                        help=("Request ID number, I use want to change from species identification."
                              " E.g. 10=all, 11=ray, 12=h20 etc."))
    parser.add_argument("--wl_min", help="Minimum Wavelength", default=False)
    parser.add_argument("--wl_max", help="Maximum Wavelength", default=False)
    parser.add_argument("-n", "--request_number", help="Tapas request number. Iterate on previous request number.",
                        type=int)
    parser.add_argument("-v", "--verbose", action="store_true", help="Turn on verbosity.")
    args = parser.parse_args()
    return args


def main(fname=("/home/jneal/Phd/data/Crires/BDs-DRACS/HD30501-1/Combined_Nods/"
                "CRIRE.2012-04-07T00:08:29.976_1.nod.ms.norm.sum.wavecal.fits"),
         listspectra=False, resolvpower=False, unit="vacuum", instrument_function="gaussian",
         sampling=10, berv=False, tapas_format="ASCII", constituents="all", output_file=False,
         request_id=10, wl_min=False, wl_max=False, request_number=None, verbose=False):
    """Fill out tapas template."""
    # verbose printing on flag raise
    v_print = print if verbose else lambda *a, **k: None

    # Dotfile for saving requst number for iteration
    home = os.path.expanduser("~")
    dot_file = os.path.join(home, ".tapas_request_num")

    # ############## Observation Settings
    ospath = os.getcwd() + "/"
    v_print("OS path {0}".format(ospath))

    path = "/".join(fname.split("/")[:-1]) + "/"
    v_print("Path obtained from fname = {0}".format(path))

    if listspectra:
        nod_dates = []
        nod_slits = []
        nod_exptimes = []
        nod_instruments = []
        with open(fname, "r") as f:
            for line in f:
                fitsname = line[:-1] + ".fits"
                v_print("The fits file to open is", fitsname)
                try:
                    v_print("Trying location{0}".format(path + "../Raw_files/" + fitsname))
                    header = fits.getheader(path + "../Raw_files/" + fitsname)
                except:
                    try:
                        v_print("Trying location  = {0}".format(path + "Raw_files/" + fitsname))
                        header = fits.getheader(path + "Raw_files/" + fitsname)
                    except:
                        try:
                            header = fits.getheader(path + fitsname)
                        except:
                            v_print("Having issues finding files from list.")
                            raise

                obsdate = header["DATE-OBS"]
                nod_dates.append(obsdate)
                nod_slits.append(header["HIERARCH ESO INS SLIT1 WID"])
                nod_exptimes.append(header["EXPTIME"])
                nod_instruments.append(header["INSTRUME"])
            v_print("Values obtained from the list files")
            v_print("The nod date-obs = {0} ".format(nod_dates))
            v_print("The nod Slit Widths = {0}".format(nod_slits))
            v_print("The nod Exposure times = {0} seconds".format(nod_exptimes))
            v_print("The listed Instruments = {0} ".format(nod_instruments))

            if len(set(nod_instruments)) == 1:
                instrument = nod_instruments[0]
            else:
                # raise error as list intruments are not unique
                raise NodError("Nods in list were not taken on the same instrument")
            if len(set(nod_exptimes)) == 1:
                exptime = nod_exptimes[0]
            else:
                # raise error as list exptime are not unique
                raise NodError("Nods in the list do not have same exposure times")
            if len(set(nod_slits)) == 1:
                slit_width = nod_slits[0]
            else:
                # raise error as list exptime are not unique
                raise NodError("Nods in list do not have same slit widths")
            # Get average time of observation (adding exposure time)
            jd_days = []
            for nod_date in nod_dates:
                t = Time(nod_date, format="fits") + TimeDelta(exptime / 2, format="sec")
                t.format = 'jd'
                jd_days.append(t.value)
            v_print("Nod dates converted into JD {0}".format(jd_days))
            if max(jd_days) - min(jd_days) > len(jd_days) * 2 * exptime / 86400.:
                raise NodError("Issue with time for the nod observations took longer then twice exposure time for each"
                               " exposure (maybe need to add a lower limit for quick observations) ")
            mean_obs_time = np.mean(jd_days)
            v_print("Mean obs time averged in JD {0}".format(mean_obs_time))
            v_print("Median obs time in JD {0}".format(np.median(jd_days)))
            date = Time(mean_obs_time, format="jd")
            date.format = "fits"  # Turn back into fits format
            date = date.value[:-5] + "Z"    # Going from Time object back to string and adding Z to end for Tapas
            target_ra = header["RA"]   # of the last observation in the list
            target_dec = header["DEC"]
            obs_wl_min = header["HIERARCH ESO INS WLEN MIN"]
            obs_wl_max = header["HIERARCH ESO INS WLEN MAX"]
    else:
        header = fits.getheader(fname)

        date = header["DATE-OBS"][:-1] + "Z"
        target_ra = header["RA"]
        target_dec = header["DEC"]
        obs_wl_min = header["HIERARCH ESO INS WLEN MIN"]
        obs_wl_max = header["HIERARCH ESO INS WLEN MAX"]
        slit_width = header["HIERARCH ESO INS SLIT1 WID"]

    # ###### Observatory Settings
    instrument = header["INSTRUME"]
    telescope = header["TELESCOP"]
    if "VLT" in telescope:
        obs_name = "ESO Paranal Chile (VLT)"
    else:
        v_print("Don't have TAPAS telecope location name (Using the Obervation value)")
        obs_name = telescope

    obs_long = header["HIERARCH ESO TEL GEOLON"]
    obs_lat = header["HIERARCH ESO TEL GEOLAT"]
    obs_alt = header["HIERARCH ESO TEL GEOELEV"]

    # ###### Target Settings
    ra_angle = Angle(target_ra, u.degree)
    # Extra str to get rid of u"string" which failed in template
    ra_j2000 = str(ra_angle.to_string(unit=u.hour, sep=':', precision=0, pad=True))
    dec_angle = Angle(target_dec, u.deg)
    # Extra str to get rid of u"string" which failed in template
    dec_j2000 = str(dec_angle.to_string(unit=u.degree, sep=':', precision=0, pad=True))
    v_print("RA degrees = {0}, RA Angle = {1}, RA for tapas = {2}".format(target_ra, ra_angle, ra_j2000))
    v_print("DEC degrees = {0}, DEC Angle = {1}, DEC for tapas = {2}".format(target_dec, dec_angle, dec_j2000))

    if wl_min:
        spec_range_min = int(wl_min)
    else:
        spec_range_min = int(obs_wl_min - 10)
    if wl_max:
        spec_range_max = int(wl_max)
    else:
        spec_range_max = int(obs_wl_max + 10)

    # Check for TAPAS consistency
    # Tapas wavelength range is from 350 to 2500 nm in vacuum.
    if spec_range_min < 350:
        v_print(("Lower wavelength bound of {0} was below tapas minimum."
                " Setting to tapas minimum of 350 nm").format(spec_range_min))
        spec_range_min = 350
    elif spec_range_min > 2500:
        v_print(("Wavelength Lower bound of {0} is above Tapas maximum wavelength of 2500 nm."
                "Check your wavelength units").format(spec_range_min))
        raise("Wavelength bounds Error")
    if spec_range_max > 2500:
        v_print(("Upper wavelength bound of {0} was above tapas maximum."
                " Setting to tapas maximum of 2500 nm").format(spec_range_max))
        spec_range_max = 2500
    elif spec_range_max < 350:
        v_print(("Wavelength Upper bound of {0} is below Tapas minimum wavelength of 350 nm."
                 " Check your wavelength units").format(spec_range_max))
        raise("Wavelength bounds Error")

    # ########## Tapas Specifications
    if tapas_format in ["ascii", "fits", "vetcdf", "vo"]:
        tapas_format = tapas_format.upper()
    else:
        v_print("TAPAS format was not correct, using default of ASCII")
        tapas_format = "ASCII"

# open save file, find request id, add 1
# get request number from previous xml file if not given.

    if not request_number:
        # Load previous tapas request number
        try:
            with open(dot_file, "r") as req:
                old_number = req.readline().split("=")[1]
                v_print("Previous request_number = {0}".format(old_number))
                request_number = int(old_number) + 1
                v_print("New request_number number = {0}".format(request_number))
        except:
            request_number = 0
            print(("Could not read request number from previous request in {0}. The default value of {1} will need"
                   " to be manually changed when submitting the request.").format(dot_file, request_number))
    else:
        # Manually given request number
        v_print("Manually given request number={0} will be used.\n".format(request_number))

    # #### Specify Species in tapas spectra
    # Could possibly make this a binary AND operation
    if constituents == "all":
        species_map = [1, 1, 1, 1, 1, 1, 1]
        species_id = 10
    elif constituents == "ray":
        species_map = [1, 0, 0, 0, 0, 0, 0]
        species_id = 11
    elif constituents == "h2o":
        species_map = [0, 1, 0, 0, 0, 0, 0]
        species_id = 12
    elif constituents == "o3":
        species_map = [0, 0, 1, 0, 0, 0, 0]
        species_id = 13
    elif constituents == "o2":
        species_map = [0, 0, 0, 1, 0, 0, 0]
        species_id = 14
    elif constituents == "co2":
        species_map = [0, 0, 0, 0, 1, 0, 0]
        species_id = 15
    elif constituents == "ch4":
        species_map = [0, 0, 0, 0, 0, 1, 0]
        species_id = 16
    elif constituents == "n2o":
        species_map = [0, 0, 0, 0, 0, 0, 1]
        species_id = 17
    elif constituents == "not_h2o":
        species_map = [1, 0, 1, 1, 1, 1, 1]
        species_id = 18
    else:
        species_id = 9999
        raise SpeciesError("Choice of constituents was not in valid range.")

    # Request id baised on the species present unless manually specified
    if request_id:
        Request_ID = request_id
    else:
        Request_ID = species_id

    if species_map[0]:
        ray = "YES"
    else:
        ray = "NO"
    if species_map[1]:
        h20 = "YES"   # alternate yes no for h20 for scaling
    else:
        h20 = "NO"   # alternate yes no for h20 for scaling
    if species_map[2]:
        o2 = "YES"
    else:
        o2 = "NO"
    if species_map[3]:
        o3 = "YES"
    else:
        o3 = "NO"
    if species_map[4]:
        co2 = "YES"
    else:
        co2 = "NO"
    if species_map[5]:
        ch4 = "YES"
    else:
        ch4 = "NO"
    if species_map[6]:
        n2o = "YES"
    else:
        n2o = "NO"

    if "air" in unit.lower():
        spectral_choice = "Standard Wavelength (nm)"
    if "vac" in unit.lower():
        spectral_choice = "Vacuum Wavelength (nm)"
    if "wave" in unit.lower():
        spectral_choice = "Wavenumber (cm-1)"

    instrument_function = instrument_function.lower()
    if "none" in instrument_function:
        ilsf_choice = -1    # -1, 0, 1
    elif "gaussian" in instrument_function:
        ilsf_choice = 1    # -1, 0, 1
    else:
        v_print("Instrument function not specifid correctly\n Valid choices are none or gaussian, default is gaussian.")

    sampling_ratio = sampling  # defualt 10

    if resolvpower:
        resolving_power = int(resolvpower)
        v_print("Resolving power manually specified at {0}".format(resolving_power))
    else:
        if "CRIRES" in instrument:
            v_print("Resolving Power\nUsing the rule of thumb equation from the CRIRES manual. \nWarning!"
                    "The use of adpative optics is not checked for!!")
            R = 100000 * 0.2 / slit_width

            resolving_power = int(R)
            v_print("Slit width was {0} inches.\nTherefore the resolving_power is set = {1}".format(slit_width,
                                                                                                    resolving_power))

            # TODO: Specify other instruments in here?

        else:
            v_print("Resolving power not defined")

    if berv:
        apply_berv = "YES"
    else:
        apply_berv = "NO"

    # ##### Atmosphere Parameters
    # hours(date[11:13]) only seem to work if multiple of 06 hours
    # assuming 00 is for 00->06 hours
    hour = int(date[11:13])
    if 0 <= hour < 6:
        hour_out = "00"
    elif 6 <= hour < 12:
        hour_out = "06"
    elif 12 <= hour < 18:
        hour_out = "12"
    elif 18 <= hour < 24:
        hour_out = "18"
    else:
        raise HourError("Error with the arletty timing. The request will fail")

    file_date = date[0:4] + date[5:7] + date[8:10] + hour_out
    arletty_file = "canr_" + file_date + ".arl"
    ecmwf_file = "canr_" + file_date + "_qo3.txt"

    v_print("arletty_file", arletty_file)
    v_print("ecmwf_file", ecmwf_file)

    d = {"request_number": request_number, "Request_ID": Request_ID, "tapas_format": tapas_format, "ray": ray,
         "h20": h20, "o3": o3, "o2": o2, "co2": co2, "ch4": ch4, "n2o": n2o, "date": date, "obs_name": obs_name,
         "obs_long": obs_long, "obs_lat": obs_lat, "obs_alt": obs_alt, "ra_j2000": ra_j2000, "dec_j2000": dec_j2000,
         "spectral_choice": spectral_choice, "spec_range_min": spec_range_min, "spec_range_max": spec_range_max,
         "ilsf_choice": ilsf_choice, "resolving_power": resolving_power, "sampling_ratio": sampling_ratio,
         "apply_berv": apply_berv, "arletty_file": arletty_file, "ecmwf_file": ecmwf_file}

    template = """<?xml version="1.0" encoding="UTF-8"?>
<tapas Id="Ether_TAPAS_$request_number">
<request Id="$Request_ID">
<preferences>
<format valid="VO,ASCII,FITS,NETCDF">$tapas_format</format>
<rayleigh_extinction valid="YES,NO">$ray</rayleigh_extinction>
<h2o_extinction valid="YES,NO">$h20</h2o_extinction>
<o3_extinction valid="YES,NO">$o3</o3_extinction>
<o2_extinction valid="YES,NO">$o2</o2_extinction>
<co2_extinction valid="YES,NO">$co2</co2_extinction>
<ch4_extinction valid="YES,NO">$ch4</ch4_extinction>
<n2o_extinction valid="YES,NO">$n2o</n2o_extinction>
</preferences>
<observation>
<date>$date</date>
<observatory>
<name>$obs_name</name>
<longitude min="-180" max="180">$obs_long</longitude>
<latitude min="-90" max="90">$obs_lat</latitude>
<altitude min="0" max="10000">$obs_alt</altitude>
</observatory>
<los>
<ra_j2000>$ra_j2000</ra_j2000>
<dec_j2000>$dec_j2000</dec_j2000>
</los>
<instrument>
<spectral_choice valid="Vacuum Wavelength (nm),Standard Wavelength (nm),Wavenumber (cm-1)">$spectral_choice</spectral_choice>
<spectral_range>$spec_range_min $spec_range_max</spectral_range>
<ilsf_choice valid="-1,0,1">$ilsf_choice</ilsf_choice>
<resolving_power min="0">$resolving_power</resolving_power>
<sampling_ratio min="0">$sampling_ratio</sampling_ratio>
<appli_berv valid="YES,NO">$apply_berv</appli_berv>
</instrument>
</observation>
<atmosphere>
<reference valid="0,1,2,3,4,5,6">0</reference>
<arletty_file>/data/tapas///arletty/$arletty_file</arletty_file>
<ecmwf_file>/data/tapas///ecmwf/$ecmwf_file</ecmwf_file>
</atmosphere>
</request>
</tapas> """

# ECMWF = ARLETTY
# [{"text": "ECMWF", "value": 0}, {"text": "TROPICAL", "value": 1}, {"text": "MEDIUM_LATITUDE_SUMMER", "value": 2},
# {"text": "MEDIUM_LATITUDE_WINTER", "value": 3}, {"text": "SUBARCTIC_SUMMER", "value": 4}, {"text": "SUBARCTIC_WINTER",
# "value": 5}, {"text": "US_STANDARD_1976", "value": 6}]

# ###################### TODO - Missing all in one transmission

# tapasTexts["NM_STANDARD"] = "Standard Wavelength (nm)";
#            tapasTexts["NM_VACUUM"] = "Vacuum Wavelength (nm)";
#            tapasTexts["CM"] = "Wavenumber (cm-1)";
#
#            tapasTexts["NONE"] = "None";
#            tapasTexts["GAUSSIAN"] = "Gaussian";
#            tapasTexts["ECMWF"] = "ARLETTY";
#            tapasTexts["TROPICAL"] = "Tropical";
#            tapasTexts["MEDIUM_LATITUDE_SUMMER"] = "Average latitude summer";
#            tapasTexts["MEDIUM_LATITUDE_WINTER"] = "Average latitude winter";
#            tapasTexts["SUBARCTIC_SUMMER"] = "Subarctic summer";
#            tapasTexts["SUBARCTIC_WINTER"] = "Subarctic winter";
#            tapasTexts["US_STANDARD_1976"] = "Standard US 1976";

    from string import Template
    s = Template(template)

    v_print("Xml Template object = {0}".format(s))

    sub = s.substitute(d)
    # print(sub)

    # Save xml ouptut for copying to tapas
    # TO TRY in future - submit straight to tapas

    # output_file = "/home/jneal/Phd/Codes/UsefulModules/Tapas_xml_request_file.xml"

    if not output_file:
        output_file = "tapas_request_{0}.xml".format(request_number)
    else:
        # Add request number
        split = output_file.split(".")
        if len(split) == 2:
            output_file = split[0] + "_{0}.".format(request_number) + split[1]
        else:
            print("Warning output filename may not be the expected format 'some_name.xml'")
        print("Split output file", split)

    try:
        with open(output_file, "w") as out:
            out.write(sub)
        print("Saved tapas xml request to \t {0}".format(output_file))
    except:
        print("Failed to save xml request to \t {0}. \nHere is a printed version.".format(output_file))
        print("\n{0}\n".format(sub))

    try:
        with open(dot_file, "w") as req:
            req.write("request_number = {0}".format(request_number))
        v_print("Stored current request number {0} in {1} for later iteration".format(request_number, dot_file))
    except:
        v_print("Failed to store request number {0} in {1} .".format(request_number, dot_file))


class NodError(Exception):
    """A Error class for nod errors."""

    pass


class SpeciesError(Exception):
    """A Error class for species errors."""

    pass


class HourError(Exception):
    """A Error class for hour errors."""

    pass


if __name__ == "__main__":
    args = vars(_parser())
    fname = args.pop('fname')
    opts = {k: args[k] for k in args}

    main(fname, **opts)
