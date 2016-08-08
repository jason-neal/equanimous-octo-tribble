#obs_info.py

# script to extract observational info from the raw nod files.
# e.g. start time, end time, exposure time, # obs
# average airmass,

# generate some values for paper

from __future__ import print_function, division
import os
import datetime
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
    parser.add_argument("obslist", help='Input file name', type=str)
    parser.add_argument("-o","--output",help="Specify output filename", default=False)
    parser.add_argument("-s","--separator",help="Separator for time section of fits", default=":")

    args = parser.parse_args()
    return args


def main(fname, output=False, separator=":", verbose=True):

    # verbose printing on flag raise
    v_print = print if verbose else lambda *a, **k: None

    ############### Observation Settings
    ospath = os.getcwd() + "/"
    
    path = "/".join(fname.split("/")[:-1]) #+ "/"
    print("Path = {}",path)

    nod_dates = []
    nod_slits = []
    nod_exptimes = []
    nod_instruments = []
    nod_object = []
    nod_filter = []
    nod_airmass = []

    with open(fname, "r") as f:
        for line in f:

            #Replace ":"" in list file with "-"
            if separator != ":":
                line_parts = line.split(":")
                line = separator.join(line_parts)

            fitsname = line[:-1] + ".fits"
            v_print("The fits file to open is" , fitsname)
            try:
                v_print("Trying location = {}".format(path + "../Raw_files/" + fitsname))
                header = fits.getheader(path + "../Raw_files/" + fitsname)
            except:
                try:
                    v_print("Trying location = {}".format(path + "Raw_files/" + fitsname))
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
            nod_object.append(header["OBJECT"])
            nod_filter.append(header["ESO INS FILT1 NAME"])

            # AIRMASS values
            airmass_start = header["HIERARCH ESO TEL AIRM START"]
            airmass_end = header["HIERARCH ESO TEL AIRM END"]
            nod_mean_airmass = round((airmass_start + airmass_end) / 2 , 4)
            nod_airmass.append(nod_mean_airmass)
            #ESO TEL AIRM END: 1.204
            #ESO TEL AIRM START: 1.214

        
    v_print("Values obtained from the list files")
    v_print("The nod date-obs = {} ".format(nod_dates))
    v_print("The nod Slit Widths = {}".format(nod_slits))
    v_print("The nod Exposure times = {} seconds".format(nod_exptimes))
    v_print("The listed Instruments = {} ".format(nod_instruments))
        
       
    if len(set(nod_instruments)) == 1:
        instrument = nod_instruments[0]
    else:
        #raise error as list intruments are not unique
        raise NodError("Nods in list were not taken on the same instrument")
    if len(set(nod_exptimes)) == 1:
        exptime = nod_exptimes[0]
    else:
        #raise error as list exptime are not unique
        raise NodError("Nods in the list do not have same exposure times")
    if len(set(nod_slits)) == 1:
        slit_width = nod_slits[0]
    else:
        #raise error as list exptime are not unique
        raise NodError("Nods in list do not have same slit widths")
    # Get average time of observation (adding exposure time)
    if len(set(nod_object)) == 1:
        Object = nod_object[0]
    else:
        #raise error as list exptime are not unique
        raise NodError("Nods in list are not of the Object")
    
    if len(set(nod_filter)) == 1:
        nod_filter = nod_filter[0]
    else:
        #raise error as list exptime are not unique
        raise NodError("Nods in list are not of the Object")
  
    

    jd_days = []
    for nod_date in nod_dates:
        t = Time(nod_date, format="fits") + TimeDelta(exptime/2, format="sec")
        t.format= 'jd'
        jd_days.append(t.value)
    v_print("Nod dates converted into JD {}".format(jd_days))
    if max(jd_days)-min(jd_days) > len(jd_days)*2*exptime/86400.:
        raise NodError("Issue with time for the nod observations took longer " \
            "then twice exposure time for each exposure (maybe need to add a lower" \
            " limit for quick observations) ")
    
    mean_obs_time = np.mean(jd_days)
    median_obs_time = np.median(jd_days)

    obs_begin_time = Time(nod_dates[0], format="fits")
    print(" obs_begin_time",  obs_begin_time, "nod_dates[0]", nod_dates[0])
    obs_end_time = Time(nod_dates[-1], format="fits") + TimeDelta(exptime, format="sec")
    
    mean_obs_t = Time(mean_obs_time, format="jd")
    mean_obs_t.format = "fits"
    med_obs_t = Time(median_obs_time, format="jd")
    med_obs_t.format = "fits"

    obs_begin_time2 = Time(obs_begin_time)
    obs_end_time2 = Time(obs_end_time)
    obs_begin_time2.format = "jd"
    obs_end_time2.format = "jd"
    obs_total = obs_end_time2-obs_begin_time2         
    obs_total_minutes = obs_total * 24*60
   
    Middle_obs = obs_begin_time2 + TimeDelta(60*obs_total_minutes/2, format="sec")
    Middle_obs.format = "fits"
    Middle_obs2 = Time(Middle_obs)
    Middle_obs2.format = "jd"

    v_print("Mean obs time averged in JD {}".format(mean_obs_time))
    v_print("Median obs time in JD {}".format(median_obs_time))
    v_print("Begin obs time in JD {}".format(obs_begin_time))
    v_print("End obs time in JD {}".format(obs_end_time))

    
    target_ra = header["RA"]   # of the last observation in the list
    target_dec = header["DEC"]

    obs_wl_min = header["HIERARCH ESO INS WLEN MIN"]
    obs_wl_max = header["HIERARCH ESO INS WLEN MAX"] 

    # Print header values to find the ones most interested in.
    ##for key in header:
        #print(str(key) + ": " + str(header[key]))
    
    
    Obs_name = header["ESO OBS NAME"]
    
    Nabcycles = header["ESO SEQ NABCYCLES"] # ESO SEQ NABCYCLES: 4
    NEXP = header["ESO TPL NEXP"] # ESO TPL NEXP: 8
    
    # Observation IDs
    obs_pi_id= header["ESO OBS PI-COI ID"]
    prosal_id= header["ESO OBS PROG ID"]
    tpl_id = header["ESO TPL ID"]  #CRIRES_spec_obs_AutoNodOnSlit
    ao_loop_state = header["ESO AOS RTC LOOP STATE"]


    #Airmass calculation
    mean_airmass = np.mean(nod_airmass)
    airmass_range = np.max(nod_airmass)-np.min(nod_airmass)

    ####### Observatory Settings 
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

    ####### Target Settings
    ra_angle = Angle(target_ra, u.degree)
    ra_j2000 = str(ra_angle.to_string(unit=u.hour, sep=':', precision=0, pad=True))  # Extra str to get rid of u"string" which failed in template
    dec_angle = Angle(target_dec, u.deg)
    dec_j2000 = str(dec_angle.to_string(unit=u.degree, sep=':', precision=0, pad=True)) # Extra str to get rid of u"string" which failed in template
    v_print("RA degrees = {0}, RA Angle = {1}, RA for tapas = {2}".format(target_ra, ra_angle, ra_j2000))
    v_print("DEC degrees = {0}, DEC Angle = {1}, DEC for tapas = {2}".format(target_dec, dec_angle, dec_j2000))

        
    # Resolving power rule of thumb
    R = 100000*0.2 / slit_width
            
    resolving_power = int(R)
    print("Resolving Power = {}".format(resolving_power))

    



    if output:
       output_file = output
    else: 
        output_file = "key_observation_values.txt"
    try:
        with open(output_file, "w") as out:
############################################################
            #Put the useful parameters I need to know here
            
            out.write("KEY CRIRES OBSERVATION PARAMETERS:\n------------------\n")
            out.write("Target Name           = {}\n".format(Object))
            out.write("Observation Name      = {}\n".format(Obs_name))
            out.write("RA                    = {}, {} \n".format(ra_angle, ra_j2000))
            out.write("DEC                   = {}, {} \n".format(dec_angle, dec_j2000))
            
            out.write("\nOBS IDs:\n------------------\n")
            out.write("Proposal ID           = {}\n".format(prosal_id))
            out.write("PI-COI ID             = {}\n".format(obs_pi_id))
            out.write("Obs Template ID       = {}\n".format(tpl_id))
            out.write("AO State(OPEN=noAo)   = {}\n".format(ao_loop_state))

            out.write("\nDetector Setup:\n------------------\n")
            out.write("Telescope             = {}\n".format(telescope))
            out.write("Instrument            = {}\n".format(instrument))
            out.write("Spec Filter           = {}\n".format(nod_filter))
            out.write("Min wavelength        = {} nm\n".format(obs_wl_min))
            out.write("Max Wavelength        = {} nm\n".format(obs_wl_max))
            out.write("Slit width            = {} \n".format(slit_width))
            out.write("Resolving Power       = {}\n".format(resolving_power))
            out.write("Exposure time         = {}\n".format(exptime))   
            out.write("Number AB cycles      = {}\n".format(Nabcycles))
            out.write("Number of Exposures   = {}\n".format(NEXP))
                 
            out.write("\nTime Calculations:\n------------------\n")
            # FITS Format
            out.write("## FITS format\n")
            out.write("Obs Start time        = {}\n".format(obs_begin_time))
            out.write("Average obs time      = {}\n".format(mean_obs_t))
            out.write("Median obs time       = {}\n".format(med_obs_t))
            #out.write("Middle of obs         = {} \n".format(Middle_obs))
            out.write("Obs End time          = {}\n".format(obs_end_time))
            obs_begin_time.format = "jd"
            obs_end_time.format = "jd"
            # JD Format
            out.write("## JD format\n")
            out.write("Obs start time (JD)   = {}\n".format(obs_begin_time2))
            out.write("Average obs time (JD) = {}\n".format(mean_obs_time))      
            out.write("Medaian obs time (JD) = {}\n".format(median_obs_time))
            #out.write("Middle of obs (JD)    = {}\n".format(Middle_obs2))
            out.write("Obs End time (JD)     = {}\n".format(obs_end_time2))
            
            out.write("Total observation time  [inc. exptime] (JD)  = {} \n".format(obs_total))
            out.write("Total observation time  [inc. exptime] (min) = {} \n".format(obs_total_minutes))
            out.write("Total intergration time [NEXP*exptime] (min) = {} \n".format(exptime*NEXP/60))
      
            out.write("\nAirmass values:\n------------------\n")
            out.write("Nod airmasses = {}\n".format(nod_airmass))
            out.write("Mean airmass valu       = {}\n".format(mean_airmass))
            out.write("Airmass range [max-min] = {}\n".format(airmass_range))
            
            # DATE CREATED
            out.write("\n\nThis file was generated on {}\n".format(datetime.datetime.now()))
            

#############################################################
        print("Saved Observation information to \t {0}".format(output_file))
    except:
        print("Failed to save Observation information to \t {0}.".format(output_file))
        raise

if __name__ == '__main__':
    args = vars(_parser())
    fname = args.pop('obslist')
    opts = {k: args[k] for k in args}

    main(fname, **opts)