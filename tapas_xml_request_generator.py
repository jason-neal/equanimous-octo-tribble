#!/usr/bin/env python3
# -*- coding: utf8 -*-

""" Take CRIRES fits file and generate a tapas xml request to submit on tapas

"""
import argparse
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import Angle

def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Xml generator for tapas submission')
    parser.add_argument("fname", help='Input file name', type=str)
    parser.add_argument("-o", "--output_file", help='Output file', type=str)
    parser.add_argument("-l", "--listspectra", help="Was filename a list of spectra for observation", action="store_true")
    parser.add_argument("-r", "--resolvpower", help="Specify Instrument resolution power, defined as λ/FWHM for convolution", default=False) 
    parser.add_argument("-s", "--sampling", type=int, help="Sampling ratio - This is the number of points per FHWM interval on which the convolved transmission will be sampled.", default=10)     
    parser.add_argument("-i", "--instrument_function", help="Instrument function - gaussian or none", type=str, default="gaussian", choices=["none","gaussian"])     
    parser.add_argument("-u", "--unit", help="Spectra Unit", choices=["air", "vacuum", "wavenumber"], type=str, default="air")
    parser.add_argument("-b", "--berv", help="Have BERV RV correction applied to the Tapas spectra", action="store_true") 
    parser.add_argument("-f", "--tapas_format", help="Tapas file format", type=str, default="ascii", choices=["fits","ascii","netcdf","vo"])     
    parser.add_argument("-c","--constituents", help="Atmospheric constituents for spectra", type=str, default="all", choices=[ "all","ray", "h2o", "o2", "o3", "co2", "ch4", "n2o", "not_h2o"])    
    parser.add_argument("--request_id", help="Request ID number", type=int, default="10000")    
    parser.add_argument("--wl_min", help="Minimum Wavelength", default=False)    
    parser.add_argument("--wl_max", help="Maximum Wavelength", default=False)    
    parser.add_argument("-n", "--request_number", help="Request number, iterate this", default=0, type=int)    

    args = parser.parse_args()
    return args


def main(fname="/home/jneal/Phd/data/Crires/BDs-DRACS/HD30501-1/Combined_Nods/CRIRE.2012-04-07T00:08:29.976_1.nod.ms.norm.sum.wavecal.fits", 
        listspectra=False, resolvpower=False, unit="air", instrument_function="gaussian",
        sampling=10, berv=False, tapas_format="ASCII", constituents="all", output_file=False, request_id=10000, wl_min=False, wl_max=False,request_number=0):

    output_file = "/home/jneal/Phd/Codes/UsefulModules/Tapas_xml_request_file.xml"  # For testing
    

    ############### Observation Settings
    if listspectra:
        pass
        # Calcualte averages from loading all headers
    #MJD-OBS =       56024.00590250 / Obs start 2012-04-07T00:08:29.976
    #DATE-OBS= '2012-04-07T00:08:29.9764' / Observing date
    #EXPTIME =          180.0000000 / Integration time
        #target_ra = header["RA"]
        #target_dec = header["DEC"]
        header = fits.getheader(firstlist_name)  # for telescope info later
    else:
        header = fits.getheader(fname)
    
        date = header["DATE-OBS"][:-1]+"Z"     
        target_ra = header["RA"]
        target_dec = header["DEC"]
        obs_wl_min = header["HIERARCH ESO INS WLEN MIN"]
        obs_wl_max = header["HIERARCH ESO INS WLEN MAX"] 
        slit_width = header["HIERARCH ESO INS SLIT1 WID"]


    ####### Observatory Settings 
    instrument = header["INSTRUME"]
    telescope = header["TELESCOP"]
    if "VLT" in telescope:
    	obs_name = "ESO Paranal Chile (VLT)"
    else:
    	print("Don't have TAPAS telecope location name (Using the Obervation value)")
    	obs_name = telescope

    obs_long = header["HIERARCH ESO TEL GEOLON"]
    obs_lat = header["HIERARCH ESO TEL GEOLAT"]
    obs_alt = header["HIERARCH ESO TEL GEOELEV"]

    ####### Target Settings
    ra_angle = Angle(target_ra, u.degree)
    ra_j2000 = str(ra_angle.to_string(unit=u.hour, sep=':',precision=0, pad=True))  # Extra str to get rid of u"string" which failed in template
    dec_angle = Angle(target_dec, u.deg)
    dec_j2000 = str(dec_angle.to_string(unit=u.degree, sep=':', precision=0, pad=True)) # Extra str to get rid of u"string" which failed in template
    print("RA Angle =",ra_angle, "RA deg =", target_ra, "ra_2000 h:m:s=", ra_j2000)
    print("Dec Angle =",dec_angle, "DEC deg =", target_dec, "dec_2000 d:m:s=", dec_j2000)

    #spectral_range = "{0} {1}".format(min_range, max_range)
    if wl_min:
        spec_range_min = int(wl_min)
    else:
        spec_range_min = int(obs_wl_min - 10)
    if wl_max:
        spec_range_max = int(wl_max)
    else:
        spec_range_max = int(obs_wl_max + 10) 
    # Check for TAPAS consistency    
    #Tapas wavelength range is from 350 to 2500 nm in vacuum.
    if spec_range_min < 350:
        print("Lower wavelength bound of {} was below tapas minimum. Setting to tapas minimum of 350 nm".format(spec_range_min))
        spec_range_min = 350
    elif spec_range_min > 2500:
        print("Wavelength Lower bound of {} is above Tapas maximum wavelength of 2500 nm. Check your wavelength units".format(spec_range_min))
        raise("Wavelength bounds Error")
    if spec_range_max > 2500:
        print("Upper wavelength bound of {} was above tapas maximum. Setting to tapas maximum of 2500 nm".format(spec_range_max))
        spec_range_max = 2500
    elif spec_range_max < 350:
        print("Wavelength Upper bound of {} is below Tapas minimum wavelength of 350 nm. Check your wavelength units".format(spec_range_max))
        raise("Wavelength bounds Error")

    ########### Tapas Specifications
    if tapas_format in ["ASCII","FITS","NETCDF","VO"]:
          tapas_format = tapas_format
    else:
        print("TAPAS format was not correct, using default of ASCII")
        tapas_format = "ASCII"

# open save file, find request id, add 1
# get request number from previous xml file if not given.
    if not request_number:
        # get from prvious
        #request_number = ######
        pass
    #try:
    #    with open(output_file, "r") as f:
     #       for line in f:
    #            if line.startswith("<request_Id"):
    #                ID_num = line.split('"')[1]
    #                print("ID_number from last file", ID_num)
    #                Request_ID = ID_num + 1  # change request_id from default if it is suposed to iterate
    #                break
   #except:
    Request_ID = request_id    # Iterate on previous ID if found

# Specify different tapas spectra
    #Could possibly make this a 6bit AND operation
    if constituents == "all":
        species_map = [1,1,1,1,1,1,1]
    elif constituents == "ray":
        species_map = [1,0,0,0,0,0,0]
    elif constituents == "h2o":
        species_map = [0,1,0,0,0,0,0]
    elif constituents == "o3":
        species_map = [0,0,1,0,0,0,0]
    elif constituents == "o2":
        species_map = [0,0,0,1,0,0,0]
    elif constituents == "co2":
        species_map = [0,0,0,0,1,0,0]
    elif constituents == "ch4":
        species_map = [0,0,0,0,0,1,0]
    elif constituents == "n2o":
        species_map = [0,0,0,0,0,0,1]
    elif constituents == "not_h2o":
        species_map = [1,0,1,1,1,1,1]
    
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
        print("Instrument function not specifid correctly\n Valid choices are none or gaussian, The default is gaussian.")

    sampling_ratio = sampling # defualt 10

    if resolvpower:
       resolving_power = resolvpower
    else:
        if "CRIRES" in instrument:
            if slit_width == 0.400:
                print("Slit width was 0.4, \nSetting resolving_power = 50000")
                resolving_power = 50000    # if Crires take the two slit width/Resolution values
            elif slit_width == 0.200:
                print("Slit width was 0.2, \nSetting resolving_power = 100000")
                resolving_power = 100000    # if Crires take the two slit width/Resolution values
            else:
                print("Slitwidth of CRIRES does not match the two fixed settings")
                print("Is this older CRIRES data?")
        else:
            print("Resolving power not defined")

    if berv:
        apply_berv = "YES"
    else:
        apply_berv = "NO"


    ###### Atmosphere Parameters
    arletty_file = "canr_"+ date[0:4] + date[5:7] + date[8:10] + date[11:13] + ".arl"
    ecmwf_file = "canr_" + date[0:4] + date[5:7] + date[8:10] + date[11:13] + "_qo3.txt"
    
    print("arletty_file", arletty_file)
    print("ecmwf_file", ecmwf_file)

    d = {"request_number":request_number, "Request_ID":Request_ID, "tapas_format":tapas_format, "ray":ray, "h20":h20, "o3":o3, "o2":o2, "co2":co2, "ch4":ch4, "n2o":n2o, "date":date, "obs_name":obs_name, "obs_long":obs_long, "obs_lat":obs_lat, "obs_alt":obs_alt, "ra_j2000":ra_j2000, "dec_j2000":dec_j2000, "spectral_choice":spectral_choice, "spec_range_min":spec_range_min, "spec_range_max":spec_range_max, "ilsf_choice":ilsf_choice, "resolving_power":resolving_power, "sampling_ratio":sampling_ratio, "apply_berv":apply_berv, "arletty_file":arletty_file, "ecmwf_file":ecmwf_file}

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
#[{"text":"ECMWF","value":0},{"text":"TROPICAL","value":1},{"text":"MEDIUM_LATITUDE_SUMMER","value":2},{"text":"MEDIUM_LATITUDE_WINTER","value":3},{"text":"SUBARCTIC_SUMMER","value":4},{"text":"SUBARCTIC_WINTER","value":5},{"text":"US_STANDARD_1976","value":6}]

#######################TODO - Missing all in one transmission

#tapasTexts["NM_STANDARD"] = "Standard Wavelength (nm)";
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

    print(s)
    
    sub = s.substitute(d)
    print(sub)

    # Save xml ouptut for copying to tapas 
    #TOTRY in future - submit straight to tapas

    #output_file = "/home/jneal/Phd/data/Tapas/Tapas_xml_request_file.xml"
    output_file = "/home/jneal/Phd/Codes/UsefulModules/Tapas_xml_request_file.xml"
    with open(output_file, "w") as out:
    	out.write(sub)
    print("Saved tapas xml request to copy at \t {}".format(output_file))



if __name__ == "__main__":
    args = vars(_parser())
    fname = args.pop('fname')
    opts = {k: args[k] for k in args}
    
#    fname="/home/jneal/Phd/data/Crires/BDs-DRACS/HD30501-1/Combined_Nods/CRIRE.2012-04-07T00:08:29.976_1.nod.ms.norm.sum.wavecal.fits"  # for testing
    print(opts)

    main(fname, **opts)
    #main()
