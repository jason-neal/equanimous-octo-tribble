#!/usr/bin/env python3
# -*- coding: utf8 -*-


""" Take CRIRES fits file and generate a tapas xml request to submit on tapas


"""

from astropy.io import fits

def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Telluric Removal')
    parser.add_argument('fname', help='Input fits file')
    
    # name  with a flag to specify if the name is to a listspectra.txt file to do whole observation!


# Many flags for tapas settings  
# Try get most from fits files
    #parser.add_argument('-x', '--export', default=False,
 #                       help='Export result to fits file True/False')
    #parser.add_argument('-o', '--output', default=False,
#                        help='Ouput Filename')
    #parser.add_argument('-t', '--tellpath', default=False,
 #                       help='Path to find the telluric spectra to use.')
    #parser.add_argument('-k', '--kind', default="linear",
 #                       help='Interpolation order, linear, quadratic or cubic')
    #parser.add_argument('-m', '--method', default="scipy",
 #                       help='Interpolation method numpy or scipy')
    #parser.add_argument("-s", "--show", default=True,
  #                      help="Show plots") #Does not wokwithout display though for some reason
    args = parser.parse_args()
    return args


def main(fname="/home/jneal/Phd/data/Crires/BDs-DRACS/HD30501-1/Combined_Nods/CRIRE.2012-04-07T00:08:29.976_1.nod.ms.norm.sum.wavecal.fits"):


    header = fits.getheader(fname)

    
    print(header["DATE-OBS"])
     # Observation Specifications
    # Date
    #MJD-OBS =       56024.00590250 / Obs start 2012-04-07T00:08:29.976
    #DATE-OBS= '2012-04-07T00:08:29.9764' / Observing date
    #EXPTIME =          180.0000000 / Integration time

    date = header["DATE-OBS"][:-1]+"Z"     

    print("TODO !!! - Time correct for tapas timing issue")
    print("date", date)
    #date = "2013-11-17T00:41:00.000Z"
    if "VLT" in header["TELESCOP"]:
    	obs_name = "ESO Paranal Chile (VLT)"
    else:
    	print("Don't have TAPAS telecope location name (just the fits one)")
    	obs_name = header["TELESCOP"]

    obs_long = header["HIERARCH ESO TEL GEOLON"]
    obs_lat = header["HIERARCH ESO TEL GEOLAT"]
    obs_alt = header["HIERARCH ESO TEL GEOELEV"]

    # Target 
    # Procees cooords
    from astropy import units as u
    from astropy.coordinates import Angle
    ra_angle = Angle(header["RA"], u.degree)
    ra_j2000 = str(ra_angle.to_string(unit=u.hour, sep=':',precision=0, pad=True))  # Extra str to get rid of u"string" which failed in template
    print("Angle =",ra_angle, "header ra =",header["RA"], "ra_2000 ::=",ra_j2000)
    dec_angle = Angle(header["DEC"], u.deg)
    dec_j2000 = str(dec_angle.to_string(unit=u.degree, sep=':', precision=0, pad=True)) # Extra str to get rid of u"string" which failed in template
    print("Angle =",dec_angle, "header DEC =", header["DEC"], "dec_2000 ::=", dec_j2000)

    print(dec_j2000)



    #spectral_range = "{0} {1}".format(min_range, max_range)
    spec_range_min = round(header["HIERARCH ESO INS WLEN MIN"] - 10)
    spec_range_max = round(header["HIERARCH ESO INS WLEN MAX"] + 10) 


    if "CRIRES" in header["INSTRUME"]:
        if header["HIERARCH ESO INS SLIT1 WID"] == 0.400:
            print("Slit width was 0.4, \nSetting resolving_power = 50000")
            resolving_power = 50000    # if Crires take the two slit width/Resolution values
        elif header["HIERARCH ESO INS SLIT1 WID"] == 0.200:
            print("Slit width was 0.2, \nSetting resolving_power = 100000")
            resolving_power = 100000    # if Crires take the two slit width/Resolution values
        else:
        	print("Slitwidth of CRIRES does not match 2 fixed settings")
        	print("Is this older data?")

# Things want out of observations files 
#["DATE-OBS"]
#TELESCOP= 'ESO-VLT-U1'
#RA      =            71.409028          / 04:45:38.1 RA (J2000) pointing (deg)
#DEC     =            -50.07734 
#HIERARCH ESO INS SLIT1 WID   =        0.400 
#INSTRUME= 'CRIRES  '
#HIERARCH ESO TEL GEOELEV     =         2648. / Elevation above sea level (m)
#HIERARCH ESO TEL GEOLAT      =     -24.6276 / Tel geo latitute (+=North) (deg)
#HIERARCH ESO TEL GEOLON      =     -70.4051 / Tel geo longitude (+=East) (deg)
#HIERARCH ESO INS WLEN MAX    =     2160.205 / Wavelength Limit + [nm].
#HIERARCH ESO INS WLEN MIN    =     2121.056 / Wavelength Limit - [nm].

#Test_file = 
    #Template_xml_file = "/home/jneal/Phd/Codes/UsefulModules/tapas_xml_template.txt"
    
    # Tapas Specifications 
    Format = "ASCII"
    Resquest_ID = "10000"    # Iterate on previous ID if found
    ray = "YES"
    h20 = "NO"   # alternate yes no for h20 for scaling
    o3 = "YES"
    o2 = "YES"
    co2 = "YES"
    ch4 = "YES"
    n2o = "YES"
    #"if Air:
    spectral_choice = "Standard Wavelength (nm)" 
    #"if vac"
    #"if wavenumber"
    ilsf_choice = 1    # -1, 0, 1
   
    sampling_ratio = 10
    apply_berv = "YES"


    #arletty_file = "canr_2013111700.arl"
    #ecmwf_file = "canr_2013111700_qo3.txt"
    #date = "2013-11-17T00:41:00.000Z"
    arletty_file = "canr_"+ date[0:4] + date[5:7] + date[8:10] + date[11:13] + ".arl"
    ecmwf_file = "canr_" + date[0:4] + date[5:7] + date[8:10] + date[11:13] + "_qo3.txt"
    
    print("arletty_file", arletty_file)
    print("ecmwf_file", ecmwf_file)

    d = {"Resquest_ID":Resquest_ID, "Format":Format, "ray":ray, "h20":h20, "o3":o3, "o2":o2, "co2":co2, "ch4":ch4, "n2o":n2o, "date":date, "obs_name":obs_name, "obs_long":obs_long, "obs_lat":obs_lat, "obs_alt":obs_alt, "ra_j2000":ra_j2000, "dec_j2000":dec_j2000, "spectral_choice":spectral_choice, "spec_range_min":spec_range_min, "spec_range_max":spec_range_max, "ilsf_choice":ilsf_choice, "resolving_power":resolving_power, "sampling_ratio":sampling_ratio, "apply_berv":apply_berv, "arletty_file":arletty_file, "ecmwf_file":ecmwf_file}

    template = """‹?xml version="1.0" encoding="UTF-8"?›
<tapas Id="Ether_TAPAS_999"›
<request Id="$Resquest_ID"›
<preferences>
<format valid="VO,ASCII,FITS,NETCDF">$Format</format>
<rayleigh_extinction valid="YES,NO">$ray</rayleigh_extinction>
<h2o_extinction valid="YES,NO">$h20</h2o_extinction>
<o3_extinction valid="YES,NO">$o3</o3_extinction>
<o2_extinction valid="YES,NO">$o3</o2_extinction>
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
#    args = vars(_parser())
#    fname = args.pop('fname')
#    opts = {k: args[k] for k in args}

#    main(fname, **opts)
    main()
