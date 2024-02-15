# equanimous-octo-tribble

[![Python package](https://github.com/jason-neal/equanimous-octo-tribble/actions/workflows/python-package.yml/badge.svg?branch=master)](https://github.com/jason-neal/equanimous-octo-tribble/actions/workflows/python-package.yml)[![Codacy Badge](https://api.codacy.com/project/badge/Grade/7fc5c278070a4ff188f5c88cccd6be64)](https://www.codacy.com/app/jason-neal/equanimous-octo-tribble?utm_source=github.com&utm_medium=referral&utm_content=jason-neal/equanimous-octo-tribble&utm_campaign=badger)

Useful scripts for astronomy and the like.

## Convolution
Convole a spectrum by a IP of a given Resolution

There is a mulitprocessing version and a normal (slower) version.
Calculates the IP for every pixel/wavelength value individually.
Does not need a eqidistant wavelength axis.
Assumes a gaussian IP profile


## tapas_xml_request_generator.py
Creates an xml request to submit to get tapas data. 
You must copy the text created in the .xml file into the xml box on the tapas request page (The button is three left of "Exucute"). 

The script takes a single fits spectra or a nod cycle list of spectra and determines the parameters needed for tapas. e.g. Date/Time, RA,DEC, location, resolution.


The tapas request number needs to be iterated. If your request fails try the request again but with the request number changed to one more than the number in the email you received.

Most parameters can be changed by adding extra command line arguments.
E.g.   -r --resolution, -s --sampling

To see the full list of options type 
    python tapas_xml_request_generator.py -h 

Example of useage:

    python tapas_xml_request_generator.py CRIRES.fits 

    python tapas_xml_request_generator.py list_spectra.txt -l --sampling 5 --constituents h20
    python tapas_xml_request_generator.py list_spectra.txt -l -s 15 -c all-n 3050 --verbose --unit air --berv --tapas-format fits


## tapas_spectra_namer.py

Uniquely name the spectra downloaded from tapas out of the header information.
E.g. Time, Resolution, wavelenght scale

Can also add your own prefix to personally identify file

To see the full list of options type
python tapas_spectra_namer.py -h

Example of useage:

    python tapas_spectra_name.py tapas_00001.ipac --prefix HD30501-1

    python tapas_spectra_name.py tapas_00001.ipac.gz --extract --keep-original
