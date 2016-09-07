#!/bin/bash

declare -a arr=("HD202206-1" "HD202206-2" "HD202206-3"
 "HD162020-1" "HD162020-2"
 "HD167665-1a" "HD167665-1b" "HD167665-2"
 "HD168443-1" "HD168443-2"
 "HD211847-1" "HD211847-2"
 "HD4747-1")

## now loop through the above array
for FOLDER in "${arr[@]}"
do
    echo "Request for ", $FOLDER
    #python tapas_xml_request_generator.py ../../data/Crires/BDs-DRACS/$FOLDER/list_spectra.txt -f none -u vacuum -c h2o -l -o "../../data/Crires/BDs-DRACS/"$FOLDER"/"$FOLDER"_tapas_h2o_bash_request.xml"
    python tapas_xml_request_generator.py ../../data/Crires/BDs-DRACS/$FOLDER/list_spectra.txt -f none -u vacuum -c h2o -l -o "bash_requests/"$FOLDER"_tapas_h2o_bash_request.xml"
 
done
