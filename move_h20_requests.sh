#!/bin/bash
# Move tapas request created and requested to their respected star folders

declare -a arr=("HD202206-1" "HD202206-2" "HD202206-3"
 "HD162020-1" "HD162020-2"
 "HD167665-1a" "HD167665-1b" "HD167665-2"
 "HD168443-1" "HD168443-2"
 "HD211847-1" "HD211847-2"
 "HD4747-1")

## now loop through the above array
for FOLDER in "${arr[@]}"
do
    echo "Moving ", $FOLDER, "xml files"
    #python tapas_xml_request_generator.py ../../data/Crires/BDs-DRACS/$FOLDER/list_spectra.txt -f none -u vacuum -c h2o -l -o "../../data/Crires/BDs-DRACS/"$FOLDER"/"$FOLDER"_tapas_h2o_bash_request.xml"
    mv bash_requests/$FOLDER*.xml ../../data/Crires/BDs-DRACS/$FOLDER/
done