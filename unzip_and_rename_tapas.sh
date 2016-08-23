# Extract and rename tapas files


declare -a arr=("HD202206-1" "HD202206-2" "HD202206-3"
 "HD162020-1" "HD162020-2"
 "HD167665-1a" "HD167665-1b" "HD167665-2"
 "HD168443-1" "HD168443-2"
 "HD211847-1" "HD211847-2"
 "HD4747-1")

#declare -a arr=("HD202206-2")

 for FOLDER in "${arr[@]}"
do
    echo "Moving ", $FOLDER, "xml files"
    cd ~/Phd/data/Crires/BDs-DRACS/$FOLDER/
    rm $FOLDER_*h2o_bash_request.xml

    python ~/Phd/Codes/UsefulModules/tapas_spectra_namer.py tapas_000010.ipac.gz -k -x -s all

    python ~/Phd/Codes/UsefulModules/tapas_spectra_namer.py tapas_000012.ipac.gz -k -x -s h2o

    python ~/Phd/Codes/UsefulModules/tapas_spectra_namer.py tapas_000012_3*.ipac.gz -k -x -s h2o

    python ~/Phd/Codes/UsefulModules/tapas_spectra_namer.py tapas_000018.ipac.gz -k -x -s noh2o

done