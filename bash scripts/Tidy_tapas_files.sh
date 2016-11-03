# Extract and rename tapas files


declare -a arr=("HD202206-1" "HD202206-2" "HD202206-3"
 "HD162020-1" "HD162020-2"
 "HD167665-1a" "HD167665-1b" "HD167665-2"
 "HD168443-1" "HD168443-2"
 "HD211847-1" "HD211847-2"
 "HD4747-1")

declare -a arr2=("HD30501-1" "HD30501-2a" "HD30501-2b" "HD30501-3")

for FOLDER in "${arr[@]}"
do
    echo "Tidying telluric files in", $FOLDER
    cd ~/Phd/data/Crires/BDs-DRACS/$FOLDER/

    mkdir Telluric_files
    mv tapas_*ipac Telluric_files/
    mv *tapas*request*.xml Telluric_files/

done

for FOLDER in "${arr2[@]}"
do
    echo "Tidying telluric files in", $FOLDER
    cd ~/Phd/data/Crires/BDs-DRACS/$FOLDER/

    mkdir Telluric_files
    mv Combined_Nods/tapas_*ipac Telluric_files/
    mv ~/Phd/Codes/UsefulModules/requests/$FOLDER*xml-request*.txt Telluric_files/

done