#!/bin/bash
# Loop over obs folders and generate the key_obs_values.txt files
for FOLDER in $(ls -d HD*); do
    echo $FOLDER
    cd $FOLDER
    python ~/Phd/Codes/UsefulModules/obs_info.py list_spectra.txt 
    cd ..
done