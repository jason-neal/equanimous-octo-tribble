#!/bin/bash
# Loop over observation folders and generate run the dracs quicklooks code files
for FOLDER in $(ls -d HD*); do
    echo "Remaking quicklooks in $FOLDER"
    cd $FOLDER
    python ../dracs_quicklooks2017.py
    cd ..
done
