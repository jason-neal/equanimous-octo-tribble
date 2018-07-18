#!/bin/bash
# For post reduction steps after IRAF processing
# To organize dracs data directories and automate
# some post reduction stages.
# Author: Jason Neal

mkdir Combined_Nods
mv CRIRE*_*.nod.ms.norm.*.fits  Combined_Nods/

mkdir Intermediate_steps
mv CRIRE*_*.fits  Intermediate_steps/
mv Flat*.fits Intermediate_steps/
mv MasterDark*.fits Intermediate_steps/
mv list*.txt_* Intermediate_steps/
cp list_spectra.txt Combined_Nods/

mkdir Raw_files
mv CRIRE*.fits Raw_files/
mv M.CRIRES*.fits Raw_files/

# Make Quicklook plots
mkdir images
dracs_quicklooks2017.py

# Replace : characters in reduced data in Combined_Nods data
file_renamer.py ./Combined_Nods/* -x ":" -v "-"
