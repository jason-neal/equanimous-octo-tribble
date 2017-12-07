#!/bin/bash
# Script to deorganize dracs data directory after it has been organized with my other script. 
# This is to run extra iraf scripts if needed.
# Author: Jason Neal

#mkdir Combined_Nods
mv Combined_Nods/CRIRE*_*.nod.ms.norm.*.fits .

#mkdir Intermediate_steps
mv Intermediate_steps/CRIRE*_*.fits .
mv Intermediate_steps/Flat*.fits .
mv Intermediate_steps/MasterDark*.fits .
mv Intermediate_steps/list*.txt_* .
#cp list_spectra.txt Combined_Nods/

#mkdir Raw_files
mv Raw_files/CRIRE*.fits .
mv Raw_files/M.CRIRES*.fits .
