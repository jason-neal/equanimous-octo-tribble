#!/bin/bash
# Loop over obs folders and generate the key_obs_values.txt files
for FOLDER in $(ls -d HD*); do
    echo $FOLDER
    cp $FOLDER/key_observation_values.txt  /run/media/jneal/SP\ UFD\ U2/Data/$FOLDER/
done