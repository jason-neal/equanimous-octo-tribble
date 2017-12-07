#!/bin/bash
echo "Generating list of Crires files without extensions."
for NAME in $(ls C*.fits); do
    #echo $NAME
    echo "${NAME%.*}"  

done