#!/bin/bash
echo "Generating list of names and object types and expouse times "
for NAME in $(ls C*.fits); do
    #echo $NAME
    echo "${NAME%.*}"  
    dfits $NAME | grep OBJECT 
    dfits $NAME | grep EXPTIME
    
    echo 
    #fitsort $NAME OBJECT,EXPTIME

    # could maybe auto generate lists but going to jsut copy and paste from the generated list today.
done