#!/bin/bash

ROMS1=$1
#ROMS1=$(pwd)/$1
#ROMS1=/home/nmolem/ROMS_atlbig
#ROMS2=$2
ROMS2=$(pwd)/$2
#ROMS2=/home/nmolem/ROMS_atlbig


cd $ROMS1


for i in `find . -maxdepth 5 -type f \( -name "*.py" \)`; do

        diff_output=`diff --ignore-all-space $i $ROMS2/$i`

        if  [ -n "$diff_output" ]; then
            echo $i
            diff --ignore-all-space $i $ROMS2/$i
            echo '   '
            echo '   '
            echo '   '
        fi

done


