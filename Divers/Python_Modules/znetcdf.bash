#!/bin/bash
# use it as : 
# ./znetcdf.bash files*

#Log file
COMPRESSED=$(dirname $1)/COMPRESSED

# find /shirma/* -type f \( -name "*.nc" \) > selection
date >> $COMPRESSED
echo 'The following filess have been converted to netcdf4 with zlib - using nc3tonc4 (complevel=6)' >> $COMPRESSED

netcdf3='64-bit offset'
netcdf3_classic='classic'
netcdf4='netCDF-4'
netcdf4_classic='netCDF-4 classic model'

for i in $@; do
# for i in `cat selection`; do
    echo $i
    format=`ncdump -k $i`
    if [ "$format" = "$netcdf3" -o "$format" = "$netcdf3_classic" ]; then
        echo 'format is', $format
        nc3tonc4 --chunk=1 --complevel=6 $i $i.znc > $COMPRESSED.log
        ./check_netcdf_files.py $i $i.znc $COMPRESSED.log
        echo $i >> $COMPRESSED
    fi

done


echo 'All files have been converted to netcdf4 with zlib - using nc3tonc4 (complevel=6)' >> $COMPRESSED
date >> $COMPRESSED




