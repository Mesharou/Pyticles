#!/usr/bin/env bash



local_folder=`pwd`
source_code='/home2/datahome/jgula/code/Pyticles_master'

###############


wsed=20

dfile=61 #each experiment is 30.5 days

start_file0=17280
end_file0=$(($start_file0+12*$dfile))

echo $start_file0,$end_file0

###############


for (( start_file=$start_file0; start_file<$end_file0; start_file+=$dfile )); do

    end_file=$(($start_file+$dfile+30))

    echo $start_file,$end_file

    # Define a folder for the exp.
    exp_number=`printf %06d $start_file`
    folder_out=experiment_$exp_number
    mkdir $folder_out

    ln -s $source_code/Modules $folder_out/
    ln -s $source_code/Pyticles_subroutines $folder_out/
    mkdir $folder_out/Inputs

    sed -e 's/WSED/'$(($wsed))'/g;
            s/STARTFILE/'$(($start_file))'/g;
            s/ENDFILE/'$(($end_file))'/g'  $source_code/Inputs/input_file.py.generic > $folder_out/Inputs/input_file.py

    cp -rf run_Pyticles.pbs $folder_out/
    cp -rf Pyticles.py $folder_out/

    cd $folder_out
    #qsub run_Pyticles.pbs
    cd $local_folder

done


