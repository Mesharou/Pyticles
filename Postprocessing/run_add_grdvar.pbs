#!/bin/bash
#PBS -q omp
#PBS -l ncpus=14
#PBS -l mem=65g
#PBS -l walltime=24:00:00
#  COMMENT
#  Tina ODAKA 06.03.2017
#  example for using datarmor with 'queue omp' in a personalised way
#  By default, if you use 'omp' it use 1 node with 2 openmp threads
#  here we show you example of it using 10 openmp threads (specify ncpus=10)
#  and it use 10g of memory in this example (specifly mem=10g)
#  and you can run max for 30 minutes with walltime setting
#  To optimised multi thread, you need to place a command 'omplace' before your 'exe'
#  Here, in this example, we suppose that your code 'exe' is compiled with 
#  intel compiler.
# for the users who use gnu compiler, plz use omp.10g.gnu.pbs
#  
# cd to the directory you submitted your job
cd $PBS_O_WORKDIR
qstat -f $PBS_JOBID

source /usr/share/Modules/3.2.10/init/csh
module list

module load anaconda-py2.7/4.3.13
cd /appli/anaconda/2.7/bin
source activate /home2/datawork/jgula/POLGYR/my_conda

echo $LD_LIBRARY_PATH


cd /home2/datawork/jgula/Pyticles/Pyticles_nesea/Postprocessing

export PYTHONPATH=$HOME/Python_Modules


echo "you use " $OMP_NUM_THREADS "threads for your omp jobs"
date


for file in `ls ../north/*.nc`; do
    echo $file
    python add_grdvar_to_file.py $file > output
done
date
