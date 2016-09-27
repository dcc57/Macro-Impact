#!/bin/bash

#SBATCH -o Generation1.o%j
#SBATCH --time=24:00:00
#SBATCH -N 1 -n 1 --mem=50gb


#Copy to the scratch directory and change directory
cp -r Ray-Tracing $PFSDIR
cd $PFSDIR/HPC/DataGeneration

#Load MATLAB module
module load matlab


# MATLAB Preference Setting
matlab_prefdir="/tmp/$USER/matlab/`hostname`_PID$"
test -d $matlab_prefdir || mkdir -p $matlab_prefdir
export MATLAB_PREFDIR="$matlab_prefdir"

#RUN MATLAB script
matlab -singleCompThread -nodisplay -r Generation1

quit

echo Master process running on $HOSTNAME
echo Directory is `pwd`
NPROCS=$(( $SLURM_NPROCS*$SLURM_NNODES))
echo This job is allocated $NPROCS nodes
echo `date`

#Copy Back all the output files to the present working directory
cp -r Mat1.mat $SLURM_SUBMIT_DIR
cp -r Num1.mat $SLURM_SUBMIT_DIR