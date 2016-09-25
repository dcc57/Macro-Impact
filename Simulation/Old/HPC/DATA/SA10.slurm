#!/bin/bash

#SBATCH -o Disp1.o%j
#SBATCH --time=24:00:00
#SBATCH -N 1 -n 1 --mem=40gb

cp SlowAnalysis10.m $PFSDIR
cp PointsTestMat10.mat $PFSDIR
cp PointsTestNum10.mat $PFSDIR
cd $PFSDIR/Ray-Tracing

module load matlab

# MATLAB Preference Setting
matlab_prefdir="/tmp/$USER/matlab/`hostname`_PID$"
test -d $matlab_prefdir || mkdir -p $matlab_prefdir
export MATLAB_PREFDIR="$matlab_prefdir"

#RUN MATLAB script
matlab -singleCompThread -nodisplay -r SlowAnalysis10

quit

cp -r PreDisp10.mat $SLURM_SUBMIT_DIR