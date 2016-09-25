#!/bin/bash

#SBATCH -o Disp1.o%j
#SBATCH --time=24:00:00
#SBATCH -N 1 -n 1 --mem=40gb

cp SlowAnalysis1.m $PFSDIR
cp PointsTestMat1.mat $PFSDIR
cp PointsTestNum1.mat $PFSDIR

module load matlab

# MATLAB Preference Setting
matlab_prefdir="/tmp/$USER/matlab/`hostname`_PID$"
test -d $matlab_prefdir || mkdir -p $matlab_prefdir
export MATLAB_PREFDIR="$matlab_prefdir"

#RUN MATLAB script
matlab -singleCompThread -nodisplay -r SlowAnalysis1

quit

cp -r PreDisp1.mat $SLURM_SUBMIT_DIR