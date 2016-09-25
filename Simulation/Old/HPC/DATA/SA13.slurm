#!/bin/bash

#SBATCH -o Disp1.o%j
#SBATCH --time=24:00:00
#SBATCH -N 1 -n 1 --mem=40gb

cp SlowAnalysis13.m $PFSDIR
cp PointsTestMat13.mat $PFSDIR
cp PointsTestNum13.mat $PFSDIR

module load matlab

# MATLAB Preference Setting
matlab_prefdir="/tmp/$USER/matlab/`hostname`_PID$"
test -d $matlab_prefdir || mkdir -p $matlab_prefdir
export MATLAB_PREFDIR="$matlab_prefdir"

#RUN MATLAB script
matlab -singleCompThread -nodisplay -r SlowAnalysis13

quit

cp -r PreDisp13.mat $SLURM_SUBMIT_DIR