#!/bin/bash

#SBATCH -o Disp1.o%j
#SBATCH --time=48:00:00
#SBATCH -N 1 -n 4 --mem=80gb

cp SlowAnalysis16.m $PFSDIR
cp PointsTestMat16.mat $PFSDIR
cp PointsTestNum16.mat $PFSDIR

module load matlab

# MATLAB Preference Setting
matlab_prefdir="/tmp/$USER/matlab/`hostname`_PID$"
test -d $matlab_prefdir || mkdir -p $matlab_prefdir
export MATLAB_PREFDIR="$matlab_prefdir"

#RUN MATLAB script
matlab -singleCompThread -nodisplay -r SlowAnalysis16

quit

cp -r PreDisp16.mat $SLURM_SUBMIT_DIR