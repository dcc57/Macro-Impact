#!/bin/bash

#SBATCH -o Disp1.o%j
#SBATCH --time=48:00:00
#SBATCH -N 1 -n 4 --mem=80gb

cp SlowAnalysis17.m $PFSDIR
cp PointsTestMat17.mat $PFSDIR
cp PointsTestNum17.mat $PFSDIR

module load matlab

# MATLAB Preference Setting
matlab_prefdir="/tmp/$USER/matlab/`hostname`_PID$"
test -d $matlab_prefdir || mkdir -p $matlab_prefdir
export MATLAB_PREFDIR="$matlab_prefdir"

#RUN MATLAB script
matlab -singleCompThread -nodisplay -r SlowAnalysis17

quit

cp -r PreDisp17.mat $SLURM_SUBMIT_DIR