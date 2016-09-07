#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --mem=10g
#SBATCH -J Analysis1

module load python

cp -r /home/dcc57/HPCAnalysis/HPCDataAnalysis.py $PFSDIR
cp -r /home/dcc57/HPC/PointsTestNum1Test.mat $PFSDIR
cp -r /home/dcc57/HPC/PointsTestMat1Test.mat $PFSDIR

export PYTHONPATH=$PYTHONPATH:/home/dcc57/pip/env/lib/python2.7/site-packages

cd $PFSDIR
echo "Working in: " . $PFSDIR

python HPCDataAnalysis.py 1 1.737E6 1.0E3 1.6E1 1Test

cp -r PreDisplacement1Test.mat $SLURM_SUBMIT_DIR