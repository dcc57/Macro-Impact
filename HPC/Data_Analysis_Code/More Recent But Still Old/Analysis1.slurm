#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=80g
#SBATCH -J Analysis1

module load python

cp -r /home/dcc57/HPCAnalysis/HPCDataAnalysis.py $PFSDIR
cp -r /home/dcc57/HPC/PointsTestNum1.mat $PFSDIR
cp -r /home/dcc57/HPC/PointsTestMat1.mat $PFSDIR

export PYTHONPATH=$PYTHONPATH:/home/dcc57/pip/env/lib/python2.7/site-packages

cd $PFSDIR
echo "Working in: " . $PFSDIR

python HPCDataAnalysis.py 4 1.737E6 1.0E3 1.6E1 1

cp -r PreDisplacement1.mat $SLURM_SUBMIT_DIR