#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=100g
#SBATCH -J Analysis12
 
module load python
 
cp -r /home/dcc57/HPCAnalysis/HPCDataAnalysis.py $PFSDIR
cp -r /home/dcc57/HPC/PointsTestNum12.mat $PFSDIR
cp -r /home/dcc57/HPC/PointsTestMat12.mat $PFSDIR

export PYTHONPATH=$PYTHONPATH:/home/dcc57/pip/env/lib/python2.7/site-packages

cd $PFSDIR
echo "Working in: " . $PFSDIR
 
python HPCDataAnalysis.py 4 1.737E6 1.0E3 1.6E1 12
 
cp -r PreDisplacement12.mat $SLURM_SUBMIT_DIR