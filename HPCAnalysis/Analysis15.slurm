#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=20g
#SBATCH -J Analysis15
 
module load python
 
cp -r /home/dcc57/HPCAnalysis/HPCDataAnalysis.py $PFSDIR
cp -r /home/dcc57/HPC/PointsTestNum15.mat $PFSDIR
cp -r /home/dcc57/HPC/PointsTestMat15.mat $PFSDIR
cd $PFSDIR
echo "Working in: " . $PFSDIR
 
python HPCDataAnalysis.py 4 1.737E6 1.0E3 1.6E1 15
 
cp -r PreDisplacement15.mat $SLURM_SUBMIT_DIR