#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=20g
#SBATCH -J Analysis10
 
module load python
 
cp -r /home/dcc57/HPCAnalysis/HPCDataAnalysis.py $PFSDIR
cp -r /home/dcc57/HPC/PointsTestNum10.mat $PFSDIR
cp -r /home/dcc57/HPC/PointsTestMat10.mat $PFSDIR
cd $PFSDIR
echo "Working in: " . $PFSDIR
 
python HPCDataAnalysis.py 4 1.737E6 1.0E3 1.6E1 10
 
cp -r PreDisplacement10.mat $SLURM_SUBMIT_DIR