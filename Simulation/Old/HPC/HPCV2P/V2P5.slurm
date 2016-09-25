#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=80g
#SBATCH -J Analysis3
 
module load python
 
cp -r /home/dcc57/HPCAnalysis/HPCVec2Pix.py $PFSDIR
cp -r /home/dcc57/HPC/PointsTestMat5.mat $PFSDIR

export PYTHONPATH=$PYTHONPATH:/home/dcc57/pip/env/lib/python2.7/site-packages

cd $PFSDIR
echo "Working in: " . $PFSDIR
 
python HPCVec2Pix.py 4 5
 
cp -r hppoints5.mat $SLURM_SUBMIT_DIR