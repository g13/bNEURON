#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4G
#SBATCH --time=12:00:00
#SBATCH -o mpiLib0_%A.out
#SBATCH --job-name=mpiLib0

module purge
module load anaconda
module list

date
source activate py2

echo $theme
echo $nc
python getNib_mpi.py --theme=${theme} --vid=0 -r --nc=$nc
date
hostname
