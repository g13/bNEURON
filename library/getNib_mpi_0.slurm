#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00
#SBATCH -a 1-7
#SBATCH -o mpiLib0_%A.out
#SBATCH --job-name=mpiLib0

module purge
module load python/2.7.11
module list

date

echo $theme
python getNib_mpi.py --theme=${theme} --vid=0 -r

date
hostname
