#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=5G
#SBATCH --time=08:00:00
#SBATCH -o mpiVib_%A.out
#SBATCH --job-name=mpiVib
#SBATCH --mail-type=FAIL,END

module purge
module load anaconda
module list

date
set -e
source activate py2

#python sv.py 
python testDvOnly.py
#--theme=${theme} --vid=$(($SLURM_ARRAY_TASK_ID-1)) --nc=$nc

date
hostname
