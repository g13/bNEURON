#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=9
#SBATCH --mem-per-cpu=5G
#SBATCH --time=8:00:00
#SBATCH -a 1-4
#SBATCH -o mpiLib_%A_%a.out
#SBATCH --job-name=mpiLib
#SBATCH --mail-type=FAIL,END

module purge
module load python/2.7.11
module list

date
echo $theme
echo $SLURM_ARRAY_TASK_ID

python getNib_mpi.py --theme=${theme} --vid=$(($SLURM_ARRAY_TASK_ID-1))

date
hostname
