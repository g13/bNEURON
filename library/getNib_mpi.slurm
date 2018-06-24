#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=6G
#SBATCH --time=16:00:00
#SBATCH -a 1-10
#SBATCH -o mpiLib_%A_%a.out
#SBATCH --job-name=mpiLib
#SBATCH --mail-type=FAIL,END

module purge
module load anaconda
module list

date
set -e
source activate py2
echo $theme
echo $SLURM_ARRAY_TASK_ID

python getNib_mpi.py --theme=${theme} --vid=$(($SLURM_ARRAY_TASK_ID-1))

date
hostname
