#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=5G
#SBATCH --time=8:00:00
#SBATCH -o mpiLibi_%A.out
#SBATCH --job-name=mpiLibi
#SBATCH --mail-type=FAIL,END

module purge
module load python/2.7.11
module list

date
theme=mpiTest-active
id=4
echo $theme

python getNib_mpi.py --theme=${theme} --vid=$(($id-1))

date
hostname
