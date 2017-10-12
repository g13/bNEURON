#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=08:00:00
#SBATCH --mem=8G
#SBATCH -o fLm_%A.out
#SBATCH --job-name=mpiTest
#SBATCH --mail-type=FAIL,END

module purge
module load python/2.7.11 openmpi/1.10.2
module list
date

if [ "$mode" == "0" ]; then
    python freqLimit0_mpi.py
else
    python freqLimit_mpi.py
fi

date
