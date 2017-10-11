#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH -o fLsp_%A.out
#SBATCH --job-name=fLsp
#SBATCH --mail-type=FAIL,END

module purge
module load python/2.7.11
module list
date

if [ "$mode" == "0" ]; then
    python fL0_spList.py
else
    python fL_spList.py
fi

date
