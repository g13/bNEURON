#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --time=0:15:00
#SBATCH -o brun_%A.out
#SBATCH --job-name=brun

module purge
module load python
module load gcc
module load boost
module load matlab

date

./bsN -t $theme > blog 2>&1 && matlab -r "bnsynCompare;exit;"

date
hostname
