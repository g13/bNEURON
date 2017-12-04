#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --time=0:10:00
#SBATCH -o brun_%A.out
#SBATCH --job-name=gC

module purge
module load python
module load gcc
module load boost
module load matlab

date

echo $theme

dataFn=./gainCurve -m $theme 
matlab -r "plotGainCurve('$dataFn');exit;"

date
hostname
