#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH -o brun_%A.out
#SBATCH --job-name=gC

module purge
module load python
module load gcc
module load boost
module load matlab

date

echo $theme

./gainCurve -m $theme && matlab -nodisplay -nosplash -r "plotGainCurve('input.cfg');exit;"

date
hostname
