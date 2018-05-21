#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --cpus-per-task=1
#SBATCH --time=0:30:00
#SBATCH -o readLib_%A.out

module purge
module load matlab

date

theme='active_fast'
echo $theme
matlab -nodesktop -nosplash -r "readNib4('$theme');exit;"

date
hostname
