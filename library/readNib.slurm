#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --time=0:15:00
#SBATCH -o readLib_%A.out
#SBATCH --job-name=readLib

module purge
module load matlab

date

theme='act-BT'
echo $theme
matlab -nodesktop -nosplash -r "readNib4('$theme');exit;"

date
hostname
