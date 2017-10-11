#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=40G
#SBATCH --cpus-per-task=1
#SBATCH --time=0:15:00
#SBATCH -o readLib_%A.out

module purge
module load matlab

date

#theme=stem-L13dt-dend-I
theme=mpiTest-massive
echo $theme
matlab -nodesktop -nosplash -r "readNib4('$theme');exit;"

date
hostname
