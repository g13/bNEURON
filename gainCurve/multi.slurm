#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=7
#SBATCH --time=24:00:00
#SBATCH -a 1-7
#SBATCH -o m_%A_%a.out
#SBATCH --job-name=multi
#SBATCH --mail-type=END,FAIL

module purge
module load python
module load gcc
module load boost
module load matlab
module load anaconda

set -e
date

source activate py2
echo $nMethod
method=$(($SLURM_ARRAY_TASK_ID-1))

if [ $SLURM_ARRAY_TASK_ID == $nMethod ]; then
    ./multi -d $method --tIn=1
    echo "output tIncome with method $method"
else
    ./multi -d $method
fi

date
hostname
