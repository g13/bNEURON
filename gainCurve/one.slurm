#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --mem=60G
#SBATCH --cpus-per-task=15
#SBATCH --time=4:00:00
#SBATCH -o gC_%A.out
#SBATCH --job-name=gC

module purge
module load python
module load gcc
module load boost
module load matlab
module load anaconda

date
echo $ext

source activate py2

for method in {0..5}
do
    ./one -d $method
done

matlab -nodisplay -nosplash -r "plotGainCurve('input.cfg','$ext',true,false);exit;"

date
hostname
