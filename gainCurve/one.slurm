#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --mem=60G
#SBATCH --cpus-per-task=15
#SBATCH --time=8:00:00
#SBATCH -o one_%A.out
#SBATCH --job-name=one
#SBATCH --mail-type=END,FAIL

module purge
module load python
module load gcc
module load boost
module load matlab
module load anaconda

set -e
date
echo $ext
echo $getDendV

source activate py2

for method in {0..5}
do
    ./one -d $method
done

matlab -nodisplay -nosplash -r "plotGainCurve('input.cfg','$ext',true,false,$getDendV);exit;"

date
hostname
