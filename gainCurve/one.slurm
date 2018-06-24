#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --time=0:20:00
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

source activate py2

for method in  {0..6}
do
    ./one -d $method
done

echo $ext
echo $getDendV

matlab -nodisplay -nosplash -r "plotGainCurve('input.cfg','$ext',true,true,$getDendV);exit;"

date
hostname
