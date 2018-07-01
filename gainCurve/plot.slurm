#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1
#SBATCH --time=0:05:00
#SBATCH -o plot_%A.out
#SBATCH --job-name=plot
#SBATCH --mail-type=END,FAIL

module purge 
module load matlab

set -e
date 
echo $ext
echo $getDendV

plotSub=true
plotAuto=true
plotInput=true
cfgFile=input.cfg
sizeSize=int64
matlab -nodisplay -nosplash -r "plotMultiGainCurve('$cfgFile','$ext',$plotSub,$plotInput,$getDendV,'$sizeSize',$plotAuto);exit;"

date
hostname
