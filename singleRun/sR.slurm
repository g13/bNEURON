#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --time=0:10:00
#SBATCH -o brun_%A.out
#SBATCH --job-name=sR

module purge
module load python
module load gcc
module load boost
module load matlab

date

echo $lib
#echo ${lib:3}
echo $v0
echo $theme
echo $rE
echo $rI
echo $t
echo $tref
echo $seed
./sR -m $theme -v $v0 -r $rE -r $rI -t $t --tref=$tref -s $seed -L $lib && matlab -r "bnsynCompare('$theme',$rE,$rI,'$lib',$seed,$t);exit;"

date
hostname
