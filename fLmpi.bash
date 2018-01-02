#!/bin/bash

mode=1 # 0 with vClamp on soma only
n128=1
dir0='freqLimit'
seed=26574839
theme0='test'
theme=$theme0'-n'$n128-'m'$mode
## p-preset; r-random
dir=$dir0'/'$theme

if [ -d "./$dir" ]; then
    rm -r $dir/*
else
    mkdir $dir
fi
cp $dir0/freqMpi.slurm $dir
cp $dir0/freqLimit_mpi.py $dir

if [ "$n128" == "0" ]; then
    cp base/n128.py_pas $dir/n128.py
else
    cp base/n128.py_active $dir/n128.py
fi
cp base/neuroAlter.py $dir
cp -R x86_64 $dir

cd $dir

export mode
export seed
sbatch --export=ALL freqMpi.slurm
