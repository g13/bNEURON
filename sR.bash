#!/bin/bash
n128=0
theme='massive-spike'
dir0='singleRun'
dir1=$dir0'/'$theme
if [ -d "$dir1" ]; then
    rm -r $dir1/*
else
    mkdir $dir1
fi

if [ "$n128" == "0" ]; then
    cp base/n128.py_pas $dir1/n128.py
else
    cp base/n128.py_active $dir1/n128.py
fi
cp base/neuroAlter.py $dir1
cp $dir0/sr.slurm $dir1
cp -R x86_64 $dir1

cd $dir1

sbatch --export=theme=$theme sr.slurm
