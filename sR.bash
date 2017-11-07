#!/bin/bash
n128=0
theme='massive-spike'
lib='mpiTest-massive-NEURON.mat'
rE=40
rI=60
t=1000
v0=1
tref=13
seed=193864

datafdr='../../data'
dir0='singleRun'
dir=$dir0'/'$theme
lib=$datafdr'/'$lib
if [ -d "$dir" ]; then
    rm -r $dir/*
else
    mkdir $dir
fi

if [ "$n128" == "0" ]; then
    cp base/n128.py_pas $dir/n128.py
else
    cp base/n128.py_active $dir/n128.py
fi
cp base/neuroAlter.py $dir
cp $dir0/sR.slurm $dir
cp $dir0/sR_init.cfg $dir
cp $dir0/sR $dir
cp $dir0/bnsynCompare.m $dir
cp -R x86_64 $dir

cd $dir

export lib
export v0
export theme
export rE
export rI
export t
export tref
export seed
sbatch --export=ALL sR.slurm
