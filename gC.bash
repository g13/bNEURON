#!/bin/bash
n128=0
datafdr='../../data'
dir0='gainCurve'

theme=$1

if [ "$theme" == "" ]; then
    echo "need to assign theme"
    exit 0
fi

dir=$dir0'/'$theme
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
cp $dir0/gC.slurm $dir
cp gainCurve.cfg $dir/input.cfg
cp $dir0/gainCurve $dir
cp $dir0/plotGainCurve.m $dir
cp -R x86_64 $dir

cp $dir0/levels.bin $dir
cp $dir0/inputTable.bin $dir

cd $dir
export theme
export newInput
sbatch --export=ALL gC.slurm
