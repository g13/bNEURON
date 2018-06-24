#!/bin/bash
n128=1
dir0='gainCurve'
getDendV=true

fdr=$1
newInput=$2
ext=$3
input=true
level=false

if [ "$fdr" == "" ]; then
    echo "need to assign fdr"
    exit 0
fi

dir=$dir0'/'$fdr
if [ -d "$dir" ]; then
    rm -r $dir/*
else
    mkdir $dir
fi

if [ "$n128" == "0" ]; then
    cp base/n128.py_pas $dir/n128.py
fi
if [ "$n128" == "1" ]; then
    cp base/n128.py_active $dir/n128.py
fi
if [ "$n128" == "2" ]; then
    cp base/n128.py_simple $dir/n128.py
fi

if [ "$newInput" == "1" ]; then
    cp gainCurve.cfg $dir0
    cp setInput.m $dir0
    cd $dir0
    matlab -nosplash -nodisplay -r "setInput($input,$level);exit"
    cd ..
else 
    echo no new input
fi
cp base/neuroAlter.py $dir
cp $dir0/one.slurm $dir
cp $dir0/plot.slurm $dir
cp gainCurve.cfg $dir/input.cfg
cp $dir0/one $dir
cp $dir0/plotGainCurve.m $dir
cp $dir0/read_cfg.m $dir
cp -r x86_64 $dir

cp $dir0/levels.bin $dir
cp $dir0/inputTable.bin $dir


cd $dir
export ext
export getDendV
sbatch --export=ALL one.slurm
