#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=8:00:00
#SBATCH -a 1-4
#SBATCH -o cK_%A_%a.out
#SBATCH --job-name=cK
#SBATCH --mail-type=END,FAIL

module purge
module load anaconda
module list

date
source activate py2
echo $nc
nc=1
v0=-74
dV=4
seed=238660
dt=0
fmt='eps'
#nc=1    # no clamp
so=0
#so=1    # single only
#dtRange = np.array([0,4,8,12,20,30,70,130,210],dtype='float')
python confirmK.py --dV=$dV --vmin=$v0 --vid=$(($SLURM_ARRAY_TASK_ID-1)) --seed=$seed --dt=$dt --fmt=$fmt --so=$so --nf=1 --nc=$nc

date
hostname
