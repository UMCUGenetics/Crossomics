#!/bin/bash

# code_dir="/home/cog/mkerkhofs/Crossomics_HPC/src_HPC"
code_dir="/hpc/shared/dbg_mz/marten/Crossomics_2019_09_03/src_HPC"

seed=$1
echo "Seed is $seed"

qsub -l h_rt=20:00:00 -l h_vmem=20G -N Crossomics_MetabMapper -t 1-119 -o /home/cog/$USER/Crossomics_HPC/logs -e /home/cog/$USER/Crossomics_HPC/logs -m ea -M "m.h.p.kerkhofs-7@umcutrecht.nl" $code_dir/run_MetabMapper.sh $seed $code_dir
