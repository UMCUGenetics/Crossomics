#!/bin/bash
# This script is to run the run_R.sh script which calls Rscript to run the GeneMetabMapper.R script

mock_genes_directory="/hpc/shared/dbg_mz/marten/Crossomics_2019_11_01/Results/Mock_genes"
code_dir="/hpc/shared/dbg_mz/marten/Crossomics_2019_11_01/src_HPC"
R_location="/hpc/local/CentOS7/dbg_mz/R_libs/3.6.0/bin"
thresholds="-1;1.5,-1.5;2,-3;3,-5;5"
maxrxns="8,10,12,15,17,19"
steps="0,1,2,3,4,5"

mock_files=$mock_genes_directory/*.txt

# set -- $mock_files # for testing

for seed_file in $mock_files
# for seed_file in $1 $2 # for testing
do 

    echo "queued $seed_file"

    qsub -l h_rt=04:30:00 -l h_vmem=0.5G -N Crossomics_MetabMapper -t 1-112 -o /home/cog/$USER/Crossomics_HPC/logs -e /home/cog/$USER/Crossomics_HPC/logs -m a -M "m.h.p.kerkhofs-7@umcutrecht.nl" $code_dir/run_R.sh $thresholds $maxrxns $steps $code_dir $seed_file $R_location

done
