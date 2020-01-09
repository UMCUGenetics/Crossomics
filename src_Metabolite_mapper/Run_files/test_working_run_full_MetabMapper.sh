#!/bin/bash
# This script is to run the run_R.sh script which calls Rscript to run the GeneMetabMapper.R script

mock_genes_directory="/hpc/shared/dbg_mz/marten/Crossomics_2019_12_10/Results/Mock_genes"
code_dir="/hpc/shared/dbg_mz/marten/Crossomics_2019_12_10/src_HPC"
core_dump_dir="/hpc/dbg_mz/users/Marten/Core_Dump"
R_location="/hpc/local/CentOS7/dbg_mz/R_libs/3.6.0/bin"
thresholds="-1;1.5"
maxrxns="8"
steps="0"

mock_files=$mock_genes_directory/*.txt

set -- $mock_files

# for seed_file in $mock_files
for seed_file in $1 $2
do 

    echo "queued $seed_file"

    qsub -l h_rt=00:10:00 -l h_vmem=1G -N Crossomics_MetabMapper -wd $core_dump_dir -t 1-2 -o /home/cog/$USER/Crossomics_HPC/logs -e /home/cog/$USER/Crossomics_HPC/logs -m bea -M "m.h.p.kerkhofs-7@umcutrecht.nl" $code_dir/run_R.sh $thresholds $maxrxns $steps $code_dir $seed_file $R_location

done
