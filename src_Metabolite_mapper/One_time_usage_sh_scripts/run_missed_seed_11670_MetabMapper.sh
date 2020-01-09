#!/bin/bash
# This script is to run the run_R.sh script which calls Rscript to run the GeneMetabMapper.R script

#mock_genes_directory="/hpc/shared/dbg_mz/marten/Crossomics_2019_11_01/Results/Mock_genes"
#mock_genes_directory="/hpc/shared/dbg_mz/marten/Crossomics_2019_11_01/Results/Mock_genes_missed/"
code_dir="/hpc/shared/dbg_mz/marten/Crossomics_2019_12_04/src_HPC"
R_location="/hpc/local/CentOS7/dbg_mz/R_libs/3.6.0/bin"
thresholds="-1;1.5,-1.5;2,-3;3,-5;5"
maxrxns="8,10,12,15,17,19"
steps="0,1,2,3,4,5"

mock_file_begin="mock_genes200_seed"
mock_file_end=".txt"
#for patient_number in {1..112}
patient_number=53
#do
    seed=11670
    seed_file=$mock_file_begin$seed$mock_file_end
    echo queed patient: $patient_number, seed: $seed_file
    qsub -l h_rt=00:40:00 -l h_vmem=1G -N Mock_genes_redo -o /home/cog/$USER/Crossomics_HPC/logs -e /home/cog/$USER/Crossomics_HPC/logs -m bea -M "m.h.p.kerkhofs-7@umcutrecht.nl" $code_dir/run_R_missed_seeds.sh $thresholds $maxrxns $steps $code_dir $seed_file $R_location $patient_number
#done
