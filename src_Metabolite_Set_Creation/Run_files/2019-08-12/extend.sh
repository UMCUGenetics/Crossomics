#!/bin/bash

# gene_id=$1
gene_id=$SGE_TASK_ID

code_dir="/hpc/shared/dbg_mz/marten/Crossomics/src_Metabolite_Set_Creation"
outdir="/hpc/shared/dbg_mz/marten/Crossomics/Data/2019-08-12"
Rlocation="/hpc/local/CentOS7/dbg_mz/R_libs/3.6.0"



echo "Run set $gene_id in R"

# module load R
$Rlocation/bin/R --slave --no-save --no-restore --no-environ --args $gene_id $outdir $code_dir $Rlocation< "$code_dir/extend.R"
