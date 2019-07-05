#!/bin/bash

if [ "$SGE_TASK_ID" = "undefined" ]; then
  echo "No array range defined! Specify -t option"
else
  gene_id=$SGE_TASK_ID
fi

scripts="/hpc/shared/dbg_mz/marcel/Crossomics_Build_Metabolite_Set/src"
outdir="/hpc/shared/dbg_mz/marcel/Crossomics_Build_Metabolite_Set/results"

echo "Run set $gene_id in R"

module load R
R --slave --no-save --no-restore --no-environ --args $gene_id $outdir $scripts < "$scripts/extend.R"
