#!/bin/bash

gene_id=$1

scripts="/Users/mkerkho7/DIMS2_repo/Crossomics/src_Metabolite_Set_Creation"

echo "Run set $gene_id in R"
echo "$scripts/new_extend.R"

# module load R
# R --slave --no-save --no-restore --no-environ --args $gene_id< "$scripts/new_extend.R"
