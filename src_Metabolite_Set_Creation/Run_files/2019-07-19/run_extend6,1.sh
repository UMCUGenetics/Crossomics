#!/bin/bash

scripts="/Users/mkerkho7/DIMS2_repo/Crossomics/src_Metabolite_Set_Creation/Run_files/2019-07-19"

for gene_id in {7374..7421}
do
	$scripts/extend.sh $gene_id
done
