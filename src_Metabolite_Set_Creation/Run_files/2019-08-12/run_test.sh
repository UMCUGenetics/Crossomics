#!/bin/bash

# scripts="/Users/mkerkho7/DIMS2_repo/Crossomics/src_Metabolite_Set_Creation/Run_files/2019-08-12"
code_dir="/hpc/shared/dbg_mz/marten/Crossomics/src_Metabolite_Set_Creation/Run_files"

#for gene_id in {1..2}
#do
#	$scripts/extend_test.sh $gene_id
#done

qsub -l h_rt=00:01:00 -l h_vmem=0.1G -N Crossomics_BuildMetabSet -t 1-3 -o /home/cog/$USER/Crossomics/logs -e /home/cog/$USER/Crossomics/logs -m bea -M "m.h.p.kerkhofs-7@umcutrecht.nl" $code_dir/extend_test.sh
