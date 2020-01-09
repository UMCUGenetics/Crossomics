#!/bin/bash

# scripts="/Users/mkerkho7/DIMS2_repo/Crossomics/src_Metabolite_Set_Creation/Run_files/2019-08-12"
code_dir="/hpc/shared/dbg_mz/marten/Crossomics/src_Metabolite_Set_Creation/Run_files"

#for gene_id in {1..7422}
#do
#	$scripts/extend.sh $gene_id
#done

qsub -l h_rt=03:00:00 -l h_vmem=35G -N Crossomics_BuildMetabSet -t 1-7422 -o /home/cog/$USER/Crossomics/logs -e /home/cog/$USER/Crossomics/logs -m a -M "m.h.p.kerkhofs-7@umcutrecht.nl" $code_dir/extend.sh

