#!/bin/bash
#####################################################################
scripts="/hpc/shared/dbg_mz/marcel/Crossomics_Build_Metabolite_Set/src"
inpdir="/hpc/shared/dbg_mz/marcel/Crossomics_Build_Metabolite_Set/db"
outdir="/hpc/shared/dbg_mz/marcel/Crossomics_Build_Metabolite_Set/results"

find $inpdir -iname "*.RData" | while read rdata;
 do
   echo "Bingo submitting $rdata";
   qsub -l h_rt=08:00:00 -l h_vmem=16G -N "findSurroundings" $scripts/findSurroundings.sh $rdata $outdir $scripts
 done