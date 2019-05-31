#! /bin/bash

# load(paste(src, "HGNC2Uniprot.RData", sep="/"))
# dim(HGNC2Uniprot) => 19266

scripts="/hpc/shared/dbg_mz/marcel/Crossomics_Build_Metabolite_Set/src"

qsub -l h_rt=08:00:00 -l h_vmem=16G -t 1-22968:1 -tc 50 $scripts/extend.sh
