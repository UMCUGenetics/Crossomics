#!/bin/bash
# File to exclusively run GeneMetabMapper.R with all arguments that it needs. It needs to be called via a qsub call which holds the information for all of these parameters (see run_full_MetabMapper.sh)

patient_number=$SGE_TASK_ID
thresholds=$1
maxrxns=$2
steps=$3
code_dir=$4
seed_file=$5 
R_location=$6

$R_location/Rscript $code_dir/GeneMetabMapper.R $patient_number $thresholds $maxrxns $steps $code_dir $seed_file $R_location

rm -f ./core.*