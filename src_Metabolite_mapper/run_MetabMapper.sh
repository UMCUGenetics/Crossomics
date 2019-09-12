#!/bin/bash

if [ "$SGE_TASK_ID" = "undefined" ]; then
  echo "No array range defined (patient index)! Specify -t option"
else
  patient_number=$SGE_TASK_ID
fi

seed=$1
code_dir=$2

# code_dir="/home/cog/mkerkhofs/Crossomics_HPC/src_HPC"
# code_dir="/hpc/shared/dbg_mz/marten/Crossomics_HPC/src_HPC"
R_location="/hpc/local/CentOS7/dbg_mz/R_libs/3.6.0/bin"


for threshold in 1 2 3 4
do
	for maxrxn in 8 10 12 15 17 19
    do
	    for step in 0 1 2 3 4 5
        do
            $R_location/R --slave --no-save --no-restore --no-environ --args $threshold $maxrxn $step $patient_number $code_dir $seed < "$code_dir/GeneMetabMapper.R"
        done
    done
done
