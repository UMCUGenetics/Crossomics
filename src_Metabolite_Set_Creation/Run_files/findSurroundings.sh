#!/bin/bash
file=$1
outdir=$2
scripts=$3

echo "### Inputs findSurroundings.sh #####################################################"
echo "	file:	$file"
echo "	outdir:	$outdir"
echo "	scripts:	$scripts"
echo "#############################################################################"

echo "Run file $file in R"

module load R
R --slave --no-save --no-restore --no-environ --args $file $outdir $scripts < "$scripts/findSurroundings.R"