#!/bin/bash

tool_path=$(dirname $0)
input=$1
fileout=$2
fileout_bysample=$3
step=$4

perl $tool_path/CalculateSlidingWindowsSNPdensitiesFromVCF.pl -i $input -o $fileout -s $step

cp  $fileout.by_sample $fileout_bysample
rm $fileout.by_sample

