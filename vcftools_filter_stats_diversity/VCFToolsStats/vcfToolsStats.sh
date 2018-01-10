#!/bin/bash

tool_path=$(dirname $0)

filein=$1
fileout_label=$(date "+%Y%m%d%H%M%S")
fileout_het=$2
fileout_imiss=$3
fileout_sum=$4
filelog=$5



perl $tool_path/VCFToolsStats.pl --input $filein --out $fileout_label

cp  $fileout_label.het $fileout_het ; rm $fileout_label.het
cp  $fileout_label.imiss $fileout_imiss ; rm $fileout_label.imiss
cp  $fileout_label.TsTv.summary $fileout_sum ; rm $fileout_label.TsTv.summary

cp vcftools.log $filelog
rm vcftools.log
