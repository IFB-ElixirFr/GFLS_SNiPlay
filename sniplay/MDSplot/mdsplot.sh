#!/bin/bash

tool_path=$(dirname $0)
ped=$1
map=$2
fileout_label=$(date "+%Y%m%d%H%M%S")
fileout_matrix=$3
fileout_plot=$4
fileout_log=$5

rsync -a $ped input.ped 
rsync -a $map input.map

perl $tool_path/MDSbasedOnIBSmatrix.pl --in input --out $fileout_label

rm -f input.ped input.map

cp $fileout_label.ibs_matrix.txt $fileout_matrix
cp $fileout_label.mds_plot.txt $fileout_plot
cp input.plink.log $fileout_log


rm -f $fileout_label.ibs_matrix.txt $fileout_label.mds_plot.txt input.plink.log
