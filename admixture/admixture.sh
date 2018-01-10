#!/bin/bash
bed=$1
fam=$2
bim=$3
outputs=$4
logs=$5
best_k_output=$6
best_k_logfile=$7
kmin=$8
kmax=$9
groups=${10}
threshold_group=${11}

directory=`dirname $0`
mkdir tmpdir$$
cp -rf $bed tmpdir$$/input.bed
cp -rf $fam tmpdir$$/input.fam
cp -rf $bim tmpdir$$/input.bim

 
perl $directory/Admixture.pl -i tmpdir$$/input -o $outputs -k $kmin -m $kmax -d tmpdir$$ -t $threshold_group

mv tmpdir$$/output $best_k_output
mv tmpdir$$/log $best_k_logfile
mv tmpdir$$/outputs.Q $outputs
mv tmpdir$$/logs $logs
mv tmpdir$$/groups $groups


