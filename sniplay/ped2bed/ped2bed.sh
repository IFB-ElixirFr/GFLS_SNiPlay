#!/bin/bash
ped=$1
map=$2
bed=$3
fam=$4
bim=$5
logs=$6

directory=`dirname $0`
mkdir tmpdir$$
cp -rf $ped tmpdir$$/input.ped
cp -rf $map tmpdir$$/input.map
 
plink --file tmpdir$$/input --out tmpdir$$/out --make-bed --noweb >>$logs 2>&1

mv tmpdir$$/out.bed $bed
mv tmpdir$$/out.fam $fam
mv tmpdir$$/out.bim $bim


