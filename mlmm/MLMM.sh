#!/bin/bash
geno=$1
map=$2
pheno=$3
steps=$4
method=$5
output=$6
pdf=$7
kinship=$8
rss=$9
step_table=${10}
log=${11}

directory=`dirname $0`
mkdir tmpdir$$
cp -rf $geno tmpdir$$/geno
cp -rf $map tmpdir$$/map
cp -rf $pheno tmpdir$$/pheno
 
perl $directory/MLMM.pl -g tmpdir$$/geno -i tmpdir$$/map -p tmpdir$$/pheno -s $steps -m $method -o tmpdir$$/output -d $directory/source_library >>$log 2>&1

mv tmpdir$$/output.pdf $pdf
mv tmpdir$$/output.kinship $kinship
mv tmpdir$$/output.res_asso $output
mv tmpdir$$/output.rss $rss
mv tmpdir$$/output.steptable $step_table


