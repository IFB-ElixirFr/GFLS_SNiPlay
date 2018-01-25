#!/bin/bash
hapmap=$1
trait=$2
out_hapmap=$3
out_trait=$4
stats=$5

directory=`dirname $0`
mkdir tmpdir$$
#cp -rf $input tmpdir$$/input
 
perl $directory/CheckGWASInputs.pl -h $hapmap -t $trait -o tmpdir$$/out >>$stats 2>&1

ls
ls tmpdir$$/

mv tmpdir$$/out.hapmap $out_hapmap
mv tmpdir$$/out.trait $out_trait


