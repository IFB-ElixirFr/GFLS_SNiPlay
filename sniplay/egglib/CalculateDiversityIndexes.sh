#!/bin/bash
input=$1
output=$2
log=$3

directory=`dirname $0`
 
perl $directory/CalculateDiversityIndexes.pl -i $input -o $output -d $directory >>$log 2>&1



