#!/bin/bash
vcf=$1
genome=$2
gff=$3
output=$4
html=$5
log=$6

directory=`dirname $0`
 
/usr/bin/perl $directory/SnpEff.pl -i $vcf -f $genome -g $gff -o $output -h $html >>$log 2>&1



