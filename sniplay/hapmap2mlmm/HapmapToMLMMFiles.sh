#!/bin/bash
hapmap=$1
map=$2
geno=$3

directory=`dirname $0`
 
perl $directory/HapmapToMLMMFiles.pl -h $hapmap -g $geno -m $map -p $directory



