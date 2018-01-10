#!/bin/bash

tool_path=$(dirname $0)

filein=$1
fileout_label=$(date "+%Y%m%d%H%M%S")
fileout=$2
filelog=$3
export=$4
frequency=$5
max_freq=$6
allow_missing=$7
nb_alleles_min=$8
nb_alleles_max=9
type=${10}
bound_start=${11}
bound_end=${12}


if [ "${13}" != "None" ]
then samples="--samples ${13}"
fi

if [ "${14}" != "None" ]
then chromosomes="--chromosomes ${14}"
fi

if [ "$bound_start" -gt "$bound_end" ]
then tmp=$bound_start ; bound_start=$bound_end ; bound_end=$tmp ; echo "Warning : Lower bound must be lower than greater bound!" >&2
fi

if [ "$nb_alleles_min" -gt "$nb_alleles_max" ]
then tmp=$nb_alleles_min ; nb_alleles_min=$nb_alleles_max ; nb_alleles_max=$tmp ; echo "Warning : Minimum number of alleles must be lower than maximum number of allele!" >&2
fi

perl $tool_path/VCFToolsFilter.pl --input $filein --out $fileout_label --export $export --frequency $frequency --max_freq $max_freq --allow_missing $allow_missing --nb_alleles $nb_alleles_min','$nb_alleles_max --type $type --bounds $bound_start','$bound_end $samples $chromosomes

if [ "$export" = "VCF" ]
then cp  $fileout_label.recode.vcf $fileout ; rm $fileout_label.recode.vcf
elif [ "$export" = "freq" ]
then cp  $fileout_label.frq $fileout ; rm $fileout_label.frq
else cp  $fileout_label.ped $fileout; cp $fileout_label.map ${15} ; rm $fileout_label.ped $fileout_label.map
fi

cp vcftools.log $filelog
rm vcftools.log
