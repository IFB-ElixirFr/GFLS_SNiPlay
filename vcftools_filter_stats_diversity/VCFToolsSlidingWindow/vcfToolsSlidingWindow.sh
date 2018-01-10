#!/bin/bash

tool_path=$(dirname $0)

filein=$1
fileout_label=$(date "+%Y%m%d%H%M%S")
window=$2
filelog=$3
fileout_taj=$4
fileout_tstv=$5
fileout_windowed=$6
fileout_snp=$7
group=$8
if [ "$group" != "none" ]
then fileout_fst=${9}
	fileout_fst_marker=${10}
	fileout_taj_combined=${11}
	fileout_windowed_combined=${12}
	fileout_tstv_combined=${13}
	fileout_snp_combined=${14}
fi


cp $filein ${filein}.vcf 
if [ "$group" != "none" ] 
then 
	if [ $(sort -k2,2 -t ';' -u $group | cut -f2 -d ';' |wc -l) -ge 3 ]
        then
                echo 'ERROR: it only works with 2 groups. There is too many groups in your file.' >&2
                exit 1
        else
                perl $tool_path/VCFToolsSlidingWindow.pl --input ${filein}.vcf --out $fileout_label --window $window --group $group
        fi


else perl $tool_path/VCFToolsSlidingWindow.pl --input ${filein}.vcf --out $fileout_label --window $window
fi
        
mv ${fileout_label}.vcftools.log $filelog 
mv ${fileout_label}.Tajima.D ${fileout_taj} 
mv ${fileout_label}.TsTv ${fileout_tstv} 
mv ${fileout_label}.windowed.pi ${fileout_windowed}
mv ${fileout_label}.snpden ${fileout_snp} 

if [ "$group" != "none" ] 
then mv ${fileout_label}.fst.txt ${fileout_fst} 
mv ${fileout_label}.fst.by_marker.genes.txt ${fileout_fst_marker} 
mv ${fileout_label}.combined.dtajima.txt ${fileout_taj_combined} 
mv ${fileout_label}.combined.pi.txt ${fileout_windowed_combined}
paste ${fileout_label}.Pop1.TsTv ${fileout_label}.Pop2.TsTv | column -s $'\t' -t | awk '{OFS="\t"; print $1,$2,$4,$8} ' | awk 'NR==1{OFS="\t"; $3="Pop1";$4="Pop2"}1' > ${fileout_tstv_combined}
paste ${fileout_label}.Pop1.snpden ${fileout_label}.Pop2.snpden | column -s $'\t' -t | awk '{OFS="\t"; print $1,$2,$4,$8} ' | awk 'NR==1{OFS="\t"; $3="Pop1";$4="Pop2"}1' > ${fileout_snp_combined}
fi

rm -f ${filein}.vcf ${filein}.vcf.*

