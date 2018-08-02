#!/bin/bash

# Config #######
maxsizeseq=500
maxnumseq=200
################

tool_path=$(dirname $0)

filein=$1
fileout=$2
dotfile=$3
logfile=$4
filein2=$5
groups=$6

nbline=$(sed -n '$=' $filein)
let "nbseq = $nbline / 2"
seq=$(sed -n 2p $filein)
sizeseq=${#seq}

if [ $nbseq -lt $maxnumseq ]
then
        if [ $sizeseq -lt $maxsizeseq ]
        then
	        perl $tool_path/Haplophyle.pl --input $filein --groups $groups --stats $filein2 --dot $dotfile --out $fileout --tool_path $tool_path >>$logfile 2>&1
        else
                echo "Sequence size: $sizeseq"
                echo "Input Sequences bust have a length < $maxsizeseq"
                exit 1
        fi
else
        echo "$nbseq sequences in the file"
        echo "Input file must have less than $maxnumseq sequences"
        exit 1
fi
