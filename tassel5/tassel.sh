#!/bin/bash

analyseType=$1;
out1=$2;
out2=$3;
out3=$4;
log1=$5;
galaxyOutDir=$6;

version=$(java -version 2>&1 | grep version)
if [[ ! $version =~ 1.8 ]]; then
    echo "Java found: $version. Tassel 5.0 requires java 1.8..." >&2
    exit 1
fi

mkdir $galaxyOutDir

# Suppression des 6 premiers arguments de la liste des arguments $@
shift; shift; shift; shift; shift; shift;


if [[ $analyseType == glm ]] 
then
 run_pipeline.pl $* >> $log1 2>&1
 mv "$galaxyOutDir/TASSELGLM1.txt" $out1
 if [ $? -gt 0 ] ; then exit 1 ; fi ##check the run of commands
 mv "$galaxyOutDir/TASSELGLM2.txt" $out2
fi

if [[ $analyseType == mlm ]] 
then
 run_pipeline.pl $* >> $log1 2>&1
 mv "$galaxyOutDir/TASSELMLM1.txt" $out1
 if [ $? -gt 0 ] ; then exit 1 ; fi ##check the run of commands
 mv "$galaxyOutDir/TASSELMLM2.txt" $out2
 mv "$galaxyOutDir/TASSELMLM3.txt" $out3
fi

if [[ $analyseType == ld ]] 
then
 run_pipeline.pl $* >> $log1 2>&1
fi


if [[ $analyseType == ck ]]
then
 run_pipeline.pl $* >> $log1 2>&1
 mv "$galaxyOutDir/kinship.txt" $out1
 if [ $? -gt 0 ] ; then exit 1 ; fi ##check the run of commands
fi
