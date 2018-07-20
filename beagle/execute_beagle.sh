#!/bin/bash


vcf_input=$1
phase=$2
impute=$3
out_prefix="out"

mkdir tmpdir$$

cp -rf $vcf_input tmpdir$$/input.vcf

beagle gt=$vcf_input phase-its=$phase impute-its=$impute out=$out_prefix

gunzip $out_prefix.vcf.gz

rm -rf tmpdir$$

