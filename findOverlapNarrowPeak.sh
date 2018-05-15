#!/bin/bash

# folders where raw data is stored
export PROJECT_ROOTDIR=$1
export FEATURES=$2
export CPGS=$3

# folders for intermediate results
export TEMP_CPG=$4
export TEMP_FEATURE=$5
export TEMP_OUTPUT=$6

# file names
export featurename=$7 # feature file name in temp_featureFolder
export cpgname=$8 # cpg file name
export cpgfilename=$TEMP_CPG/$cpgname.tsv

echo "Reading feature file $featurename."
echo "Reading cpg file $cpgname."

# unzip feature file to temp folder
name=$(basename "$featurename" .imputed.narrowPeak.bed.nPk.gz)
export savefeaturename=$TEMP_FEATURE/$featurename.bed
# echo "processing feature:"$featurename
gzip -dc $FEATURES/$featurename.imputed.narrowPeak.bed.nPk.gz > $savefeaturename
echo "Saved the intermediate results to "$savefeaturename

# bedtools_intersect_output
bedtools intersect -sorted -a $savefeaturename -b $cpgfilename -wa -wb > $TEMP_FEATURE/$featurename.bedtoolsIntersect.txt
echo "Saved results to "$TEMP_FEATURE/$featurename.bedtoolsIntersect.txt

# return exit code


# remove the bed file
rm $savefeaturename
