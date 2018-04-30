#!/bin/bash
export PROJECT_ROOTDIR="/home/shuang/projects/eqtm"
export FEATURES=$PROJECT_ROOTDIR/data/RoadmapFeatureSamples
export CPGS=$PROJECT_ROOTDIR/data/eqtmZscores

# for intermediate results
export TEMP_CPG=$PROJECT_ROOTDIR/data/temp/cpg
export TEMP_FEATURE=$PROJECT_ROOTDIR/data/temp/features
export TEMP_OUTPUT=$PROJECT_ROOTDIR/data/temp/output

# get the feature file
export featurename=$1 #c needs to be provided
# get the cpg file
export cpgname=$2
export cpgfilename=$TEMP_CPG/$cpgname.tsv
# check the argument 1
if [[ -n "$featurename" ]];then
  echo "Reading feature file $1."
else
  echo "Argument Error, no feature file given."
  exit
fi
# check argument 2
if [[ -n "$cpgname" ]];then
  echo "Reading cpg file $2."
else
  echo "Argument Error, no cpgname file given."
  exit
fi

# unzip feature file to temp folder
name=$(basename "$featurename" .imputed.narrowPeak.bed.nPk.gz)
export savefeaturename=$TEMP_FEATURE/$featurename.bed
echo $savefeaturename
gzip -dc $FEATURES/$featurename.imputed.narrowPeak.bed.nPk.gz > $savefeaturename

# bedtools_intersect_output
bedtools intersect -sorted -a $savefeaturename -b $cpgfilename -wa -wb > $TEMP_OUTPUT/$featurename.bedtoolsIntersect.txt

# remove the unzipped file
rm $savefeaturename
