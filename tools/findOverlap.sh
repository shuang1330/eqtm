#!/bin/bash

# folders where raw data is stored
export PROJECT_ROOTDIR=$1
export FEATURE_DIR=$2
export UNZIPPED_FEATURE_DIR=$3

# folders for intermediate results
export cpgname=$4
export cpg_filepath=$5

# names and filenames
export featurename=$6
export feature_filename=$7
export feature_type=$8

# handle excepions
export TEMP_EXCEPTION=$PROJECT_ROOTDIR/data/temp/exception_features/$feature_type

echo "Reading feature file $featurename."

# unzip feature file to temp folder
export unzipped_feature_filepath=$UNZIPPED_FEATURE_DIR/$featurename.bed
export feature_BedToolFormat_filepath=$UNZIPPED_FEATURE_DIR/$featurename.bedtoolsFormat.bed
gzip -dc $FEATURE_DIR/$feature_filename > $unzipped_feature_filepath

# bedtools_intersect_output
if [ -s $unzipped_feature_filepath ]; then
  bedtools intersect -a $unzipped_feature_filepath -b $cpg_filepath -wa -wb > $UNZIPPED_FEATURE_DIR/$featurename.bedtoolsIntersect.txt
  export status=$?

  if [ ! -d $TEMP_EXCEPTION ]; then
    mkdir -p $TEMP_EXCEPTION
  fi

  # handle exception
  if [ ! $status -eq 0 ]; then
    cp $unzipped_feature_filepath $TEMP_EXCEPTION/$featurename.bed
    echo "Exeption occurs and feature $featurename is moved to $TEMP_EXCEPTION"
  fi
else
  echo "Empty feature file: $unzipped_feature_filepath"
fi


# remove the bed file
rm $unzipped_feature_filepath
