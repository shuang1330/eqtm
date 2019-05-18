#!/bin/bash

PROJECT_ROOTDIR=$1
FEATURE_DIR=$2
feature_filename=$3
featurename=$4
input_filepath=$5
input_name=$6

# intermediate results
UNZIPPED_FEATURE_DIR=${PROJECT_ROOTDIR}/data/temp/intersect_features/${input_name}
if [ ! -d ${UNZIPPED_FEATURE_DIR} ]; then
    mkdir -p ${UNZIPPED_FEATURE_DIR}
fi

# handle excepions
TEMP_EXCEPTION=${PROJECT_ROOTDIR}/data/temp/exception_features/${input_name}
TEMP_EMPTY=${PROJECT_ROOTDIR}/data/temp/empty_features
if [ ! -d ${TEMP_EXCEPTION} ]; then
    mkdir -p ${TEMP_EXCEPTION}
fi
if [ ! -d ${TEMP_EMPTY} ]; then
    mkdir -p ${TEMP_EXCEPTION}
fi

echo "Reading feature file ${featurename}."
# unzip feature file to temp folder
unzipped_feature_filepath=${UNZIPPED_FEATURE_DIR}/${featurename}.bed
gzip -dc ${FEATURE_DIR}/${feature_filename} > ${unzipped_feature_filepath}

# bedtools_intersect_output
if [ -s ${unzipped_feature_filepath} ]; then
  bedtools intersect -a ${unzipped_feature_filepath} -b ${input_filepath} -wa -wb > ${UNZIPPED_FEATURE_DIR}/${featurename}.bedtoolsIntersect.txt
  status=$?
  # handle exception
  if [ ! ${status} -eq 0 ]; then
    cp ${unzipped_feature_filepath} ${TEMP_EXCEPTION}/${featurename}.bed
    echo "Exeption occurs and feature ${featurename} is moved to ${TEMP_EXCEPTION}"
  fi
else
  cp ${unzipped_feature_filepath} ${TEMP_EMPTY}/${featurename}.bed
  echo "Empty feature file saved in: ${unzipped_feature_filepath}"
fi


# remove the bed file
rm ${unzipped_feature_filepath}
