#!/usr/bin/env bash
module purge
module load Python
module load BEDTools
project_rootdir="/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
#cd ${project_rootdir}

eqtmFile=${project_rootdir}/data/eqtmZscores/allCpgs/cpg_test.txt
eqtmBedtoolsFile=${project_rootdir}/data/eqtmZscores/allCpgs/cpg_test_bedtoolsFormat.txt
awk 'BEGIN {OFS="\t"};NR>2 {print "chr"$1, $2-25, $2+25, $3}' ${eqtmFile} > ${eqtmBedtoolsFile}

python ${project_rootdir}/tools/create_bedtoolsIntersect_files.py \
    --project_rootdir ${project_rootdir} \
    --data_rootdir /groups/umcg-bios/tmp03/projects/2018-methylation/input/features/Roadmap \
    --feature_filetype consolidatedImputedGappedPeak \
    --input_type cpg_test \
    --cpg_filename cpg_test_bedtoolsFormat.txt
echo "Made all intersection files."

python ${project_rootdir}/tools/makeOverlapMatrix.py \
    --input_dirname cpg_test\
    --cpg_filename cpg_test_bedtoolsFormat.txt\
    --original_eqtm_filepath ${project_rootdir}/data/eqtmZscores/allCpgs/cpg_test.txt \
    --name cpg_test
echo "Made overlapMatrix and overlapRatio."