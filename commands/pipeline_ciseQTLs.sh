#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --job-name=ciseqtl
#SBATCH --output=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/ciseqtl.out
#SBATCH --error=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/ciseqtl.err

module purge
module load Python
module load BEDTools
project_rootdir="/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
#cd ${project_rootdir}

eqtmFile=${project_rootdir}/data/eqtmZscores/allCpgs/cpg_ciseQTL_topSNPs_TSS_TES_for_longest_transcript_20180531.txt
eqtmBedtoolsFile=${project_rootdir}/data/eqtmZscores/allCpgs/cpg_ciseQTL_topSNPs_TSS_TES_for_longest_transcript_20180531_bedtoolsFormat.txt
awk 'BEGIN {OFS="\t"};NR>2 {print "chr"$4, $5-25, $5+25, $3}' ${eqtmFile} > ${eqtmBedtoolsFile}

python ${project_rootdir}/tools/create_bedtoolsIntersect_files.py \
    --project_rootdir ${project_rootdir} \
    --data_rootdir /groups/umcg-bios/tmp03/projects/2018-methylation/input/features/Roadmap \
    --feature_filetype consolidatedImputedGappedPeak \
    --input_type cpg_ciseQTL_topSNPs_TSS_TES_for_longest_transcript_20180531 \
    --cpg_filename cpg_ciseQTL_topSNPs_TSS_TES_for_longest_transcript_20180531_bedtoolsFormat.txt
echo "Made all intersection files."

##todo: the following steps
#
## create the *_with_header
#
#python ${project_rootdir}/tools/makeOverlapMatrix.py \
#    --input_dirname cpg_genetic_risk_factors_cleaned_traits_standardized\
#    --cpg_filename cpg_genetic_risk_factors_cleaned_traits_standardized_bedtoolsFormat.txt\
#    --original_eqtm_filepath ${project_rootdir}/data/eqtmZscores/allCpgs/cpg_genetic_risk_factors_cleaned_traits_standardized_bedtoolsFormat_with_header.txt \
#    --name cpg_genetic_risk_factors_cleaned_traits_standardized
#echo "Made overlapMatrix and overlapRatio."
