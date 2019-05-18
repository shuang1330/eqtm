#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --job-name=variants
#SBATCH --output=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/variants.out
#SBATCH --error=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/variants.err
#SBATCH --mem=60gb
#SBATCH --nodes=1
#SBATCH --get-user-env=L

module purge
module load Python
module load BEDTools
project_rootdir="/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
#cd ${project_rootdir}

eqtmFile=/groups/umcg-gcc/tmp03/umcg-sli/variant_prioritization/data/full_data_split/full_data_splited.txt
eqtmBedtoolsFile=${project_rootdir}/data/eqtmZscores/allCpgs/cpg_variants_full_data_split_bedtoolsFormat.txt

##todo: the following steps
#
## create the *_with_header
#

#awk -F',' 'BEGIN {OFS="\t"};NR>1 {print "chr"$3, $4-25, $4+25, $1}' ${eqtmFile} > ${eqtmBedtoolsFile}

python ${project_rootdir}/tools/create_bedtoolsIntersect_files.py \
    --project_rootdir ${project_rootdir} \
    --data_rootdir /groups/umcg-bios/tmp03/projects/2018-methylation/input/features/Roadmap \
    --feature_filetype consolidatedImputedGappedPeak \
    --input_type cpg_variants_full_data_split \
    --cpg_filename cpg_variants_full_data_split_bedtoolsFormat.txt
echo "Made all intersection files."

python ${project_rootdir}/tools/makeOverlapMatrix.py \
    --input_dirname cpg_variants_full_data_split \
    --cpg_filename cpg_variants_full_data_split_bedtoolsFormat.txt \
    --original_eqtm_filepath ${project_rootdir}/data/eqtmZscores/allCpgs/cpg_variants_full_data_split_bedtoolsFormat_with_header.txt \
    --name cpg_ciseQTL_topSNPs_TSS_TES_for_longest_transcript_20180531
echo "Made overlapMatrix and overlapRatio."
