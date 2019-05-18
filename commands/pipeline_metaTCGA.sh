#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --job-name=meta_TCGA
#SBATCH --output=meta_TCGA.out

module purge
module load Python
module load BEDTools
project_rootdir="/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
cd ${project_rootdir}

eqtmFile=${project_rootdir}/data/eqtmZscores/allCpgs/meta_analysis_tcga_flipped.txt
eqtmBedtoolsFile=${project_rootdir}/data/eqtmZscores/allCpgs/metaTCGA_bedtoolsFormat.txt
awk 'BEGIN {OFS="\t"};NR>2 {print "chr"$3, $4-25, $4+25, $2}' ${eqtmFile} > ${eqtmBedtoolsFile}

python ./tools/create_bedtoolsIntersect_files.py \
    --project_rootdir ${project_rootdir} \
    --data_rootdir /groups/umcg-bios/tmp03/projects/2018-methylation/input/features/Roadmap \
    --feature_filetype consolidatedImputedGappedPeak \
    --input_type cpg_metaTCGA \
    --cpg_filename metaTCGA_bedtoolsFormat.txt
echo "Made all intersection files."

python ${project_rootdir}/tools/makeOverlapMatrix.py \
    --input_dirname cpg_metaTCGA\
    --cpg_filename metaTCGA_bedtoolsFormat.txt\
    --original_eqtm_filepath ${project_rootdir}/data/eqtmZscores/allCpgs/metaTCGA_bedtoolsFormat_with_header.txt \
    --name cpg_metaTCGA
echo "Made overlapMatrix and overlapRatio."