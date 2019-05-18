#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --job-name=combined2_nonBlood_blood
#SBATCH --output=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/out/combined2_nonBlood_blood.out
#SBATCH --error=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/out/combined2_nonBlood_blood.err
#SBATCH --mem=50gb

module purge
module load Python
module load BEDTools
project_rootdir="/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
#cd ${project_rootdir}

eqtmFile=${project_rootdir}/data/eqtmZscores/allCpgs/cpg_combined_significant_SNPs_location_2019-02-21.txt
eqtmBedtoolsFile=${project_rootdir}/data/eqtmZscores/allCpgs/cpg_combined_significant_SNPs_location_2019-02-21_bedtoolsFormat.txt
awk 'BEGIN {OFS="\t"};NR>2 {print "chr"$2, $3-25, $3+25, $1}' ${eqtmFile} > ${eqtmBedtoolsFile}

# python ${project_rootdir}/tools/create_bedtoolsIntersect_files.py \
#    --project_rootdir ${project_rootdir} \
#    --data_rootdir /groups/umcg-bios/tmp03/projects/2018-methylation/input/features/Roadmap \
#    --feature_filetype consolidatedImputedGappedPeak \
#    --input_type cpg_combined_significant_SNPs_location_2019-02-21 \
#    --cpg_filename cpg_combined_significant_SNPs_location_2019-02-21_bedtoolsFormat.txt
#echo "Made all intersection files."

python ${project_rootdir}/tools/makeOverlapmatrix_selectCellType.py \
    --input_dirname cpg_combined_significant_SNPs_location_2019-02-21\
    --cpg_filename cpg_combined_significant_SNPs_location_2019-02-21_bedtoolsFormat.txt\
    --original_eqtm_filepath ${project_rootdir}/data/eqtmZscores/allCpgs/cpg_combined_significant_SNPs_location_2019-02-21.txt \
    --name cpg_combined_significant_SNPs_location_2019-02-21 \
    --cellType blood \
    --savepath ${project_rootdir}/data/eqtmZscores/allCpgs/cpg_combined_significant_SNPs_location_2019-02-21_OverlapMatrix_blood.txt
echo "Made overlapRatio table for blood."

python ${project_rootdir}/tools/makeOverlapmatrix_selectCellType.py \
    --input_dirname cpg_combined_significant_SNPs_location_2019-02-21\
    --cpg_filename cpg_combined_significant_SNPs_location_2019-02-21_bedtoolsFormat.txt\
    --original_eqtm_filepath ${project_rootdir}/data/eqtmZscores/allCpgs/cpg_combined_significant_SNPs_location_2019-02-21.txt \
    --name cpg_combined_significant_SNPs_location_2019-02-21 \
    --cellType non_blood \
    --savepath ${project_rootdir}/data/eqtmZscores/allCpgs/cpg_combined_significant_SNPs_location_2019-02-21_OverlapMatrix_nonBlood.txt
echo "Made overlapRatio table for non blood."

python ${project_rootdir}/tools/link_back_combined.py
