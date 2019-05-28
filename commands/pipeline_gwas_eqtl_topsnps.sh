#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --job-name=gwas_eqtl_histoneAnnotation
#SBATCH --output=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/out/gwas_eqtl_histoneAnnotation.out
#SBATCH --error=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/out/gwas_eqtl_histoneAnnotation.err
#SBATCH --mem=50gb

module purge
module load Python
module load BEDTools
project_rootdir="/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
#cd ${project_rootdir}

eqtmFile=/groups/umcg-wijmenga/tmp03/projects/eQTLGen/interpretation/GWAS_SNPs_UKB_epigenetic_marks/annotate_GWAS_SNPs_UKB/results/combined_significant_stringent_SNPs_location_2019-05-16.txt
eqtmBedtoolsFile=${project_rootdir}/data/eqtmZscores/allCpgs/combined_significant_stringent_SNPs_location_2019-05-16_bedtoolFormat.txt
#awk 'BEGIN {OFS="\t"};NR>2 {print "chr"$2, $3-25, $3+25, $1}' ${eqtmFile} > ${eqtmBedtoolsFile}

#python ${project_rootdir}/tools/create_bedtoolsIntersect_files.py \
#    --project_rootdir ${project_rootdir} \
#    --data_rootdir /groups/umcg-bios/tmp03/projects/2018-methylation/input/features/Roadmap \
#    --feature_filetype consolidatedImputedGappedPeak \
#    --input_type combined_significant_stringent_SNPs_location_2019-05-16 \
#    --cpg_filename combined_significant_stringent_SNPs_location_2019-05-16_bedtoolFormat.txt
#echo "Made all intersection files."

python ${project_rootdir}/tools/makeOverlapMatrix_significant.py \
    --input_dirname combined_significant_stringent_SNPs_location_2019-05-16 \
    --cpg_filename combined_significant_stringent_SNPs_location_2019-05-16_bedtoolFormat.txt\
    --original_eqtm_filepath ${eqtmFile} \
    --name combined_significant_stringent_SNPs_location_2019-05-16
echo "Made overlapRatio."