#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --job-name=pipeline_methy_anno
#SBATCH --output=pipeline_methy_anno.out
#SBATCH --error=pipeline_methy_anno.err
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --mem=50gb

module purge
module load Python
module load BEDTools


project_rootdir="/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
echo "Start processing."

methy_probecenter="/groups/umcg-gcc/tmp03/umcg-sli/tcga/anno_files/methy_probecenter.txt"
methy_bedtools=${project_rootdir}/data/eqtmZscores/allCpgs/cpg_methy_anno_bedtoolsFormat.txt
#awk 'BEGIN {OFS="\t"};NR>2 {print "chr"$4, $5-25, $5+25, $2}' ${methy_probecenter} > ${methy_bedtools}


python ${project_rootdir}/tools/create_bedtoolsIntersect_files.py \
    --project_rootdir ${project_rootdir} \
    --data_rootdir /groups/umcg-bios/tmp03/projects/2018-methylation/input/features/Roadmap \
    --feature_filetype consolidatedImputedGappedPeak \
    --input_type cpg_methy_anno \
    --cpg_filename cpg_methy_anno_bedtoolsFormat.txt
echo "Made all intersection files."


withheader=${project_rootdir}/data/eqtmZscores/allCpgs/cpg_methy_anno_bedtoolsFormat_with_header.txt
cat ${project_rootdir}/data/eqtmZscores/allCpgs/eqtm_fibroblast_bedtoolsFormat_with_header.txt | head -1 > ${withheader}
cat ${methy_bedtools} >> ${withheader}

python ${project_rootdir}/tools/makeOverlapMatrix.py \
    --input_dirname cpg_methy_anno\
    --cpg_filename cpg_methy_anno_bedtoolsFormat.txt\
    --original_eqtm_filepath ${project_rootdir}/data/eqtmZscores/allCpgs/cpg_methy_anno_bedtoolsFormat_with_header.txt \
    --name cpg_methy_anno
echo "Made overlapMatrix and overlapRatio."


#python ${project_rootdir}/tools/add_all_features_for_bios_data.py --name meta_analysis_tcga_flipped
#echo "Added all features."

#echo "Time to see how bad the model performs..."
#python ${project_rootdir}/tools/examine_model_performance.py

