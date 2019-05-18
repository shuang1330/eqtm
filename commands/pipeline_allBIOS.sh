#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --job-name=calculateOverlapRatio_for_biosComplete
#SBATCH --output=bioComplete.out
#SBATCH --mem=20gb

module purge
module load Python
module load BEDTools

# get all cpg sites
#eqtmfile=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/ORIGIN/full_dataset.txt
#eqtmbedtoolsFile1=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/all_cpgs/cpg_biosComplete_bedtoolsFormat_withDuplicates.txt
#awk 'BEGIN {OFS="\t"};NR>2 {print "chr"$3, $4-25, $4+25, $2}' ${eqtmfile} > ${eqtmbedtoolsFile1}
#eqtmbedtoolsFile=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/all_cpgs/cpg_biosComplete_bedtoolsFormat.txt
#awk '!x[$0]++' ${eqtmbedtoolsFile1} > ${eqtmbedtoolsFile}
#rm ${eqtmbedtoolsFile1}
project_rootdir="/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
#cd ${project_rootdir}
echo "Start processing."

python ${project_rootdir}/tools/create_bedtoolsIntersect_files.py \
    --project_rootdir ${project_rootdir} \
    --data_rootdir /groups/umcg-bios/tmp03/projects/2018-methylation/input/features/Roadmap \
    --feature_filetype consolidatedImputedGappedPeak \
    --input_type cpg_biosComplete \
    --cpg_filename cpg_biosComplete_bedtoolsFormat.txt
echo "Made all intersection files."


withheader=${project_rootdir}/data/eqtmZscores/allCpgs/cpg_biosComplete_bedtoolsFormat_with_header.txt
cat ${project_rootdir}/data/eqtmZscores/allCpgs/eqtm_fibroblast_bedtoolsFormat_with_header.txt | head -1 > ${withheader}
cat ${eqtmbedtoolsFile} >> ${withheader}

python ${project_rootdir}/tools/makeOverlapMatrix.py \
    --input_dirname cpg_biosComplete\
    --cpg_filename cpg_biosComplete_bedtoolsFormat.txt\
    --original_eqtm_filepath ${project_rootdir}/data/eqtmZscores/allCpgs/cpg_biosComplete_bedtoolsFormat_with_header.txt \
    --name cpg_biosComplete
echo "Made overlapMatrix and overlapRatio."


#python ${project_rootdir}/tools/add_all_features_for_bios_data.py --name meta_analysis_tcga_flipped
echo "Added all features."

echo "Time to see how bad the model performs..."
#python ${project_rootdir}/tools/examine_model_performance.py
