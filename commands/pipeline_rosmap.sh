#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --job-name=rosmap_histoneAnnotation
#SBATCH --output=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/out/rosmap_histoneAnnotation.out
#SBATCH --error=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/out/rosmap_histoneAnnotation.err
#SBATCH --mem=50gb

module purge
module load Python
module load BEDTools
project_rootdir="/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
#cd ${project_rootdir}

eqtmFile=/groups/umcg-biogen/tmp03/umcg-sli/ciseqtm/eQTLsFDR0.05.txt
eqtmBedtoolsFile=${project_rootdir}/data/eqtmZscores/allCpgs/rosmap_significantCpgSites_bedtoolsFormat.txt
#python /groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/tools/get_CpgSite_Position.py \
#--input_path ${eqtmFile} \
#--bedtoolFormat_savepath ${eqtmBedtoolsFile}
##awk 'BEGIN {OFS="\t"};NR>2 {print "chr"$3, $4-25, $4+25, $2}' ${eqtmFile} > ${eqtmBedtoolsFile}

python ${project_rootdir}/tools/create_bedtoolsIntersect_files.py \
    --project_rootdir ${project_rootdir} \
    --data_rootdir /groups/umcg-bios/tmp03/projects/2018-methylation/input/features/Roadmap \
    --feature_filetype consolidatedImputedGappedPeak \
    --input_type rosmap_significantCpgSites \
    --cpg_filename rosmap_significantCpgSites_bedtoolsFormat.txt
echo "Made all intersection files."

python ${project_rootdir}/tools/makeOverlapMatrix_significant.py \
    --input_dirname rosmap_significantCpgSites\
    --cpg_filename rosmap_significantCpgSites_bedtoolsFormat.txt\
    --original_eqtm_filepath ${eqtmFile} \
    --name rosmap_significantCpgSites
echo "Made overlapRatio."
