#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --job-name=bios_histoneAnnotation
#SBATCH --output=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/out/bios_histoneAnnotation.out
#SBATCH --error=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/out/bios_histoneAnnotation.err
#SBATCH --mem=50gb

module purge
module load Python
module load BEDTools
project_rootdir="/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"

eqtmFile=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/ORIGIN/et_gt.txt
eqtmBedtoolsFile=${project_rootdir}/data/eqtmZscores/allCpgs/bios_bedtoolsFormat.txt
#python /groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/tools/get_CpgSite_Position.py \
#--input_path ${eqtmFile} \
#--bedtoolFormat_savepath ${eqtmBedtoolsFile}
##awk 'BEGIN {OFS="\t"};NR>2 {print "chr"$3, $4-25, $4+25, $2}' ${eqtmFile} > ${eqtmBedtoolsFile}
#head ${eqtmBedtoolsFile}

python ${project_rootdir}/tools/create_bedtoolsIntersect_files.py \
    --project_rootdir ${project_rootdir} \
    --data_rootdir /groups/umcg-bios/tmp03/projects/2018-methylation/input/features/Roadmap \
    --feature_filetype consolidatedImputedGappedPeak \
    --input_type bios \
    --cpg_filename bios_bedtoolsFormat.txt
echo "Made all intersection files."

python ${project_rootdir}/tools/makeOverlapMatrix_significant.py \
    --input_dirname bios \
    --cpg_filename bios_bedtoolsFormat.txt\
    --original_eqtm_filepath ${eqtmFile} \
    --name bios
echo "Made overlapRatio."