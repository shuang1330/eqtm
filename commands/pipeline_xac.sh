#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --job-name=xac
#SBATCH --output=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/out/cpg_xac.out
#SBATCH --error=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/out/cpg_xac.err
#SBATCH --mem=50gb

module purge
module load Python
module load BEDTools
project_rootdir="/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
#cd ${project_rootdir}

eqtmFile=${project_rootdir}/data/eqtmZscores/allCpgs/cpg_xac.txt
eqtmBedtoolsFile=${project_rootdir}/data/eqtmZscores/allCpgs/cpg_xac_bedtoolsFormat.txt
#awk 'BEGIN {OFS="\t"};NR>2 {print "chr"$2, $3-25, $3+25, $1}' ${eqtmFile} > ${eqtmBedtoolsFile}

#python ${project_rootdir}/tools/create_bedtoolsIntersect_files.py \
#    --project_rootdir ${project_rootdir} \
#    --data_rootdir /groups/umcg-bios/tmp03/projects/2018-methylation/input/features/Roadmap \
#    --feature_filetype consolidatedImputedGappedPeak \
#    --input_type cpg_xac \
#    --cpg_filename cpg_xac_bedtoolsFormat.txt
#echo "Made all intersection files."

python ${project_rootdir}/tools/makeOverlapMatrix_significant.py \
    --input_dirname cpg_xac\
    --cpg_filename cpg_xac_bedtoolsFormat.txt\
    --original_eqtm_filepath ${project_rootdir}/data/eqtmZscores/allCpgs/cpg_xac.txt \
    --name cpg_xac
echo "Made overlapRatio."
