#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --job-name=callBEDTools
#SBATCH --output=callBEDTools_createMatrixRatio.out

module purge
module load Python
module load BEDTools
project_rootdir="/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
cd ${project_rootdir}
python ./tools/create_bedtoolsIntersect_files.py \
--project_rootdir ${project_rootdir} \
--data_rootdir /groups/umcg-bios/tmp03/projects/2018-methylation/input/features/Roadmap \
--feature_filetype consolidatedImputedGappedPeak \
--input_type cpg_lung \
--cpg_filename lung_cpg_bedtoolsFormat.txt
#python ./tools/create_bedtoolsIntersect_files.py \
#--project_rootdir ${project_rootdir} \
#--data_rootdir /groups/umcg-bios/tmp03/projects/2018-methylation/input/features/Roadmap \
#--feature_filetype consolidatedImputedGappedPeak \
#--input_type cpg_fibroblast \
#--cpg_filename eqtm_fibroblast_bedtoolsFormat.txt

# original dataset
#python ./tools/create_bedtoolsIntersect_files.py \
#--project_rootdir ${project_rootdir} \
#--data_rootdir /groups/umcg-bios/tmp03/projects/2018-methylation/input/features/Roadmap \
#--feature_filetype consolidatedImputedGappedPeak \
#--input_type cpg_largerthan0.05 \
#--cpg_filename randomly60k_eqtm_FDR_larger_than_0.05_bedtoolsFormat.txt

# TCGA dataset
#python ./tools/create_bedtoolsIntersect_files.py \
#--project_rootdir ${project_rootdir} \
#--data_rootdir /groups/umcg-bios/tmp03/projects/2018-methylation/input/features/Roadmap \
#--feature_filetype consolidatedImputedGappedPeak \
#--input_type cpg_TCGA \
#--cpg_filename TCGA_all_cpgs_bedtoolsFormat.txt
echo "Finished bedtools operation."
#python tools/makeOverlapMatrix.py
