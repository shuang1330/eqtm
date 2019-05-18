#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --job-name=makeOverlapMatrix
#SBATCH --output=makeOverlapMatrix.out

module purge
module load Python
project_rootdir=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm
python ${project_rootdir}/tools/makeOverlapMatrix.py \
    --input_dirname cpg_lung\
    --cpg_filename lung_cpg_bedtoolsFormat.txt\
    --original_eqtm_filepath ${project_rootdir}/data/eqtmZscores/allCpgs/lung_cpg_bedtoolsFormat_with_header.txt \
    --name cpg_lung
# new fibroblast dataset
#python ${project_rootdir}/tools/makeOverlapMatrix.py \
#    --input_dirname cpg_fibroblast\
#    --cpg_filename eqtm_fibroblast_bedtoolsFormat.txt\
#    --original_eqtm_filepath ${project_rootdir}/data/eqtmZscores/allCpgs/eqtm_fibroblast_bedtoolsFormat_with_header.txt \
#    --name cpg_fibroblast
# TCGA dataset
#python ${project_rootdir}/tools/makeOverlapMatrix.py \
#    --input_dirname cpg_TCGA\
#    --cpg_filename TCGA_all_cpgs_bedtoolsFormat.txt\
#    --original_eqtm_filepath ${project_rootdir}/data/eqtmZscores/allCpgs/all_TCGA_cpgs.txt \
#    --name cpg_TCGA
echo "Finished."
