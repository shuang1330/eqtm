#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --job-name=makeOverlapMatrix_gene
#SBATCH --output=makeOverlapMatrix_gene.out

module purge
module load Python
project_rootdir=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm
# new fibroblast dataset
python ${project_rootdir}/tools/makeOverlapMatrix.py \
    --input_dirname cpg_fibroblast_gene\
    --cpg_filename gene_fibroblast_bedtoolsFormat.txt\
    --original_eqtm_filepath ${project_rootdir}/data/eqtmZscores/allCpgs/gene_fibroblast_bedtoolsFormat_with_header.txt \
    --name cpg_fibroblast_gene
# TCGA dataset
#python ${project_rootdir}/tools/makeOverlapMatrix.py \
#    --input_dirname cpg_TCGA\
#    --cpg_filename TCGA_all_cpgs_bedtoolsFormat.txt\
#    --original_eqtm_filepath ${project_rootdir}/data/eqtmZscores/allCpgs/all_TCGA_cpgs.txt \
#    --name cpg_TCGA
echo "Finished."
