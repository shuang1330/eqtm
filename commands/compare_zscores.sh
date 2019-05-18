#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --job-name=compare_zscores
#SBATCH --output=compare_zscores.out

module purge
module load Python
#eqtmFile=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/ORIGIN/full_dataset.txt
#eqtmZscoreFile=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/ORIGIN/full_dataset_zscores.txt
#awk 'BEGIN {OFS="\t"}; {print $2, $5, $11}' ${eqtmfile} > ${eqtmZscoreFile}
project_rootdir=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm
python ${project_rootdir}/tools/compare_zscores_different_datasets.py