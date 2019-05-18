#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --job-name=add_ratio
#SBATCH --output=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/out/make_predictions.out
#SBATCH --error=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/out/make_predictions.err
#SBATCH --mem=50gb

module purge
module load Python


project_rootdir="/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
echo "Start processing."
#python ${project_rootdir}/tools/add_overlapRatio_for_splits.py

python ${project_rootdir}/tools/make_predictions_for_public_dataset.py
