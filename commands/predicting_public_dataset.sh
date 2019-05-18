#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --job-name=make_predictions_for_public_dataset
#SBATCH --output=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/out/predict_public.out
#SBATCH --error=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/out/predict_public.err
#SBATCH --mem=20gb

module purge
module load Python
# make predictions
python /groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/tools/make_predictions_for_public_dataset.py