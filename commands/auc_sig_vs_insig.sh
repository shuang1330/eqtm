#!/usr/bin/env bash
#SBATCH --time=5:00:00
#SBATCH --job-name=sig_vs_Insig
#SBATCH --output=auc_sig_vs_insig.out
#SBATCH --mem=50gb

module purge
module load Python
project_dirpath="/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
cd ${project_dirpath}
python tools/assumption1_tensorboard.py