#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --job-name=corr
#SBATCH --output=sig_insig_clustermap.out
#SBATCH --mem=100gb

export MPLBACKEND="agg"
module purge
module load Python/3.6.3-foss-2015b
#python tools/spearmanr_methy_Sites.py
python tools/assumption1_tensorboard.py