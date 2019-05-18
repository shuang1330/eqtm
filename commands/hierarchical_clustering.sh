#!/usr/bin/env bash
#SBATCH --time=1:00:00
#SBATCH --job-name=hierarchical_clustering
#SBATCH --output=hierarchical_clustering.out
#SBATCH --mem=50gb

module purge
module load TensorFlow/1.5.0-foss-2015b-Python-3.6.3
project_dirpath="/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
cd ${project_dirpath}
python tools/assumption1_tensorboard.py
