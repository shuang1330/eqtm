#!/usr/bin/env bash
#SBATCH --time=2:00:00
#SBATCH --job-name=conv1d_test
#SBATCH --output=conv1d_test.out
#SBATCH --mem=50gb

module purge
module load TensorFlow/1.5.0-foss-2015b-Python-3.6.3
project_dirpath="/groups/umcg-gcc/tmp04/umcg-sli/development_eqtm"
cd ${project_dirpath}
python tools/test_1dconv.py