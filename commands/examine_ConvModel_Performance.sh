#!/usr/bin/env bash
#SBATCH --cpus-per-task=8
#SBATCH --time=2:00:00
#SBATCH --job-name=convModel
#SBATCH --output=Ranfor.out
#SBATCH --mem=50gb

module purge
module load TensorFlow/1.5.0-foss-2015b-Python-3.6.3
echo "Started to build ranfor model."
python tools/withDir_examine_convModel_performance.py
