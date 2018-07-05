#!/usr/bin/env bash
#SBATCH --time=2:00:00
#SBATCH --job-name=overlap
#SBATCH --output=overlap.out
#SBATCH --mem=50gb

python tools/intersect_no_name_version_makeOverlapMatrix.py
python tools/complete_with_emptyFeatureFiles.py
