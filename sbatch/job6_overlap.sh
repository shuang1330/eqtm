#!/bin/bash
#SBATCH --ntasks=2
#SBATCH --time=05:00:00
#SBATCH --job-name=gt0.5Overlap
#SBATCH --output=gt0.5Overlap.out

module load TensorFlow/1.5.0-foss-2015b-Python-3.6.3
module load BEDTools/2.25.0-foss-2015b
cd /groups/umcg-gcc/tmp03/umcg-sli/eqtm_project
# python ./boxy_tools/buildOverlapRatioTable_boxy.py --rootDir /groups/umcg-gcc/tmp03/umcg-sli/eqtm_project --featureRootDir /groups/umcg-bios/tmp03/projects/2018-methylation/input --cpg random20k_gt0.5 --bedtools findOverlap
python finalOverlapMatrix.py
