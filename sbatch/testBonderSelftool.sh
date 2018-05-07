#!/bin/bash
#SBATCH --ntasks=2
#SBATCH --time=00:30:00
#SBATCH --job-name=testjob

module load TensorFlow/1.5.0-foss-2015b-Python-3.6.3
module load BEDTools/2.25.0-foss-2015b
python ./boxy_tools/buildOverlapRatioTable_boxy.py --rootDir /groups/umcg-gcc/tmp03/umcg-sli/eqtm_project --featureRootDir /groups/umcg-wijmenga/tmp03/projects/eQTMPrediction --cpg bonder-eQTMsFDR0.0-CpGLevel-split --bedtools findOverlap
