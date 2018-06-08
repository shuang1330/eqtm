#!/bin/bash
#SBATCH --ntasks=2
#SBATCH --time=04:00:00
#SBATCH --job-name=cpgOverlap
#SBATCH --output=cpgOverlap.out

module purge
module load TensorFlow/1.5.0-foss-2015b-Python-3.6.3
module load BEDTools/2.25.0-foss-2015b
cd /groups/umcg-gcc/tmp03/umcg-sli/development_eqtm
python ./tools/build_overlapMatrix.py \
      --rootDir /groups/umcg-gcc/tmp03/umcg-sli/development_eqtm \
      --featureRootDir /groups/umcg-bios/tmp03/projects/2018-methylation/input \
      --bedtools findOverlap
