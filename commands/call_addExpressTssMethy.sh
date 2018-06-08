#!/bin/bash
#SBATCH --ntasks=2
#SBATCH --time=02:00:00
#SBATCH --job-name=addINFOModel
#SBATCH --output=addExpressTssMethyOverlap_modelRanfor.out

module load TensorFlow/1.5.0-foss-2015b-Python-3.6.3
module load BEDTools/2.25.0-foss-2015b
cd /groups/umcg-gcc/tmp03/umcg-sli/development_eqtm
python tools/equip_eqtm_withExpression_Tss_Methy.py --dataDir /groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data
python tools/finetune_conv_models.py
