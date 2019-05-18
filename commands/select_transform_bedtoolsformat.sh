#!/usr/bin/env bash
#SBATCH --time=1:00:00
#SBATCH --job-name=select_largerthan_0.5
#SBATCH --output=select_eqtm_largerthan_0.5.out
#SBATCH --mem=50gb

module purge
module load TensorFlow/1.5.0-foss-2015b-Python-3.6.3
project_dirpath="/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
cd ${project_dirpath}
python tools/select_transform_bedtoolsFormat.py \
--project_rootdir ${project_dirpath} \
--eqtm_filename full_dataset.txt \
--cond_colName FDR \
--cond_operator larger_than \
--cond_threshold 0.5
echo "Finished transforming the file into bedtools format, started to call bedtools."
#sbatch commands/call_bedtools_build_matrix_ratio.sh
