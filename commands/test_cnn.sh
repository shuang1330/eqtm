#!/usr/bin/env bash
#SBATCH --time=10:00:00
#SBATCH --job-name=test_cnn
#SBATCH --output=training_batchsize3e10.out
#SBATCH --mem=200gb

module purge
module load TensorFlow/1.5.0-foss-2015b-Python-3.6.3
project_dirpath="/groups/umcg-gcc/tmp04/umcg-sli/development_eqtm"
cd ${project_dirpath}
#python ./tools/complete_with_emptyFeatureFiles.py
#python ./tools/transform_into2d_data.py
python ./tools/testCNN.py \
--project_rootdir /groups/umcg-gcc/tmp04/umcg-sli/development_eqtm \
--eqtm_filename 2017-12-09-eQTLsFDR-gt0_withExpressionTssMethyOverlap_withGeneOverlap \
--train_iter 20000 \
--save_step 10000 \
--batch_size 200 \
--lr 0.0001
echo "Finished training."
sbatch ${project_dirpath}/commands/cnn_ranfor_2.sh
echo "Submitted job for testing iter 20000."
