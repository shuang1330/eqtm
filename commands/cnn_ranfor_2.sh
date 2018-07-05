#!/usr/bin/env bash
#SBATCH --cpus-per-task=8
#SBATCH --time=8:00:00
#SBATCH --job-name=test_cnn
#SBATCH --output=Iter2_cnnRanfor.out
#SBATCH --mem=100gb

module purge
module load TensorFlow/1.5.0-foss-2015b-Python-3.6.3
echo "Started to transform embeddings to dataframe."
python ./tools/cnn_ranfor.py \
--ckpt_filename test_CNN_trainIter_20000_lr0.000100_batchsize100 \
--trainIter_batch_lr trainIter_20000_lr0.000100_batchsize100
echo "Started to examine ranfor performance."
python ./tools/examine_convModel_performance.py \
--embedding_filename 2017-12-09-eQTLsFDR-gt0_withExpressionTssMethy_trainIter_20000_lr0.000100_batchsize100_embeddings