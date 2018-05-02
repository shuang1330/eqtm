#!/bin/bash
#SBATCH --ntasks=2
#SBATCH --time=00:30:00
#SBATCH --job-name=testjob

module load PythonPlus/2.7.10-foss-2015b-v16.11.1
cd /groups/umcg-wijmenga/tmp03/projects/eQTMPrediction/code
python overlapMatrixTest.py -output insignificantTest.csv --cpg /groups/umcg-gcc/tmp03/umcg-sli/eqtms/random20000-eQTLsFDR-gt0.05-flipped.txt -shift 1
python multOverlapsTest.py -ratio_output insignificantRes.csv -bin_output binOutput.csv -cpgs /groups/umcg-gcc/tmp03/umcg-sli/eqtms/random20000-eQTLsFDR-gt0.05-flipped.txt
