#!/bin/bash
#SBATCH --ntasks=2
#SBATCH --time=9:00:00
#SBATCH --job-name=testjob

module load TensorFlow/1.5.0-foss-2015b-Python-3.6.3
module load BEDTools/2.25.0-foss-2015b
cd /groups/umcg-gcc/tmp03/umcg-sli/eqtm_project
# old data path
# python ./boxy_tools/buildOverlapRatioTable_boxy.py --rootDir /groups/umcg-gcc/tmp03/umcg-sli/eqtm_project --featureRootDir /groups/umcg-wijmenga/tmp03/projects/eQTMPrediction --cpg random20000-eQTLsFDR-gt0.05-flipped --bedtools findOverlap
# new data path in bios folder
python ./boxy_tools/buildOverlapRatioTable_boxy.py --rootDir /groups/umcg-gcc/tmp03/umcg-sli/eqtm_project --featureRootDir /groups/umcg-bios/tmp03/projects/2018-methylation/input --cpg random20k_gt0.5 --bedtools findOverlap
python ./boxy_tools/buildOverlapRatioTable_boxy.py --rootDir /groups/umcg-gcc/tmp03/umcg-sli/eqtm_project --featureRootDir /groups/umcg-bios/tmp03/projects/2018-methylation/input --cpg random20k_gt0.4_sm0.5 --bedtools findOverlap
python ./boxy_tools/buildOverlapRatioTable_boxy.py --rootDir /groups/umcg-gcc/tmp03/umcg-sli/eqtm_project --featureRootDir /groups/umcg-bios/tmp03/projects/2018-methylation/input --cpg random20k_gt0.3_sm0.4 --bedtools findOverlap
python ./boxy_tools/buildOverlapRatioTable_boxy.py --rootDir /groups/umcg-gcc/tmp03/umcg-sli/eqtm_project --featureRootDir /groups/umcg-bios/tmp03/projects/2018-methylation/input --cpg random20k_gt0.2_sm0.3 --bedtools findOverlap
python ./boxy_tools/buildOverlapRatioTable_boxy.py --rootDir /groups/umcg-gcc/tmp03/umcg-sli/eqtm_project --featureRootDir /groups/umcg-bios/tmp03/projects/2018-methylation/input --cpg random20k_gt0.1_sm0.2 --bedtools findOverlap
python ./boxy_tools/buildOverlapRatioTable_boxy.py --rootDir /groups/umcg-gcc/tmp03/umcg-sli/eqtm_project --featureRootDir /groups/umcg-bios/tmp03/projects/2018-methylation/input --cpg random20k_gt0.05_sm0.1 --bedtools findOverlap
# locally
# python ./tools/buildOverlapRatioTable.py --rootDir /home/shuang/projects/eqtm --cpg random20000-eQTLsFDR-gt0.05-flipped --bedtools findOverlap
