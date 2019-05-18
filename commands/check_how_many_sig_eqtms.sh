#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --job-name=prepare_bios_split_datasets
#SBATCH --output=prepare_bios_split_datasets.out

module purge
module load Python
module load BEDTools
project_rootdir="/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
cd ${project_rootdir}

# concatenate every chromosome together
# bios data part1
bioseQTMsPart1=/groups/umcg-gcc/tmp03/umcg-sli/bios/cis_eqtm_part1
alleQTMsFile_part1=${project_rootdir}/data/eqtmZscores/ORIGIN/BIOS_Part1.txt
echo "Part1"
for i in {2..22}
do
    zcat ${bioseQTMsPart1}/${i}/eQTLsFDR0.05.txt.gz | wc
done
# bios data part2
bioseQTMsPart2=/groups/umcg-gcc/tmp03/umcg-sli/bios/cis_eqtm_part2
alleQTMsFile_part2=${project_rootdir}/data/eqtmZscores/ORIGIN/BIOS_Part2.txt
echo "Part2"
for i in {2..22}
do
    zcat ${bioseQTMsPart2}/${i}/eQTLsFDR0.05.txt.gz | wc
done

#Part1
#    640   16638  201167
#    369    9592  115646
#    181    4704   56615
#    702   18250  220768
#   7650  198898 2414147
#    665   17288  208206
#    272    7070   85467
#    184    4782   57675
#    415   10788  130477
#    566   14714  177610
#    633   16456  199586
#     85    2208   26702
#    312    8110   98537
#    255    6628   80241
#    437   11360  136935
#    723   18796  228263
#     99    2572   30850
#   3817   99240 1210123
#    142    3690   44318
#    248    6446   77885
#    256    6654   80833
#Part2
#   1155   30028  351193
#    792   20590  240672
#    345    8968  103886
#   1542   40090  468635
#   8438  219386 2565362
#   1391   36164  420822
#    599   15572  181585
#    493   12816  149697
#    807   20980  246352
#   1509   39232  460128
#   1314   34162  399811
#    184    4782   55639
#    964   25062  296641
#    427   11100  128945
#   1095   28468  334144
#   1989   51712  606740
#    156    4054   46627
#   7025  182648 2133833
#    316    8214   95884
#    401   10424  121709
#    565   14688  172039