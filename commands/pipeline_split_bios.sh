#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --job-name=prepare_bios_split_datasets
#SBATCH --output=prepare_bios_split_datasets.out

module purge
module load Python
module load BEDTools
project_rootdir="/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
cd ${project_rootdir}

# concatenate every chromosome together
# bios data part1
#bioseQTMsPart1=/groups/umcg-gcc/tmp03/umcg-sli/bios_resplit/cis_eqtm_part1
#alleQTMsFile_part1=${project_rootdir}/data/eqtmZscores/ORIGIN/BIOS_resplit_Part1.txt
#zcat ${bioseQTMsPart1}/1/eQTLsFDR0.05.txt.gz > ${alleQTMsFile_part1}
#for i in {2..22}
#do
#    echo "chr"${i}
#    zcat ${bioseQTMsPart1}/${i}/eQTLsFDR0.05.txt.gz | tail -n +2 >> ${alleQTMsFile_part1}
#done
## bios data part2
#bioseQTMsPart2=/groups/umcg-gcc/tmp03/umcg-sli/bios_resplit/cis_eqtm_part2
#alleQTMsFile_part2=${project_rootdir}/data/eqtmZscores/ORIGIN/BIOS_resplit_Part2.txt
#zcat ${bioseQTMsPart2}/1/eQTLsFDR0.05.txt.gz > ${alleQTMsFile_part2}
#for i in {2..22}
#do
#    echo "chr"${i}
#    zcat ${bioseQTMsPart2}/${i}/eQTLsFDR0.05.txt.gz | tail -n +2 >> ${alleQTMsFile_part2}
#done
#
## create bedtoolsFormat files
bedtoolsFormat_Part1_2=${project_rootdir}/data/eqtmZscores/allCpgs/bios_resplit_part1_2_bedtoolsFormat_with_duplicates.txt
#awk 'BEGIN {OFS="\t"};NR>2 {print "chr"$3, $4-25, $4+25, $2}' ${alleQTMsFile_part1} > ${bedtoolsFormat_Part1_2}
#awk 'BEGIN {OFS="\t"};NR>2 {print "chr"$3, $4-25, $4+25, $2}' ${alleQTMsFile_part2} >> ${bedtoolsFormat_Part1_2}
#echo "Produced bedtoolsFormat files."
#
## remove duplicates in bedtoolsFormat Files
bedtoolsFormat_Part1_2_no_duplicates=${project_rootdir}/data/eqtmZscores/allCpgs/bios_resplit_part1_2_bedtoolsFormat.txt
#awk '!x[$0]++' ${bedtoolsFormat_Part1_2} > ${bedtoolsFormat_Part1_2_no_duplicates}
#echo "Produced bedtoolsFormat file without duplicates."

# create_intersection files
python ./tools/create_bedtoolsIntersect_files.py \
    --project_rootdir ${project_rootdir} \
    --data_rootdir /groups/umcg-bios/tmp03/projects/2018-methylation/input/features/Roadmap \
    --feature_filetype consolidatedImputedGappedPeak \
    --input_type cpg_bios_resplit \
    --cpg_filename bios_resplit_part1_2_bedtoolsFormat.txt
echo "Intersection files created!"

# create header for bedtoolsformat files
bedtoolsFormat_Part1_2_no_duplicates_with_header=${project_rootdir}/data/eqtmZscores/allCpgs/bios_resplit_part1_2_bedtoolsFormat_with_header.txt
touch ${bedtoolsFormat_Part1_2_no_duplicates_with_header}
echo -e "chr\tstartSite\tendSite\tSNPName" > ${bedtoolsFormat_Part1_2_no_duplicates_with_header}
awk 'NR>2' ${bedtoolsFormat_Part1_2_no_duplicates} >> ${bedtoolsFormat_Part1_2_no_duplicates_with_header}

# make overlapRatio matrix
python ${project_rootdir}/tools/makeOverlapMatrix.py \
    --input_dirname cpg_bios_resplit\
    --cpg_filename bios_resplit_part1_2_bedtoolsFormat.txt\
    --original_eqtm_filepath ${project_rootdir}/data/eqtmZscores/allCpgs/bios_resplit_part1_2_bedtoolsFormat_with_header.txt \
    --name cpg_bios_resplit
echo "Made overlapRatio matrix."

# add all features to eqtm files
python ${project_rootdir}/tools/add_all_features_for_bios_data.py \
    --name BIOS_resplit_Part1 \
    --overlapRatio_filename cpg_bios_resplit_overlapRatio
python ${project_rootdir}/tools/add_all_features_for_bios_data.py \
    --name BIOS_resplit_Part2 \
    --overlapRatio_filename cpg_bios_resplit_overlapRatio

# flip the direction

# modeling and analysis
#python ${project_rootdir}/tools/withDir_examine_convModel_performance.py
