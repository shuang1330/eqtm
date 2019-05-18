#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --job-name=prepare_bios_split_datasets
#SBATCH --output=prepare_bios_split_datasets.out

module purge
module load Python
module load BEDTools
project_rootdir="/groups/umcg-gcc/tmp03/umcg-sli/replication_output"
cd ${project_rootdir}

# concatenate every chromosome together
bioseQTMs=${project_rootdir}/eqtm_found_in_cordBlood_replicate_in_bios/
alleQTMsFile=${project_rootdir}/concat_replicate_results/cordBlood_replicates_in_bios.txt
zcat ${bioseQTMs}/1/eQTLs.txt.gz > ${alleQTMsFile}
for i in {2..22}
do
    echo "chr"${i}
    zcat ${bioseQTMs}/${i}/eQTLs.txt.gz | tail -n +2 >> ${alleQTMsFile}
done

livereQTMs=${project_rootdir}/eqtm_found_in_liver_metaTCGA_replicate_in_bios/
allLivereQTMsFile=${project_rootdir}/concat_replicate_results/liver_metaTCGA_replicates_in_bios.txt
zcat ${livereQTMs}/1/eQTLs.txt.gz > ${allLivereQTMsFile}
for i in {2..22}
do
    echo "chr"${i}
    zcat ${livereQTMs}/${i}/eQTLs.txt.gz | tail -n +2 >> ${allLivereQTMsFile}
done
