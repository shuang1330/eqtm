#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --job-name=add_all_features
#SBATCH --output=add_all_features.out

module purge
module load Python
project_rootdir=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm
for eqtm_name in /groups/umcg-gcc/tmp03/umcg-sli/tcga/cis_results/*
do
    name=$(basename ${eqtm_name})
    echo ${name}
    python ${project_rootdir}/tools/add_all_features_to_eqtm_file.py --name ${name}
done
#python ${project_rootdir}/tools/newData_add_all_features_to_eqtm_file.py