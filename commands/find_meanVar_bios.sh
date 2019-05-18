#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --job-name=find_meanVar
#SBATCH --output=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/out/meanVar_bios.out
#SBATCH --error=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/out/meanVar_bios.err
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --mem=50gb

module purge
module load Python

project_rootdir="/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"

python ${project_rootdir}/tools/find_meanVar_bios.py

