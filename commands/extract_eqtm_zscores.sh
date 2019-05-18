#!/usr/bin/env bash


eqtmFile=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/ORIGIN/full_dataset.txt
eqtmZscoreFile=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/ORIGIN/full_dataset_zscores.txt
awk 'BEGIN {OFS="\t"}; {print $2, $5, $11}' ${eqtmfile} > ${eqtmZscoreFile}