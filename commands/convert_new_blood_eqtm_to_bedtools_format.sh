#!/usr/bin/env bash

eqtm_file=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/new_blood_eqtm/eqtm_fibroblast.txt
eqtm_bedtools=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/new_blood_eqtm/eqtm_fibroblast_bedtoolsFormat.txt
awk 'BEGIN {OFS="\t"};NR>2 {print "chr"$3, $5-25, $5+25, $1}' ${eqtm_file} > ${eqtm_bedtools}

eqtm_gene_bedtools=/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/new_blood_eqtm/gene_fibroblast_bedtoolsFormat.txt
awk 'BEGIN {OFS="\t"};NR>2 {print "chr"$4, $6-100, $6, $2}' ${eqtm_file} > ${eqtm_gene_bedtools}
