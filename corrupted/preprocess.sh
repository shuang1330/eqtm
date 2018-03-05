#!/bin/bash

export PROJECT_ROOTDIR="/home/shuang/projects/eqtl/corrupted"

# tmp folder for all unzipped narrowPeak files
if [ -d $PROJECT_ROOTDIR/tmp ]; then
  echo "deleted earlier files in tmp folder"
  rm -f tmp/*
fi

if [ ! -d $PROJECT_ROOTDIR/tmp ]; then
  echo "created tmp folder"
  mkdir tmp
fi

# unzip the narrowPeak files
for filename in $PROJECT_ROOTDIR/ori_files/dna_files/*.narrowPeak.gz;
do
  newname=$(basename "$filename" .narrowPeak.gz)
  yes n | gzip -dck $filename >> $PROJECT_ROOTDIR/tmp/$newname.narrowPeak
  awk 'BEGIN {OFS="\t"};{print $1,$2,$3,FILENAME}' $PROJECT_ROOTDIR/tmp/$newname.narrowPeak > $PROJECT_ROOTDIR/tmp/$newname.with_filename.narrowPeak
done
# concatenate all narrowPeak files$
cat $PROJECT_ROOTDIR/tmp/*.with_filename.narrowPeak > concated_dna.txt

# # delete everything in the tmp folder
# rm tmp/*

# # deleted the chr for narrowPeak file and sort it
# awk 'BEGIN {OFS="\t"};{$1=substr($1,4); print $0}' concated_dna.txt > dna.txt
# sort -k1 -k2 -V -s dna.txt > sorted_dna.txt
# awk 'BEGIN {OFS="\t"};{$1="chr"$1; print $0}' sorted_dna.txt > dna_final.txt
# rm -r concated_dna.txt dna.txt sorted_dna.txt
#
# # add start and end position to the wg files
# awk 'NR>1' $PROJECT_ROOTDIR/ori_files/eQTMsUniqueCGs-FDR0.05.txt > cgs_no_header.txt
# awk 'BEGIN {OFS="\t"};{print $2,$3-=25,$3+=25,$NF}' cgs_no_header.txt > cgs_tmp.txt
# sort -k1 -k2 -n cgs_tmp.txt > cgs_sorted.txt
# awk 'BEGIN {OFS="\t"};{$1="chr"$1; print $0}' cgs_sorted.txt > cgs_final.txt
# rm -r cgs_no_header.txt cgs_tmp.txt cgs_sorted.txt
#
# # bedtools
# bedtools intersect -sorted -a dna_final.txt -b cgs_final.txt -wa -wb > overlap.txt
