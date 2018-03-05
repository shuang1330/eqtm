#!/bin/bash

module load BEDTools/2.25.0-foss-2015b
echo "module loaded"

export PROJECT_ROOTDIR="/groups/umcg-wijmenga/tmp03/projects/eQTMPrediction"

# tmp folder for all unzipped narrowPeak files
if [ -d $PROJECT_ROOTDIR/tmp ]; then
  echo "deleted earlier files in tmp folder"
  rm -f tmp/*
fi

if [ ! -d $PROJECT_ROOTDIR/tmp ]; then
  echo "created tmp folder"
  mkdir tmp
fi

echo "tmp folder is empty now"

# unzip the narrowPeak files
for filename in $PROJECT_ROOTDIR/features/Roadmap/consolidatedNarrowPeak/*.H3K27ac.narrowPeak.gz;
do
  name=$(basename "$filename" .gz)
  yes n | gzip -dc $filename >> $PROJECT_ROOTDIR/tmp/$(basename "$filename" .gz)
  awk 'BEGIN {OFS="\t"};{print $1,$2,$3,FILENAME}' $PROJECT_ROOTDIR/tmp/$name.narrowPeak > $PROJECT_ROOTDIR/tmp/$name.with_filename.narrowPeak
  rm $PROJECT_ROOTDIR/tmp/$name.narrowPeak
done

# concatenate all narrowPeak files
cat $PROJECT_ROOTDIR/tmp/* > concated_dna.txt | head -5
echo "all dna files have been concatenated"

# delete everything in the tmp folder
rm tmp/*

# deleted the chr for narrowPeak file and sort it
awk 'BEGIN {OFS="\t"};{$1=substr($1,4); print $0}' concated_dna.txt > dna.txt
sort -k1 -k2 -V -s dna.txt > sorted_dna.txt
awk 'BEGIN {OFS="\t"};{$1="chr"$1; print $0}' sorted_dna.txt > $PROJECT_ROOTDIR/output/dna_final.txt
rm -r concated_dna.txt dna.txt sorted_dna.txt

cat $PROJECT_ROOTDIR/output/dna_final.txt | head -5
echo "sorted the dna files"

# add start and end position to the wg files
awk 'NR>1' $PROJECT_ROOTDIR/eqtms/2017-12-09-eQTMs-1mb/eQTMsUniqueCGs-FDR0.05.txt > cgs_no_header.txt
awk 'BEGIN {OFS="\t"};{print $2,$3-=25,$3+=25,$1,$NF}' cgs_no_header.txt > cgs_tmp.txt
sort -k1 -k2 -n cgs_tmp.txt > cgs_sorted.txt
awk 'BEGIN {OFS="\t"};{$1="chr"$1; print $0}' cgs_sorted.txt > $PROJECT_ROOTDIR/output/cgs_final.txt
rm -r cgs_no_header.txt cgs_tmp.txt cgs_sorted.txt

cat $PROJECT_ROOTDIR/output/cgs_final.txt | head -5
echo "sorted cgs files"

# bedtools
bedtools intersect -sorted -a $PROJECT_ROOTDIR/output/dna_final.txt -b $PROJECT_ROOTDIR/output/cgs_final.txt -wa -wb > $PROJECT_ROOTDIR/output/bedtools_intersect_output.txt
cat $PROJECT_ROOTDIR/output/bedtools_intersect_output.txt | head -5
echo "get the overlapping result"
