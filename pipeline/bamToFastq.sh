#!/bin/bash


cd /groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10_data/quantseq/
ml bedtools
files=$( ls *.bam )
counter=0
  for j in $files ; do
        echo "$j"
        let counter=$counter+1
        bedtools bamtofastq -i "$j" -fq "$j".fastq &>"$j".log

####        bam2fastq -o "$j".fq "$j" &>"$j".log    
done
