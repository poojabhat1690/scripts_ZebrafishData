#!/bin/bash


#### converting multiple files in Andi's dataset (only quantSeq) to fastq

ml bedtools

cd /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/2016_12_QuantSeq_Annotaton/raw_data

for f in ./*/*.bam
do

bamToFastq -i "$f" -fq "$f".fq 

gzip "$f".fq

done
