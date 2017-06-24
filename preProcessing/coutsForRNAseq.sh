#!/bin/bash
# -pe smp 6-24


##/ this counts counts from RNAseq data
ml subread


featureCounts -B -s 1 -a /groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10/GTFfiles/Danio_rerio.GRCz10.89.gtf -o /groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/rnaSeqReads_featureCounts.txt /clustertmp/pooja/STARmapping/allMapped/*.bam
