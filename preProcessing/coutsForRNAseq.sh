
#!/bin/bash
# -pe smp 6-24


##/ this counts counts from RNAseq data
ml subread


featureCounts -T 15 -B  -p -g gene_id -a  /groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10/GTFfiles/Danio_rerio.GRCz10.89.gtf -o /groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/rnaSeqReads_featureCounts_strandSpecific_allGenes.txt /clustertmp/pooja/STARmapping/allMapped/*.bam
