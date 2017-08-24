#!/bin/bash

#$ -t 1-48

#$ -pe smp 6-24

########$ -wd /home/<your_UCL_id>/Scratch/output

# 8. Parse parameter file to get variables.
number=$SGE_TASK_ID
paramfile=/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/scripts_ZebrafishData/ucscBrowser/factors_normalize.txt

index=`sed -n ${number}p $paramfile | awk '{print $1}'`
variable1=`sed -n ${number}p $paramfile | awk '{print $2}'`
FAC=`sed -n ${number}p $paramfile | awk '{print $2}'`

USENAME=`sed -n ${number}p $paramfile | awk '{print $1}'`
####variable2=`sed -n ${number}p $paramfile | awk '{print $3}'`
##variable3=`sed -n ${number}p $paramfile | awk '{print $4}'`





ml deeptools/2.5.0.1-python2.7.3
ml pysam/0.10.0


cd  /clustertmp/pooja/STARmapping/allMapped/

mkdir /clustertmp/pooja/STARmapping/allMapped/splitStrand/
cd /clustertmp/pooja/STARmapping/allMapped/splitStrand/

mkdir /groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/ucscTracks/
OUTDIR=//groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/ucscTracks/

ml samtools


### getting the correct flags from the forward strand

samtools view -b -f 163 /clustertmp/pooja/STARmapping/allMapped/"$index" > "$index"_fwd1.bam
samtools view -b -f 83 /clustertmp/pooja/STARmapping/allMapped/"$index" > "$index"_fwd2.bam
samtools merge "$index"_fwd.bam "$index"_fwd1.bam "$index"_fwd2.bam

samtools index "$index"_fwd.bam "$index"_fwd.bam.bai

### getting the correct files from reverse strand


samtools view -b -f 99 /clustertmp/pooja/STARmapping/allMapped/"$index" > "$index"_rv1.bam
samtools view -b -f 147 /clustertmp/pooja/STARmapping/allMapped/"$index" > "$index"_rv2.bam
samtools merge "$index"_rv.bam "$index"_rv1.bam "$index"_rv2.bam

samtools index "$index"_rv.bam "$index"_rv.bam.bai

bamCoverage  -b "$index"_fwd.bam -o "$OUTDIR"/"$index"_plus.bg  --scaleFactor "$FAC" --binSize 1 -of bedgraph

bamCoverage  -b "$index"_rv.bam -o "$OUTDIR"/"$index"_minus.bg --scaleFactor "$FAC" --binSize 1 -of bedgraph


awk '{OFS="\t"; print $1,$2,$3,"-"$4}' "$OUTDIR"/"$index"_minus.bg >"$OUTDIR"/"$index"_neg_minus.bg


ml kent-ucsc/3.8f6f5e0a1cb75

bedGraphToBigWig "$OUTDIR"/"$index"_neg_minus.bg /groups/ameres/bioinformatics/references/danio_rerio/dr10/danRer10.chrom.sizes  "$OUTDIR"/"$index"_neg_minus.bg.bigWig


bedGraphToBigWig "$OUTDIR"/"$index"_plus.bg /groups/ameres/bioinformatics/references/danio_rerio/dr10/danRer10.chrom.sizes  "$OUTDIR"/"$index"_plus.bg.bigWig






