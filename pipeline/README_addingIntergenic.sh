#!/bin/bash
# -pe smp 6-24

CELLLINE="dr"
PIPELINE="/groups/ameres/Pooja/Projects/zebrafishAnnotation/pipeline/"
ANNOBASE="/groups/ameres/Pooja/Projects/zebrafishAnnotation/"
GBASE="/groups/ameres/bioinformatics/references/danio_rerio/dr10//"
BIN="//groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10_data/"
GBUILD="dr10"
GET_ANNOTATION="FALSE"
BOUT="/clustertmp/bioinfo/pooja/SLAMannotation/$CELLLINE/output"
 ENSEMBL="$ANNOBASE/$GBUILD/ensembl_dr10_Ensembl_Genes_88/"
         UCSC="$ANNOBASE/$GBUILD/refSeq_dr10_GRCh38_20170504/"
         GENOME="$GBASE/danRer10.fa"

ml R/3.3.0
ml bedtools
Rscript --vanilla -e "BIn='$BIN'; BOut='$BOUT'; ucscDir='$UCSC'; ensemblDir='$ENSEMBL'; source('$PIPELINE/intergenicPeaks/addingIntergenicPeaks.R')"
