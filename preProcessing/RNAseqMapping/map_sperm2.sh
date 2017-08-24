#!/bin/bash
#$ -pe smp 6-24

ml star
mkdir /clustertmp/pooja/STARmapping//Sperm2/

#STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Sperm2/zf_sperm2_longRNA_2_F.fq.gz /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Sperm2/zf_sperm2_longRNA_2_R.fq.gz  --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/STARmapping/Sperm2//zf_sperm2_longRNA_2 --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1



ml samtools 

samtools index /clustertmp/pooja/STARmapping/Sperm2//zf_sperm2_longRNA_2Aligned.sortedByCoord.out.bam /clustertmp/pooja/STARmapping/Sperm2//zf_sperm2_longRNA_2Aligned.sortedByCoord.out.bam.bai
