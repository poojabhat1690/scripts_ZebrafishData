#!/bin/bash
#$ -pe smp 6-24

ml star
mkdir /clustertmp/pooja/STARmapping/Bud/

##STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Bud/Bud_1_1.fq.gz /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Bud/Bud_1_2.fq.gz  --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/STARmapping/Bud/Bud_1 --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1

##3STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Bud/Bud_2_1.fq.gz /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Bud/Bud_2_2.fq.gz  --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/STARmapping/Bud/Bud_2 --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1

ml samtools 

samtools index /clustertmp/pooja/STARmapping/Bud/Bud_1Aligned.sortedByCoord.out.bam /clustertmp/pooja/STARmapping/Bud/Bud_1Aligned.sortedByCoord.out.bam.bai
