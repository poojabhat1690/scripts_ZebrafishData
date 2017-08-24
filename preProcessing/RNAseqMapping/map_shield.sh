#!/bin/bash
#$ -pe smp 6-24

ml star
mkdir /clustertmp/pooja/STARmapping//Shield/

#STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Shield/Shield_1_1.fq.gz /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Shield/Shield_1_2.fq.gz  --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/STARmapping/Shield/Shield_1 --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1

##STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Shield/Shield_2_1.fq.gz /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Shield/Shield_2_2.fq.gz  --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/STARmapping/Shield/Shield_2 --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1

##STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Shield/Shield_3_1.fq.gz /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Shield/Shield_3_2.fq.gz  --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/STARmapping/Shield/Shield_3 --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1


ml samtools

samtools index /clustertmp/pooja/STARmapping/Shield/Shield_1Aligned.sortedByCoord.out.bam /clustertmp/pooja/STARmapping/Shield/Shield_1Aligned.sortedByCoord.out.bam.bai

samtools index /clustertmp/pooja/STARmapping/Shield/Shield_2Aligned.sortedByCoord.out.bam /clustertmp/pooja/STARmapping/Shield/Shield_2Aligned.sortedByCoord.out.bam.bai

samtools index /clustertmp/pooja/STARmapping/Shield/Shield_3Aligned.sortedByCoord.out.bam  /clustertmp/pooja/STARmapping/Shield/Shield_3Aligned.sortedByCoord.out.bam.bai
