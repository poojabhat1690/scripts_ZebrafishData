#!/bin/bash
#$ -pe smp 6-24

ml star
mkdir /clustertmp/pooja/STARmapping//Testis2/

#STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Testis2/zf_testis_2_5_F.fq.gz /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Testis2/zf_testis_2_5_R.fq.gz  --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/STARmapping/Testis2//zf_testis_2_5 --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1

#STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Testis2/zf_testis_2_6_F.fq.gz /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Testis2/zf_testis_2_6_R.fq.gz  --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/STARmapping/Testis2//zf_testis_2_6 --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1

#STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Testis2/zf_testis_2_7_F.fq.gz /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Testis2/zf_testis_2_7_R.fq.gz  --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/STARmapping/Testis2//zf_testis_2_7 --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1

ml samtools
samtools index /clustertmp/pooja/STARmapping/Testis2//zf_testis_2_5Aligned.sortedByCoord.out.bam /clustertmp/pooja/STARmapping/Testis2//zf_testis_2_5Aligned.sortedByCoord.out.bam.bai  

samtools index /clustertmp/pooja/STARmapping/Testis2//zf_testis_2_6Aligned.sortedByCoord.out.bam /clustertmp/pooja/STARmapping/Testis2//zf_testis_2_6Aligned.sortedByCoord.out.bam.bai

samtools index /clustertmp/pooja/STARmapping/Testis2//zf_testis_2_7Aligned.sortedByCoord.out.bam /clustertmp/pooja/STARmapping/Testis2//zf_testis_2_7Aligned.sortedByCoord.out.bam.bai
