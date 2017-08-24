
#!/bin/bash
#$ -pe smp 6-24

mkdir /clustertmp/pooja/STARmapping/1KCell/
ml star

##STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/1KCell/1KCell_1_1.fq.gz /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/1KCell/1KCell_1_2.fq.gz  --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/STARmapping/1KCell/1KCell_1 --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1

##STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/1KCell/1KCell_2_1.fq.gz /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/1KCell/1KCell_2_2.fq.gz  --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/STARmapping/1KCell//KCell_2 --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 

ml samtools 

samtools index /clustertmp/pooja/STARmapping/1KCell//KCell_2Aligned.sortedByCoord.out.bam /clustertmp/pooja/STARmapping/1KCell//KCell_2Aligned.sortedByCoord.out.bam.bai
samtools index /clustertmp/pooja/STARmapping/1KCell/1KCell_1Aligned.sortedByCoord.out.bam /clustertmp/pooja/STARmapping/1KCell/1KCell_1Aligned.sortedByCoord.out.bam.bai




