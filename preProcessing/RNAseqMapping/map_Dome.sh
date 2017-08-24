
#!/bin/bash
#$ -pe smp 6-24

ml star
mkdir /clustertmp/pooja/STARmapping//Dome/

##STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Dome/Dome_1_1.fq.gz /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Dome/Dome_1_2.fq.gz  --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/STARmapping/Dome/Dome_1 --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1

##STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Dome/Dome_2_1.fq.gz /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Dome/Dome_2_2.fq.gz  --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/STARmapping/Dome/Dome_2 --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1

ml samtools
samtools index /clustertmp/pooja/STARmapping/Dome/Dome_1Aligned.sortedByCoord.out.bam /clustertmp/pooja/STARmapping/Dome/Dome_1Aligned.sortedByCoord.out.bam.bai
