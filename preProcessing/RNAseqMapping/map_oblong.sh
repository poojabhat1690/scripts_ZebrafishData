

#!/bin/bash
#$ -pe smp 6-24

ml star

mkdir /clustertmp/pooja/STARmapping//Oblong/

##STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Oblong/Oblong_1.fq.gz /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Oblong/Oblong_2.fq.gz  --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/STARmapping/Oblong/Oblong_1 --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --limitBAMsortRAM 65324197591

ml samtools

samtools index /clustertmp/pooja/STARmapping/Oblong/Oblong_1Aligned.sortedByCoord.out.bam /clustertmp/pooja/STARmapping/Oblong/Aligned.sortedByCoord.out.bam.bai



