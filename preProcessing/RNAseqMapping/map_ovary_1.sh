
#!/bin/bash
#$ -pe smp 6-24

ml star
mkdir /clustertmp/pooja/STARmapping//Ovary1/

#STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Ovary1/zf_ovary_1_5_F.fq.gz /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Ovary1/zf_ovary_1_5_R.fq.gz  --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/STARmapping/Ovary1//zf_ovary_1_5 --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1
#STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Ovary1/zf_ovary_1_6_F.fq.gz /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Ovary1/zf_ovary_1_6_R.fq.gz  --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/STARmapping/Ovary1//zf_ovary_1_6 --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1
#STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Ovary1/zf_ovary_1_7_F.fq.gz /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Ovary1/zf_ovary_1_7_R.fq.gz  --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/STARmapping/Ovary1//zf_ovary_1_7 --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1



ml samtools 

samtools index /clustertmp/pooja/STARmapping/Ovary1//zf_ovary_1_7Aligned.sortedByCoord.out.bam  /clustertmp/pooja/STARmapping/Ovary1//zf_ovary_1_7Aligned.sortedByCoord.out.bam.bai

samtools index /clustertmp/pooja/STARmapping/Ovary1//zf_ovary_1_6Aligned.sortedByCoord.out.bam /clustertmp/pooja/STARmapping/Ovary1//zf_ovary_1_6Aligned.sortedByCoord.out.bam.bai

samtools index  /clustertmp/pooja/STARmapping/Ovary1//zf_ovary_1_5Aligned.sortedByCoord.out.bam  /clustertmp/pooja/STARmapping/Ovary1//zf_ovary_1_5Aligned.sortedByCoord.out.bam.bai
