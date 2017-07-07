#!/bin/bash
#$ -pe smp 6-24

ml star
mkdir /clustertmp/pooja/STARmapping//Ovary2/

STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Ovary2/zf_ovary_2_5_F.fq.gz /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Ovary2/zf_ovary_2_5_R.fq.gz  --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/STARmapping/Ovary2//zf_ovary_2_5 --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1
STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Ovary2/zf_ovary_2_6_F.fq.gz /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Ovary2/zf_ovary_2_6_R.fq.gz  --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/STARmapping/Ovary2//zf_ovary_2_6 --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1
STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Ovary2/zf_ovary_2_7_F.fq.gz /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/Ovary2/zf_ovary_2_7_R.fq.gz  --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/STARmapping/Ovary2//zf_ovary_2_7 --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1
