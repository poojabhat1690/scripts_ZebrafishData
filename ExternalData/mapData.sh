#!/bin/bash
#$ -pe smp 6-24

mkdir /clustertmp/pooja/mashimaTomariData/map/
ml star
cd /clustertmp/pooja/mashimaTomariData/
FASTQ=`ls *.gz`

#for i in "$FASTQ"

#do

#STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /clustertmp/pooja/mashimaTomariData/"$i" --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/mashimaTomariData/map/"$i"_mapped --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1


#done



#STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /clustertmp/pooja/mashimaTomariData/SRR2138286.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/mashimaTomariData/map/SRR2138286.fastq.gz --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1


#STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /clustertmp/pooja/mashimaTomariData/SRR2138287.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/mashimaTomariData/map/SRR2138287.fastq.gz --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1


#STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /clustertmp/pooja/mashimaTomariData/SRR2138288.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/mashimaTomariData/map/SRR2138288.fastq.gz --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1



#STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /clustertmp/pooja/mashimaTomariData/SRR2138289.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/mashimaTomariData/map/SRR2138289.fastq.gz --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1


#STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /clustertmp/pooja/mashimaTomariData/SRR2138290.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/mashimaTomariData/map/SRR2138290.fastq.gz --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1


#STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /clustertmp/pooja/mashimaTomariData/SRR2138291.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/mashimaTomariData/map/SRR2138291.fastq.gz --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1



#STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /clustertmp/pooja/mashimaTomariData/SRR2138292.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/mashimaTomariData/map/SRR2138292.fastq.gz --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1


#STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /clustertmp/pooja/mashimaTomariData/SRR2138293.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/mashimaTomariData/map/SRR2138293.fastq.gz --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1


#STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /clustertmp/pooja/mashimaTomariData/SRR2138294.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/mashimaTomariData/map/SRR2138294.fastq.gz --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1



#STAR --genomeDir /groups/ameres/bioinformatics/references/danio_rerio/dr10/ --readFilesIn /clustertmp/pooja/mashimaTomariData/SRR2138295.fastq.gz --readFilesCommand zcat --twopassMode Basic --outFileNamePrefix /clustertmp/pooja/mashimaTomariData/map/SRR2138295.fastq.gz --runThreadN 5 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1


#### now counting gene-specifically from the above created BAM files. 


##/ this counts counts from RNAseq data
ml subread


mkdir  /clustertmp/pooja/mashimaTomariData/count/

featureCounts -T 15 -B  -p -g gene_id -a  /groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10/GTFfiles/Danio_rerio.GRCz10.89.gtf -o /clustertmp/pooja/mashimaTomariData/count/countsRNAseq_mashimaTomari.txt   /clustertmp/pooja/mashimaTomariData/map//*.bam



