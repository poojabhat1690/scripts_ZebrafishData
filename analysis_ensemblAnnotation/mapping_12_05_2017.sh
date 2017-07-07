#!/bin/bash


#$ -t 1-24

#$ -pe smp 6-24
#$ -l hostname=!compute-6-3&!compute-6-14&!compute-6-18&!compute-6-20&!compute-6-3*&!compute-6-4*&!compute-6-5*


########$ -wd /home/<your_UCL_id>/Scratch/output

# 8. Parse parameter file to get variables.
number=$SGE_TASK_ID
paramfile=/groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10_data/quantseq/sampleInfoIndex.txt
 
index=`sed -n ${number}p $paramfile | awk '{print $1}'`
variable1=`sed -n ${number}p $paramfile | awk '{print $2}'`
USENAME=`sed -n ${number}p $paramfile | awk '{print $1}'`
####variable2=`sed -n ${number}p $paramfile | awk '{print $3}'`
##variable3=`sed -n ${number}p $paramfile | awk '{print $4}'`

ml cutadapt



mkdir -p /clustertmp/pooja/mapping_dr10_12052017/adapterTrimmed/
OUTDIR=/clustertmp/pooja/mapping_dr10_12052017/adapterTrimmed/

##cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG -o "$OUTDIR"/"$index"_adapterTrimmed.fastq -m 18 --trim-n //groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10_data/quantseq/allData/quantseq/"$index"



module load python
pip install --user IntervalTree
module load joblib
module load pysam
module load R/3.2.2
module load samtools/1.3.1



/groups/ameres/bioinformatics/tools/slamdunk-new-1/slamdunk/bin/slamdunk all -r /groups/ameres/bioinformatics/references/danio_rerio/dr10/danRer10.fa -o "$OUTDIR" -b /groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/analysis_ensemblAnnotation/countingWindoes_ensemblAnnotation.bed  -fb  /groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10/ensembl_dr10_Ensembl_Genes_88/proteinCoding_annotatedUTRs.bed -a 5 -n 100 -t 15 -mq 0 -mi 0.95 -m -l -rl 95 "$OUTDIR"/"$index"_adapterTrimmed.fastq 
