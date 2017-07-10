#!/bin/bash
#$ -S /bin/bash
#$ -q public.q
#$ -cwd
#$ -pe smp 6-24
#$ -l hostname=!compute-6-3&!compute-6-14&!compute-6-18&!compute-6-20&!compute-6-3*&!compute-6-4*&!compute-6-5*
#$ -v GENOME
#$ -v QUANT_ALIGN
#$ -v PIPELINE
### mapping

hostname

#set dirs
PARAMETER="$QUANT_ALIGN/sampleInfo.txt"


INPUT=$QUANT_ALIGN
OUTDIR=$INPUT/n_100_global_a0/
mkdir -p $OUTDIR

#get input file (index)
COUNTER=0
for index in `cat $PARAMETER`; do #parameter file = list of files
    COUNTER=$((COUNTER+1))
    if [ $COUNTER = $SGE_TASK_ID ]; then
	break
    fi
done

#cat "$OUTDIR"/*_polyAreads_polyAremoved.fastq >"$OUTDIR"/polyAreads_polyAremoved_pooled.fastq

module purge
module load minimal
module load biopython/1.65-python2.7.3
module load joblib/0.9.4-python2.7.3
module load numpy/1.11.0
module load cmake/3.2.2
module load pandas/0.18.1
module load pysam/0.8.3


### just  to compare
#mkdir /clustertmp/pooja/polyAmapping_allTimepoints//n_100_global_a0/


#GLOBAL alignment in that version default!!!! (change if version changes)

$PIPELINE/slamdunk/bin/slamdunk -v

$PIPELINE/slamdunk/bin/slamdunk map -r $GENOME -o $OUTDIR/ -n 100 -5 0 -t $NSLOTS -a 0 -e $INPUT/${index}_5primetrimmed_trimmed_sizefiltered_polyAreads_polyAremoved.fastq

$PIPELINE/slamdunk/bin/slamdunk filter -o $OUTDIR -mq 0 -mi 0.95 -t $NSLOTS $OUTDIR/${index}_5primetrimmed_trimmed_sizefiltered_polyAreads_polyAremoved_slamdunk_mapped.bam
