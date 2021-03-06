### Manually started. Used to align all quantseq reads globally. No polyA filtering (grep).


module load gridengine

### just  to compare
#mkdir /clustertmp/pooja/polyAmapping_allTimepoints//n_100_global_a0/


#GLOBAL alignment in that version default!!!! (change if version changes)

cd $QUANT_ALIGN

#PIPELINE="/groups/bioinfo/shared/Ameres.GRP/Veronika.Herzog/SLAMannotation/pipeline/"

mkdir allReads
for f in *_5primetrimmed_trimmed.fastq; do
    echo $f

    qsub -q public.q -cwd -b y -shell y -N cutadapt -sync y -l 'hostname=!compute-6-3&!compute-6-14&!compute-6-18&!compute-6-20&!compute-6-3*&!compute-6-4*&!compute-6-5*' "module load cutadapt; cutadapt --no-indels -m 18 -e 0 -a 'A{1000}'  -o 'allReads/'${f/.fastq/_polyAremoved.fastq}  $f"
    f="allReads/"${f/.fastq/_polyAremoved.fastq}
    
    qsub -q public.q -cwd -b y -shell y -N map -pe smp 6-24 -l 'hostname=!compute-6-3&!compute-6-14&!compute-6-18&!compute-6-20&!compute-6-3*&!compute-6-4*&!compute-6-5*' "module purge; module load minimal; module load biopython/1.65-python2.7.3; module load joblib/0.9.4-python2.7.3; module load numpy/1.11.0; module load cmake/3.2.2; module load pandas/0.18.1; module load pysam/0.8.3; $PIPELINE/slamdunk/bin/slamdunk map -r $GENOME -o allReads/ -n 100 -5 0 -t \$NSLOTS -a 0 -e $f; $PIPELINE/slamdunk/bin/slamdunk filter -o allReads -mq 0 -mi 0.95 -t \$NSLOTS ${f/.fastq/_slamdunk_mapped.bam}"

done
