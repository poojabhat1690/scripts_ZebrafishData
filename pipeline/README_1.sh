#CELLLINE="mESC"
#CELLLINE="MEF"
#CELLLINE="hESC"
#CELLLINE="Hela"
CELLLINE="dr"


PIPELINE="/groups/ameres/Pooja/Projects/zebrafishAnnotation/pipeline/"
ANNOBASE="/groups/ameres/Pooja/Projects/zebrafishAnnotation/"
GBASE="/groups/ameres/bioinformatics/references/danio_rerio/dr10/"

if [ $CELLLINE = "dr" ]; then
    BIN="//groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10_data/quantseq/allData/"
    GBUILD="dr10"
    GET_ANNOTATION="FALSE"
    BOUT="/clustertmp/bioinfo/pooja/SLAMannotation/$CELLLINE//output/"
else if [ $CELLLINE = "MEF" ]; then
	 BIN="/groups/bioinfo/shared/Ameres.GRP/Veronika.Herzog/SLAMannotation/$CELLLINE/input"
	 GBUILD="mm10"
	 GET_ANNOTATION="FALSE"
	 BOUT="/clustertmp/bioinfo/thomas/SLAMannotation/$CELLLINE/output"
     else if [ $CELLLINE = "hESC" ]; then
	      BIN="/groups/bioinfo/shared/Ameres.GRP/Veronika.Herzog/SLAMannotation/$CELLLINE/input"
	      GBUILD="hg38"
	      GET_ANNOTATION="FALSE"
	      BOUT="/clustertmp/bioinfo/thomas/SLAMannotation/$CELLLINE/output"
	  else if [ $CELLLINE = "Hela" ]; then
		   BIN="/groups/bioinfo/shared/Ameres.GRP/Veronika.Herzog/SLAMannotation/$CELLLINE/input"
		   GBUILD="hg38"
		   GET_ANNOTATION="FALSE"
		   BOUT="/clustertmp/bioinfo/thomas/SLAMannotation/$CELLLINE/output"
	       fi
	  fi
     fi
fi

mkdir -p $BOUT


#download UCSC annotation as in /groups/bioinfo/shared/Ameres.GRP/incoming/3PrimeEndAnnotation/annotations_mm10/gettingAnnotations.txt
if [ $GBUILD = "mm10" ]; then

    ENSEMBL="$ANNOBASE/$GBUILD/ensembl_mm10_Ensembl_Genes_87_06-march-2017/"
    UCSC="$ANNOBASE/$GBUILD/refSeq_mm10_GRCm38_06-march-2017/"
    GENOME="$GBASE/mmu_mm10_whole_genome.fa"
    
    if [ $GET_ANNOTATION = "TRUE" ]; then
	rmd="$ANNOBASE/$GBUILD/getAnnotations.Rmd"
	mkdir -p $ENSEMBL
	mkdir -p $UCSC/processed
	Rscript --slave -e "ucscDir='$UCSC'; ensemblDir='$ENSEMBL';rmarkdown::render('$rmd', output_file='$ANNOBASE/$GBUILD/getAnnotations.html')"
    fi

else if [ $GBUILD = "dr10" ]; then

	 ENSEMBL="$ANNOBASE/$GBUILD/ensembl_dr10_Ensembl_Genes_88/"
	 UCSC="$ANNOBASE/$GBUILD/refSeq_dr10_GRCh38_20170504/"
	 GENOME="$GBASE/danRer10.fa"
    
	 if [ $GET_ANNOTATION = "TRUE" ]; then
	     rmd="$ANNOBASE/$GBUILD/getAnnotations.Rmd"
	     mkdir -p $ENSEMBL
	     mkdir -p $UCSC/processed
	     Rscript --slave -e "ucscDir='$UCSC'; ensemblDir='$ENSEMBL';rmarkdown::render('$rmd', output_file='$ANNOBASE/$GBUILD/getAnnotations.html')"
	 fi
     fi
     

fi



#Raw fastq.gz; Preprocess/map/postprocess
QUANT=$BIN/quantseq
QUANT_ALIGN=$BOUT/polyAmapping_allTimepoints
mkdir $QUANT_ALIGN
cd $QUANT_ALIGN

ls $QUANT/*.fastq.gz | perl -pe "s#$QUANT/##" > $QUANT_ALIGN/sampleInfo.txt
TASKS=`wc -l $QUANT_ALIGN/sampleInfo.txt | awk '{print $1}'`

export PIPELINE
export QUANT
export QUANT_ALIGN
export GENOME

qsub -sync y -t 1:$TASKS $PIPELINE/pre-processing/beforeMapping.new.sh
qsub -sync y -t 1:$TASKS $PIPELINE/primingSites/mapping_n100_global.sh
qsub -sync y $PIPELINE/primingSites/afterMapping.sh

#cd n_100_global_a0
#samtools merge merged.bam *filtered.bam
#samtools index merged.bam
#cd ..
#source $PIPELINE/mapAllReads.sh
#cd allReads
#samtools merge merged.bam *filtered.bam
#samtools index merged.bam
#cd ..

mkdir logs
mv *.sh.* logs
mv map.[eop]* logs
mv cutadapt.[eo]* logs

cd $BOUT
