#CELLLINE="mESC"
#CELLLINE="MEF"
#CELLLINE="hESC"
#CELLLINE="Hela"
CELLLINE="dr"


PIPELINE="/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/scripts_ZebrafishData/pipeline/"
ANNOBASE="/groups/ameres/Pooja/Projects/zebrafishAnnotation/"
GBASE="/groups/ameres/bioinformatics/references/danio_rerio/dr10//"

if [ $CELLLINE = "dr" ]; then
    BIN="//groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10_data/"
    GBUILD="dr10"
    GET_ANNOTATION="FALSE"
    BOUT="/clustertmp/bioinfo/pooja/SLAMannotation/$CELLLINE/output"
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

ls $QUANT/*.fq.gz | perl -pe "s#$QUANT/##" > $QUANT_ALIGN/sampleInfo.txt
TASKS=`wc -l $QUANT_ALIGN/sampleInfo.txt | awk '{print $1}'`

export PIPELINE
export QUANT
export QUANT_ALIGN
export GENOME

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

#extract +1 to +20 A fraction
module unload R
module load R/3.3.0

TMP="$QUANT_ALIGN"/n_100_global_a0/

MINUSSUM=`awk '{ sum += $1 } END { print sum } '  "$TMP"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique.bed`
PLUSSUM=`awk '{ sum += $1 } END { print sum } '  "$TMP"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique.bed`
TOTALSUM=`expr $PLUSSUM + $MINUSSUM`
THRESHOLD=`expr  $"$TOTALSUM" / 1000000`



rmd="$PIPELINE/OverlappingPrimingSitesWithAnnotations/sequencesForNucleotideProfile.R"
INPUT=$QUANT_ALIGN/n_100_global_a0/
Rscript --vanilla -e "InPath='$INPUT';threshold="$THRESHOLD";source('$rmd')"


mkdir -p $BOUT/PASplot
OUTPUT=$BOUT/PASplot
rmd="$PIPELINE/OverlappingPrimingSitesWithAnnotations/nucleotideProfiles_markdown.new.R"

Rscript --slave -e "PPath='$PIPELINE'; InPath='$INPUT'; OutPath='$OUTPUT'; ucscDir='$UCSC'; ensemblDir='$ENSEMBL';source('$rmd')"


OUTPUT=$BOUT/intergenicPeaks
mkdir -p $OUTPUT

rmd="$PIPELINE/intergenicPeaks/getLongestUTR.R"
Rscript --slave -e "PPath='$PIPELINE'; InPath='$INPUT'; OutPath='$OUTPUT'; ucscDir='$UCSC'; ensemblDir='$ENSEMBL';source('$rmd')"


source $PIPELINE/intergenicPeaks/getClosestGene.sh



#make sure reads are on opposite strand "-S" !!!
mkdir $BOUT/ExtendingINtergenicRegions
mkdir $BOUT/coverage
Rscript --vanilla -e "BIn='$BIN'; BOut='$BOUT'; ucscDir='$UCSC'; ensemblDir='$ENSEMBL'; source('$PIPELINE/intergenicPeaks/addingIntergenicPeaks.R')"


rmd="$PIPELINE/90PercentFiltering_merging_countingWindows/assignToUTRs.R"
Rscript --slave -e "BIn='$BIN'; BOut='$BOUT'; ucscDir='$UCSC'; ensemblDir='$ENSEMBL';source('$rmd')"

rmd="$PIPELINE/90PercentFiltering_merging_countingWindows/90PercentFiltering.R"
Rscript --slave -e "BIn='$BIN'; BOut='$BOUT'; ucscDir='$UCSC'; ensemblDir='$ENSEMBL';source('$rmd')"

Rscript --vanilla -e "BIn='$BIN'; BOut='$BOUT'; ucscDir='$UCSC'; ensemblDir='$ENSEMBL'; fai='$GENOME.fai'; source('$PIPELINE/90PercentFiltering_merging_countingWindows/mergingCounting.new.R')"

find $BOUT* | egrep ".bed$" | egrep -v "final90percent" | xargs -P 1 gzip

mkdir $BIN/../output
rsync -rva  --exclude "*.bam" --exclude "*.fastq" --exclude "*.fq" $BOUT/* $BIN/../output
