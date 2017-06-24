#### looking at RNAseq data

### 

RNAseqCounts = read.table("/Volumes/groups-1/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/rnaSeqReads_featureCounts.txt",stringsAsFactors = F,header=T)
sampleInfo = read.delim("/Volumes/groups-1/ameres/Pooja/Projects/zebrafishAnnotation/ucscBrowser/RNAseqDatasets_sampleInfo.txt",sep="\t",header=F)
RNAseqCounts$Length = RNAseqCounts$Length/1000
reads_RNAseq = RNAseqCounts[,c(7:ncol(RNAseqCounts))]
reads_RNAseq = reads_RNAseq/RNAseqCounts$Length
scalingFactor = apply(reads_RNAseq,2,sum)/1000000
reads_RNAseq  = t(reads_RNAseq)
reads_RNAseq = reads_RNAseq/scalingFactor
reads_RNAseq = as.data.frame(t(reads_RNAseq))
colnames(reads_RNAseq) = sampleInfo$V3
corrlation_timepoints = cor(reads_RNAseq)
corrplot::corrplot(corrlation_timepoints,type = "lower")
