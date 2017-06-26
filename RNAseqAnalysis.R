#### looking at RNAseq data
library(SummarizedExperiment)

### 

RNAseqCounts = read.table("/Volumes/groups-1/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/rnaSeqReads_featureCounts.txt",stringsAsFactors = F,header=T)
sampleInfo = read.delim("/Volumes/groups-1/ameres/Pooja/Projects/zebrafishAnnotation/ucscBrowser/RNAseqDatasets_sampleInfo.txt",sep="\t",header=F,stringsAsFactors = F)
RNAseqCounts$Length = RNAseqCounts$Length/1000
reads_RNAseq = RNAseqCounts[,c(7:ncol(RNAseqCounts))]
#colnames(reads_RNAseq) = sampleInfo$V3
reads_RNAseq = reads_RNAseq/RNAseqCounts$Length
scalingFactor = apply(reads_RNAseq,2,sum)/1000000
reads_RNAseq  = t(reads_RNAseq)
reads_RNAseq = reads_RNAseq/scalingFactor
reads_RNAseq = as.data.frame(t(reads_RNAseq))
colnames(reads_RNAseq) = sampleInfo$V3
rownames(reads_RNAseq) = RNAseqCounts$Geneid

reads_RNAseq_5TPM = reads_RNAseq[which(apply(reads_RNAseq,1,max)>5),]
reads_RNAseq_5TPM$`256cell` = NULL
corrlation_timepoints = cor(reads_RNAseq)

corrplot::corrplot(corrlation_timepoints,type = "lower")


#### RNAseq and quantSeq correlation - I have a problem here .. the quantSeq data only has gene names but the RNAseq has ENS ids

sampleInfoPath = read.delim("/Volumes/groups-1//ameres/Pooja/Projects/zebrafishAnnotation/slamSeq/dataDescription.txt",header = F,stringsAsFactors = F)

filePaths = file.path("/Volumes/groups-1/ameres/Pooja/Projects/zebrafishAnnotation/slamSeq/counts/",sampleInfoPath$V1)

quantSeqReads = lapply(filePaths,function(x) read.table(x,header=T,stringsAsFactors = F,sep="\t"))


CPM_quantSeq = lapply(quantSeqReads,function(x) x$ReadsCPM)
CPM_quantSeq = do.call(cbind,CPM_quantSeq)
#CPM_quantSeq = as.data.frame(CPM_quantSeq)
metadata = quantSeqReads[[1]][,c(1,2,3,6,4)]
#CPM_quantSeq = cbind.data.frame(metadata,CPM_quantSeq)
metadata_granges = GRanges(metadata)
colnames(CPM_quantSeq) = sampleInfoPath$V4

CPM_quantSeq_summarizedExpt = SummarizedExperiment(assays=list(counts=CPM_quantSeq),rowData=metadata_granges, colData=sampleInfoPath)
CPM_quantSeq = cbind.data.frame(metadata,CPM_quantSeq)
library(dplyr)
###### from bioconductor get the gene names
CPM_quantSeq_split = split(CPM_quantSeq,CPM_quantSeq$Name,T)
CPM_quantSeq_split = lapply(CPM_quantSeq_split, function(x) colSums(x[,c(6:ncol(x))]))
CPM_quantSeq_split = do.call(rbind,CPM_quantSeq_split)
CPM_quantSeq_split = as.data.frame(CPM_quantSeq_split)

##### converting RNAseq enssembl Ids to the gene names : 
library(biomaRt)

listMarts(host = "www.ensembl.org")

ensembl  = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="drerio_gene_ensembl",host="www.ensembl.org")
geneNames_id = getBM(attributes = c("ensembl_gene_id","external_gene_name","gene_biotype"),values = RNAseqCounts$Geneid,filters="ensembl_gene_id",mart = ensembl)
geneNames_id = geneNames_id[which(geneNames_id$gene_biotype=="protein_coding"),]
reads_RNAseq_merged = merge(reads_RNAseq,y =geneNames_id,by.x="row.names",by.y="ensembl_gene_id")

quantSeqAndRnaSeq = merge(CPM_quantSeq_split,reads_RNAseq_merged,by.x="row.names",by.y="external_gene_name",all=T)
quantSeqAndRnaSeq = quantSeqAndRnaSeq[-1,] 
quantSeqAndRnaSeq$Row.names = NULL
quantSeqAndRnaSeq$gene_biotype = NULL
quantSeqAndRnaSeq  = quantSeqAndRnaSeq[complete.cases(quantSeqAndRnaSeq),]
quantSeqAndRnaSeq_5cpm = quantSeqAndRnaSeq[which(apply(quantSeqAndRnaSeq,1,max)>100),]
allCorrelations = cor(quantSeqAndRnaSeq_5cpm,use = "complete.obs",method = "pearson")
corrplot::corrplot(allCorrelations)
