library(reshape)
library(ggplot2)
library(SummarizedExperiment)
library(biomaRt)
sampleInfoPath = read.delim("/Volumes/groups-1//ameres/Pooja/Projects/zebrafishAnnotation/slamSeq/dataDescription.txt",header = F)
sampleInfoPath = sampleInfoPath[which(sampleInfoPath$V7 == "yes"),]

filePaths = file.path("/Volumes/groups-1//ameres/Pooja/Projects/zebrafishAnnotation/slamSeq/counts/",sampleInfoPath$V1)

data_counts = lapply(filePaths,function(x) read.table(x,stringsAsFactors = F,sep="\t",header = T))
names(data_counts) = sampleInfoPath$V6

### get the CPMS of all 

RPM_control= do.call(cbind,lapply(data_counts,function(x) x$ReadsCPM))
RPM_control = as.data.frame(RPM_control)

RPM_control  = as.data.frame(t(RPM_control))
RPM_control$condition = rownames(RPM_control)
RPM_control_split = split(RPM_control,RPM_control$condition,T)
RPM_control_split = lapply(RPM_control_split,function(x) colMeans(x[,c(1:(ncol(x))-1)]))
RPM_control_combined = do.call(rbind,RPM_control_split)
RPM_control_combined = (t(RPM_control_combined))

###
metadata = data_counts[[1]][,c(1:6)]

listMarts(host = "www.ensembl.org")

ensembl  = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="drerio_gene_ensembl",host="www.ensembl.org")

geneName_ensemblIds = getBM(attributes = c("external_gene_name","ensembl_gene_id"),filters = "external_gene_name",values = metadata$Name,mart = ensembl)
colnames(geneName_ensemblIds) = c("Name","ensembl_gene_id")
library(plyr)
metadata_combine = join(metadata,geneName_ensemblIds,match="first")

metadata = makeGRangesFromDataFrame(metadata_combine,keep.extra.columns = T,ignore.strand = F,starts.in.df.are.0based = T,seqnames.field = "Chromosome",start.field = "Start",end.field = "End",strand.field = "Strand")
rownames(RPM_control_combined) = metadata$ensembl_gene_id
     
RPM_control_combined = SummarizedExperiment(assays = RPM_control_combined,rowData = (metadata))
library(dplyr)

RPM_control_combined_split = split(as.data.frame(assay(RPM_control_combined)),rownames(RPM_control_combined),T)
RPM_control_combined_split = lapply(RPM_control_combined_split,function(x) colMeans(x))
RPM_control_combined_perGene = do.call(rbind,RPM_control_combined_split)

max_RPM = which(apply((RPM_control_combined_perGene),1,max)>5)
RPM_control_combined = RPM_control_combined_perGene[max_RPM,]
correlation_samples = cor((RPM_control_combined_perGene))
corrplot::corrplot(correlation_samples)


#### reading in the rnaSeq counts produced using featureCounts


RNAseqCounts = read.table("/Volumes/groups-1/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/rnaSeqReads_featureCounts_strandSpecific_allGenes.txt",sep="\t",stringsAsFactors = F,header = T)
colnames(RNAseqCounts)
colnames(RNAseqCounts) = c(colnames(RNAseqCounts[1:6]),"RNAseq1KCell","RNAseq2-4cell","RNAseq2-4cell","RNAseq256cell","RNAseq28hpf","RNAseq28hpf","RNAseq2dpf","RNAseq2dpf","RNAseq5dpf","RNAseq5dpf","RNAseqBud","RNAseqBud","RNAseqDome","RNAseqDome","RNAseq1Kcell","RNAseqOblong","RNAseqShield","RNAseqShield","RNAseqShield","RNAseqOocyte","RNAseqOocyte","RNAseqOocyte","RNAseqOocyte","RNAseqOocyte","RNAseqOocyte","RNAseqOvary","RNAseqOvary","RNAseqOvary","RNAseqSperm","RNAseqSperm","RNAseqTestis","RNAseqTestis","RNAseqTestis","RNAseqTestis","RNAseqTestis","RNAseqTestis")
RNAseq_onlyCounts = RNAseqCounts[,c(7:ncol(RNAseqCounts))]
RNAseqCounts$Length = RNAseqCounts$Length/1000



tpm <- function(counts, lengths) {  
  rate <- counts /lengths
  rate / sum(rate) * 1e6
}

RNAseq_TPM = apply(RNAseq_onlyCounts,2,function(x) tpm(x,RNAseqCounts$Length))
metadata_RNAseq = RNAseqCounts[,c(1,6)]
RNAseq_TPM = as.data.frame(t(RNAseq_TPM))
RNAseq_TPM$timepoints = c("RNAseq1KCell","RNAseq2-4cell","RNAseq2-4cell","RNAseq256cell","RNAseq28hpf","RNAseq28hpf","RNAseq2dpf","RNAseq2dpf","RNAseq5dpf","RNAseq5dpf","RNAseqBud","RNAseqBud","RNAseqDome","RNAseqDome","RNAseq1Kcell","RNAseqOblong","RNAseqShield","RNAseqShield","RNAseqShield","RNAseqOocyte","RNAseqOocyte","RNAseqOocyte","RNAseqOocyte","RNAseqOocyte","RNAseqOocyte","RNAseqOvary","RNAseqOvary","RNAseqOvary","RNAseqSperm","RNAseqSperm","RNAseqTestis","RNAseqTestis","RNAseqTestis","RNAseqTestis","RNAseqTestis","RNAseqTestis")
RNAseq_TPM_split = split(RNAseq_TPM,RNAseq_TPM$timepoints,T)
RNAseq_TPM_split = lapply(RNAseq_TPM_split,function(x) colMeans(x[,c(1:(ncol(x))-1)]))
RNAseq_TPM_split_combined = do.call(rbind,RNAseq_TPM_split)
RNAseq_TPM_split_combined = (t(RNAseq_TPM_split_combined))
rownames(RNAseq_TPM_split_combined) = RNAseqCounts$Geneid

quantSeq_RNAseqCombined = merge(RPM_control_combined,RNAseq_TPM_split_combined,by="row.names")

quantSeq_RNAseqCombined$Row.names = NULL
quantSeq_RNAseqCombined_cor = cor(quantSeq_RNAseqCombined,method = "spearman")
corrplot::corrplot(quantSeq_RNAseqCombined_cor)
