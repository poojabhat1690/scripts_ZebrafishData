library(checkmate)  
library(ggplot2)
library(reshape)
### reading in reference sequences :
## the refSeq annotation is created in merging250nts.Rmd
## the UTR exons are merged and also duplicates removed
refSeq_unmerged = read.table("/Volumes/groups//ameres/Pooja/Projects/zebrafishAnnotation/dr10/refSeq_dr10_GRCh38_20170504//processed/refSeq_mrna_utrsPresent.bed",stringsAsFactors = F)
refFlat <- read.table("/Volumes/groups//ameres/Pooja/Projects/zebrafishAnnotation/dr10/refSeq_dr10_GRCh38_20170504//refFlat.txt", stringsAsFactors = F)

chromosome_mm = paste0("chr",c(1:25))
refFlat = refFlat[is.element(set = chromosome_mm,el = refFlat$V3),]

refSeqTranscr <- merge(refFlat, refSeq_unmerged, by.x = "V2", by.y = "V7")
nrow(refSeq_unmerged)
nrow(refSeqTranscr)
refSeqTranscr <- unique(refSeqTranscr[,c("V3.x", "V5.x", "V6.x", "V4.y", "V4.x", "V2")])
colnames(refSeqTranscr) <- c("chr", "start", "end", "gid", "strand", "tid")


ensembl = read.delim("/Volumes/groups//ameres/Pooja/Projects/zebrafishAnnotation/dr10//ensembl_dr10_Ensembl_Genes_88//proteinCoding_annotatedUTRs.bed",stringsAsFactors = F,header = F)
ensemblTrans <- read.delim("//Volumes/groups//ameres/Pooja/Projects/zebrafishAnnotation/dr10//ensembl_dr10_Ensembl_Genes_88/transcriptStartsAndEnds_all.txt", stringsAsFactors = F, header = T)
head(ensemblTrans)
ensemblTrans <- merge(ensemblTrans, ensembl, by.x = "ensembl_transcript_id", by.y = "V7")
head(ensemblTrans)
ensemblTrans <- unique(ensemblTrans[,c("chromosome_name", "transcript_start", "transcript_end", "external_gene_name", "strand", "ensembl_transcript_id")])
colnames(ensemblTrans) <- c("chr", "start", "end", "gid", "strand", "tid")

allTrans = rbind(refSeqTranscr,ensemblTrans)
totalAnnotations= allTrans


##### doing the same for ends before filtering 

ends_all= read.delim("/Volumes/clustertmp/bioinfo/pooja/SLAMannotation/dr/output/final90percent/ends_greater90percent_intergenic_n100.bed",stringsAsFactors = F,header=F)


allEnds = ends_all


allEnds_plus = allEnds[which(allEnds$V6== "+"),]
allEnds_minus = allEnds[which(allEnds$V6 == "-"),]


distance_plus = c()
tmp = c()
for(i in 1:nrow(allEnds_plus)){
  tmp_annotation  = totalAnnotations[which(totalAnnotations$gid == allEnds_plus$V4[i]),]
  diffs = allEnds_plus[i,3] - tmp_annotation$end
  distance_tmp = diffs[which.min(abs(diffs))]
  distance_plus = c(distance_plus,distance_tmp)
  cat(i)
}


distance_minus = c()
for(i in 1:nrow(allEnds_minus)){
  tmp_annotation  = totalAnnotations[which(totalAnnotations$gid == allEnds_minus$V4[i]),]
  diffs =  tmp_annotation$start - allEnds_minus[i,2] ### this is inverted because if not the shorted end has a higher coordinate than the UTR and this is opposite for the positive strand
  
  distance_tmp = diffs[which.min(abs(diffs))]
  distance_minus = c(distance_minus,distance_tmp)
  cat(i)
}


allEnds_plus$distanceMin = distance_plus
allEnds_minus$distanceMin = distance_minus

allEnds = rbind(allEnds_plus,allEnds_minus)


### 

within5nt = allEnds[which(abs(allEnds$distanceMin)<=5),]
allEnds_order = allEnds[order(abs(allEnds$distanceMin),decreasing = T),]
allEnds_short = cbind.data.frame(allEnds$V4,allEnds$distanceMin,stringsAsFactors=F)
colnames(allEnds_short) = c("V1","V2")
a_split = split( allEnds_short,allEnds_short$V1,T)
allEnds_total = do.call(rbind,a_split)
allEnds_total = allEnds_total[order(allEnds_total$V2,decreasing = F),]
tabulated = as.data.frame(table(allEnds_total$V1))
allEnds_total_numberEnds = merge(allEnds_total,tabulated,by.x="V1",by.y="Var1")
allEnds_total_numberEnds= allEnds_total_numberEnds[order(allEnds_total_numberEnds$Freq , abs(allEnds_total_numberEnds$V2)),]
#allEnds_total_numberEnds$V2 =as.factor(allEnds_total_numberEnds$V2)
#allEnds_total_numberEnds= allEnds_total_numberEnds[order(allEnds_total_numberEnds$Freq ),]

allEnds_total_numberEnds$Freq = as.factor(allEnds_total_numberEnds$Freq)
allEnds_total_numberEnds$order = seq(length(allEnds_total_numberEnds$V1),1,by = -1)

allEnds_total_numberEnds_beforeFilter= allEnds_total_numberEnds

##### 
allEnds_total_numberEnds_beforeFilter.sepFreq = split(x = allEnds_total_numberEnds_beforeFilter,f = allEnds_total_numberEnds_beforeFilter$Freq,drop = T)
allEnds_total_numberEnds_beforeFilter.sepFreq = melt(lapply(allEnds_total_numberEnds_beforeFilter.sepFreq,function(x) length(unique(x$V1))))
library(ggplot2)
pdf("/Volumes/groups///ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/distanceOfEndsFromAnnotations.pdf")
p = ggplot(allEnds_total_numberEnds_beforeFilter,aes(y=reorder(V1, order),x=V2,group=Freq)) + geom_point(aes(group=Freq),colour="black",alpha=1/10)+theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
p = p+ xlab("Minimum distance of end from annotated refSeq and Ensembl 3' ends") + ylab("Genes") + ggtitle(paste0("allEnds n=",nrow(allEnds_total_numberEnds_beforeFilter))) + xlim(c(-10000,10000))
print(p)

p = ggplot(allEnds_total_numberEnds_beforeFilter,aes(y=reorder(V1, order),x=V2,group=Freq)) + geom_point(aes(group=Freq),colour="black",alpha=1/10)+theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) 
p = p+ xlab("Minimum distance of end from annotated refSeq and Ensembl 3' ends") + ylab("Genes") + geom_vline(xintercept = c(-5,5),linetype="longdash") + xlim(-50,50)+ ggtitle(paste0("allEnds zoomin n=",nrow(allEnds_total_numberEnds_beforeFilter)))
print(p)
dev.off()





refSeqTpm = read.delim("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/countingWindowsGreaterThan5cpm.bed",header = T,stringsAsFactors = F)
allEnds_total_numberEnds= allEnds_total_numberEnds[-which(allEnds_total_numberEnds$V1 == ""),]
allEnds_total_numberEnds_tpm10 = merge(allEnds_total_numberEnds,refSeqTpm,by.y="Name",by.x="V1")

library(reshape)
allEnds_total_numberEnds_beforeFilter.sepFreq = split(x = allEnds_total_numberEnds_tpm10,f = allEnds_total_numberEnds_tpm10$Freq,drop = T)
allEnds_total_numberEnds_beforeFilter.sepFreq = melt(lapply(allEnds_total_numberEnds_beforeFilter.sepFreq,function(x) length(unique(x$V1))))

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots//distanceAfterFIlter_tmpgreaterThan10.pdf")
p = ggplot(allEnds_total_numberEnds_tpm10,aes(y=reorder(V1, order),x=V2,col=Freq,group=Freq)) + geom_point(aes(group=Freq),colour="black",alpha=1/10)+theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
p = p+ xlab("Minimum distance of end from annotated refSeq and Ensembl 3' ends (CPM>10)") + ylab("Genes")
print(p)
dev.off()
