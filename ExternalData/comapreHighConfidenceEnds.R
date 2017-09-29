### I now want to compare the high confidence ends we get from zebrafish data with ends obtained from ulitsky etal. 


options(scipen=999)
### igor ulitsky dataset 


ends_ulitsky = read.delim("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/externalDatat/3PrimeEnds_ulitskyEtal.txt",sep="\t",stringsAsFactors = F)

chrs = lapply(strsplit(ends_ulitsky$X3P.Tag.cluster.position,":",T),function(x) x[1])
leftOverl = lapply(strsplit(ends_ulitsky$X3P.Tag.cluster.position,":",T),function(x) x[2])
chrs = do.call(c,chrs)
ends_ulitsky$leftOverl = do.call(c,leftOverl)
start = lapply(strsplit(ends_ulitsky$leftOverl,"-",T),function(x) x[1])
start = do.call(c,start)
leftOverl = lapply(strsplit(ends_ulitsky$leftOverl,"-",T),function(x) x[2])
ends_ulitsky$leftOverl = do.call(c,leftOverl)

ends = lapply(strsplit(ends_ulitsky$leftOverl," (",T),function(x) x[1])
ends = do.call(c,ends)

strand = substr(ends_ulitsky$X3P.Tag.cluster.position,start = nchar(ends_ulitsky$X3P.Tag.cluster.position)-1,stop =  nchar(ends_ulitsky$X3P.Tag.cluster.position)-1)


totalBed = cbind.data.frame(chrs,as.numeric(start),as.numeric(ends),ends_ulitsky$Category,0,strand)
totalBed = totalBed[complete.cases(totalBed),]
colnames(totalBed) = c("V1","start","end","Name","Score","V6")
#colnames(totalBed) = c("chrom","chromStart","chromEnd","name","score","strand")
totalBed$Name = paste0("peak",c(1:nrow(totalBed)))
id = paste0(totalBed$V1,totalBed$start,totalBed$end,totalBed$V6)
write.table(totalBed,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/externalDatat/3primeEnds_convertedTobed_dr9.bed",sep="\t",row.names = F,quote = F,col.names = F)
#############################################
totalBed = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/externalDatat/3primeEnds_convertedTobed_dr10.bed",sep="\t",stringsAsFactors = F)
totalBed$id = paste0(totalBed$V1,totalBed$V2,totalBed$V3,totalBed$V6)
totalBed = totalBed[!duplicated(totalBed$id),]
### ends using QuantSeq data : 
colnames(totalBed) =  c("V1","start","end","name","score","V6")
library(dplyr)
library(plyr)
library(GenomicRanges)
### the stages are : 
# 1dpf
# 256 cell stage
# 2 cell 
# 2dpf
# 4pdf
# bud
# dome
# oocyte
# sphere
# testis

path_allStages = "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/2016_12_QuantSeq_Annotaton/raw_data/"
stages = c("1dpf","256cell","2cell","2dpf","4dpf","bud","dome","oocyte","sphere","testis")
completePath = paste0(path_allStages,stages,"/output/final90percent/")

## getting the number of filtered reads 



##### comparing high confidence ends from the different stages

completePath_highConfidenceEnds = paste0(completePath,"ends_greater90percent_intergenic_n100.bed")

highConfidenceEnds = lapply(completePath_highConfidenceEnds,function(x) read.delim(x,stringsAsFactors = F,header = F))
names(highConfidenceEnds) = stages

for(i in 1:length(highConfidenceEnds)){
  temp_highConfideneEnds = highConfidenceEnds[[i]]
  temp_highConfideneEnds_plus = temp_highConfideneEnds %>% filter(V6 == "+")
  temp_highConfideneEnds_plus$start = temp_highConfideneEnds_plus$V3 -5 
  temp_highConfideneEnds_plus$end = temp_highConfideneEnds_plus$V3 + 5
  
  temp_highConfideneEnds_minus =  temp_highConfideneEnds %>% filter(V6 == "-")
  temp_highConfideneEnds_minus$end = temp_highConfideneEnds_minus$V2 + 5
  temp_highConfideneEnds_minus$start = temp_highConfideneEnds_minus$V2 - 5 
  highConfidenceEnds[[i]] = rbind(temp_highConfideneEnds_plus,temp_highConfideneEnds_minus)
}

id_highConfidenceEnds = lapply(highConfidenceEnds,function(x) paste0(x$V1,x$V2,x$V3,x$V4,x$V6))
allIntersections  = Reduce(intersect, id_highConfidenceEnds)

stages = c(stages,"ulitsky")
allCombinations = expand.grid(stages,stages)
allCombinations$intersection = NA

intersections_allCombinations =  as.data.frame(matrix(NA , nrow = length(stages),ncol = length(stages)))
allCombinations$overlaps = NA
allCombinations$queryOverlaps = NA
allCombinations$subjectOverlaps = NA
allCombinations$querySubject = NA
totalBed$start = as.numeric(as.character(totalBed$start))
totalBed$end = as.numeric(as.character(totalBed$end))
totalBed = totalBed[complete.cases(totalBed),]
highConfidenceEnds = list(highConfidenceEnds$`1dpf`,highConfidenceEnds$`256cell`,highConfidenceEnds$`2cell`,highConfidenceEnds$`2dpf`,highConfidenceEnds$`4dpf`,highConfidenceEnds$bud,highConfidenceEnds$dome,highConfidenceEnds$oocyte,highConfidenceEnds$sphere,highConfidenceEnds$testis,totalBed)
names(highConfidenceEnds) = stages
for(i in 1:nrow(allCombinations)){
  
  sample1 =  highConfidenceEnds[which(stages == allCombinations[i,1])]
  sample2 =  highConfidenceEnds[which(stages == allCombinations[i,2])]
  granges_1 = makeGRangesFromDataFrame(df = sample1,keep.extra.columns = T,seqnames.field = "V1",start.field = "start",end.field = "end",strand.field = "V6",ignore.strand = F,starts.in.df.are.0based = T )
  
  granges_2 = makeGRangesFromDataFrame(df = sample2,keep.extra.columns = T,seqnames.field = "V1",start.field = "start",end.field = "end",strand.field = "V6",ignore.strand = F,starts.in.df.are.0based = T )
  
  allCombinations$overlaps[i]= length(findOverlaps(granges_1,granges_2))
  sample1_overlap = sample1[[1]][queryHits(findOverlaps(granges_1,granges_2)),]
  sample1_overlap= sample1_overlap[!duplicated(sample1_overlap),]
  sample2_overlap = sample2[[1]][subjectHits(findOverlaps(granges_1,granges_2)),]
  sample2_overlap= sample2_overlap[!duplicated(sample2_overlap),]
  allCombinations$queryOverlaps[i] = nrow(sample1_overlap)
  allCombinations$subjectOverlaps[i] = nrow(sample2_overlap)
  allCombinations$querySubject[i] = paste(paste(names(sample1),"=",nrow(sample1_overlap)),paste(names(sample2),"=",nrow(sample2_overlap)),sep = "\n")
}
library(ggplot2)
findOverlaps = ggplot(allCombinations, aes(Var1, Var2)) + geom_tile(aes(fill = overlaps))+geom_text(aes(label = round(overlaps, 1))) + scale_fill_gradient(low = "red", high = "white") + xlab("Stage") + ylab("Stage")
print(findOverlaps)
jpeg("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/externalDataComparison/overlapWithUlitskyPeaks.jpg",height=1000,width=1200)
ggplot(allCombinations, aes(Var1, Var2)) + geom_tile(aes(fill = overlaps))+geom_text(aes(label = (querySubject))) + scale_fill_gradient(low = "white", high = "red") + xlab("Stage") + ylab("Stage")
dev.off()

