#### creating list for counting windows greater than 5 cpm

sampleInfoPath = read.delim("/Volumes/groups//ameres/Pooja/Projects/zebrafishAnnotation/slamSeq/dataDescription.txt",header = F)
library(DESeq2)

filePaths = file.path("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/slamSeq/counts/",sampleInfoPath$V1)

data_counts = lapply(filePaths,function(x) read.table(x,stringsAsFactors = F,sep="\t",header = T))
metaData = data_counts[[1]][,c(1:6)]

onlyCounts = lapply(data_counts,function(x) x$ReadsCPM)

onlyCounts = do.call(cbind,onlyCounts)
onlyCounts = cbind.data.frame(metaData,onlyCounts)
library(GGally)

colnames(onlyCounts) = c(colnames(metaData),as.character(sampleInfoPath$V4))
onlyCounts = onlyCounts[-which(apply(onlyCounts[,c(7:ncol(onlyCounts))],1,max)<5),]



write.table(onlyCounts,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/countingWindowsGreaterThan5cpm.bed",sep="\t",row.names = F,quote = F)

