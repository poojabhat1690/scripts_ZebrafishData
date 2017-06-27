## Including Plots

sampleInfoPath = read.delim("/Volumes/groups-1//ameres/Pooja/Projects/zebrafishAnnotation/slamSeq/dataDescription.txt",header = F)
sampleInfoPath = sampleInfoPath[1:13,]

filePaths = file.path("/Volumes/groups-1//ameres/Pooja/Projects/zebrafishAnnotation/slamSeq/counts/",sampleInfoPath$V1)

data_counts = lapply(filePaths,function(x) read.table(x,stringsAsFactors = F,sep="\t",header = T))
names(data_counts) = sampleInfoPath$V4

### separating into injection and control

data_counts_injection = data_counts[grep("injection",names(data_counts))]
data_counts_control = data_counts[grep("control",names(data_counts))]

### get the CPMS and tc counts of the injection

RPM_injection = do.call(cbind,lapply(data_counts_injection,function(x) x$ReadsCPM))
tcConversions_injection = do.call(cbind,lapply(data_counts_injection,function(x) x$ConversionRate))

##### creating summarized experiments of both these :

library(SummarizedExperiment)


metadata = data_counts_injection[[1]][,c(1:6)]

data_counts_injection_readsCPM  = SummarizedExperiment(assays = RPM_injection,rowData = GRanges(metadata))
data_counts_injection_tcConversions  = SummarizedExperiment(assays = tcConversions_injection,rowData = GRanges(metadata))

index_greaterthan5 = which(apply(assay(data_counts_injection_readsCPM),1,max)>5)
data_counts_injection_readsCPM = data_counts_injection_readsCPM[index_greaterthan5,]
data_counts_injection_tcConversions = data_counts_injection_tcConversions[index_greaterthan5,]



#### clustering the expression data : 
### to do this, i first need to normalize. each gene will be normalised to the time point at which it is namimally exoressed.
### this will help maintain shape of expression. 


max_gene = apply(assay(data_counts_injection_readsCPM),1,max)
assay(data_counts_injection_readsCPM) = assay(data_counts_injection_readsCPM)/max_gene
nclus = 3
clusterRPMs = kmeans(assay(data_counts_injection_readsCPM),centers = nclus)
names_clus = paste0("clus",c(1:nclus))
diff_clusters= vector("list",nclus)
names(diff_clusters) = names_clus

for(i in 1:nclus){
   
  diff_clusters[[i]] <- data_counts_injection_readsCPM[which(clusterRPMs$cluster==i),]
  }


#### it looks like there are 3 clusters based on thw quantSeq expression levels - 1. maternally provided transcripts that will devrease
#### transcripts that can be maternal or zygotic - so the expression should remain constant 3- increase in zygotically expressed transcripts. 





data_counts_injection = data_counts[grep("injection",names(data_counts))]
data_counts_control = data_counts[grep("control",names(data_counts))]
data_counts_injection_1 = data_counts[grep("injection",names(data_counts))]
### get the CPMS of the injection

data_counts_injection = do.call(cbind,lapply(data_counts_injection,function(x) x$ReadsCPM))
lessTHan5cpm = which(apply(data_counts_injection,1,max)<5)

data_conversion_injection_1 = do.call(cbind,lapply(data_counts_injection_1,function(x) x$ConversionRate))
data_conversion_injection_1 = data_conversion_injection_1[-lessTHan5cpm,]

max_tcCounts = apply(data_conversion_injection_1,1,max)
data_conversion_injection_1 = data_conversion_injection_1/max_tcCounts
data_conversion_injection_1= data_conversion_injection_1[complete.cases(data_conversion_injection_1),]

clusterRPMs =kmeans(data_conversion_injection_1,centers = 20,iter.max = 20)

clus1 = data_conversion_injection_1[which(clusterRPMs$cluster==1),]
clus2 = data_conversion_injection_1[which(clusterRPMs$cluster==2),]
clus3 = data_conversion_injection_1[which(clusterRPMs$cluster==3),]
clus4 = data_conversion_injection_1[which(clusterRPMs$cluster==4),]
clus5 = data_conversion_injection_1[which(clusterRPMs$cluster==5),]
clus6 = data_conversion_injection_1[which(clusterRPMs$cluster==6),]
clus7 = data_conversion_injection_1[which(clusterRPMs$cluster==7),]
clus8 = data_conversion_injection_1[which(clusterRPMs$cluster==8),]
clus9 = data_conversion_injection_1[which(clusterRPMs$cluster==9),]
clus10 = data_conversion_injection_1[which(clusterRPMs$cluster==10),]
clus11 = data_conversion_injection_1[which(clusterRPMs$cluster==11),]
clus12 = data_conversion_injection_1[which(clusterRPMs$cluster==12),]
clus13 = data_conversion_injection_1[which(clusterRPMs$cluster==13),]
clus14 = data_conversion_injection_1[which(clusterRPMs$cluster==14),]
clus15 = data_conversion_injection_1[which(clusterRPMs$cluster==15),]
clus16 = data_conversion_injection_1[which(clusterRPMs$cluster==16),] 
clus17 = data_conversion_injection_1[which(clusterRPMs$cluster==17),] 
clus18 = data_conversion_injection_1[which(clusterRPMs$cluster==18),] 
clus19 = data_conversion_injection_1[which(clusterRPMs$cluster==19),] 
clus20 = data_conversion_injection_1[which(clusterRPMs$cluster==20),] 



matplot(t(clus2),type="l")
matplot(t(clus1),type="l")
matplot(t(clus3),type="l")
matplot(t(clus4),type="l")
matplot(t(clus5),type="l")
matplot(t(clus6),type="l")
matplot(t(clus7),type="l")
matplot(t(clus8),type="l")
matplot(t(clus9),type="l")
matplot(t(clus10),type="l")
matplot(t(clus11),type="l")
matplot(t(clus12),type="l")
matplot(t(clus13),type="l")
matplot(t(clus14),type="l")
matplot(t(clus15),type="l")
matplot(t(clus16),type="l")
matplot(t(clus17),type="l")
matplot(t(clus18),type="l")
matplot(t(clus19),type="l")
matplot(t(clus20),type="l")



### trying to do this in a backgtound subtracted
```{r}


data_counts_injection = data_counts[grep("injection",names(data_counts))]
data_counts_control = data_counts[grep("control",names(data_counts))]
data_counts_injection_1 = data_counts[grep("injection",names(data_counts))]
### get the CPMS of the injection

data_counts_injection = do.call(cbind,lapply(data_counts_injection,function(x) x$ReadsCPM))
lessTHan5cpm = which(apply(data_counts_injection,1,max)<5)

data_conversion_injection_1 = do.call(cbind,lapply(data_counts_injection_1,function(x) x$ConversionRate))
data_conversion_injection_1 = data_conversion_injection_1[-lessTHan5cpm,]
data_conversion_injection_1 = data_conversion_injection_1 - data_conversion_injection_1[,1]
data_conversion_injection_1[which(data_conversion_injection_1<0)] <- 0

max_tcCounts = apply(data_conversion_injection_1,1,max)
data_conversion_injection_1 = data_conversion_injection_1/max_tcCounts
data_conversion_injection_1= data_conversion_injection_1[complete.cases(data_conversion_injection_1),]

clusterRPMs =kmeans(data_conversion_injection_1,centers = 20,iter.max = 20)

```