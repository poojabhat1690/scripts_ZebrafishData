## Including Plots
library(reshape)
library(ggplot2)
library(SummarizedExperiment)

sampleInfoPath = read.delim("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10_data/quantseq/sampleInfoIndex.txt",sep="\t",header = F,stringsAsFactors = F)
sampleInfoPath = sampleInfoPath[11:24,]
sampleInfoPath$countFile = paste0(sampleInfoPath$V1,"_adapterTrimmed_slamdunk_mapped_filtered_tcount.tsv")
filePaths = file.path("/Volumes//groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/comparingAnnotationsBetweenStages/count/",sampleInfoPath$countFile)

data_counts = lapply(filePaths,function(x) read.table(x,stringsAsFactors = F,sep="\t",header = T))
names(data_counts) = sampleInfoPath$V3

### separating into injection and control

data_counts_injection = data_counts[grep("4su",names(data_counts))]
data_counts_control = data_counts[grep("con",names(data_counts))]

### get the CPMS and tc counts of the injection

RPM_injection = do.call(cbind.data.frame,lapply(data_counts_injection,function(x) x$ReadsCPM))
names(RPM_injection) = names(data_counts_injection)

RPM_injection = cbind.data.frame(RPM_injection$`128_4su`,RPM_injection$`30min_4su`,RPM_injection$`65min_4su`,RPM_injection$`105min_4su`,RPM_injection$`140min_4su`,RPM_injection$`170min_4su`)
colnames(RPM_injection) = c("0min","30min","65min","105min","140min","170min")
tcConversions_injection = do.call(cbind.data.frame,lapply(data_counts_injection,function(x) x$ConversionRate))
names(tcConversions_injection) = names(data_counts_injection)
tcConversions_control = do.call(cbind.data.frame,lapply(data_counts_control,function(x) x$ConversionRate))
names(tcConversions_control) = names(data_counts_control)


### rearrange the samples based on time : 

tcConversions_control = cbind.data.frame(tcConversions_control$`128_con`,tcConversions_control$`30min_con`,tcConversions_control$`65min_con`,tcConversions_control$`105min_con`,tcConversions_control$`140min_con`,tcConversions_control$`170min_con`)
colnames(tcConversions_control) = c("0min","30min","65min","105min","140min","170min")

tcConversions_injection = cbind.data.frame(tcConversions_injection$`128_4su`,tcConversions_injection$`30min_4su`,tcConversions_injection$`65min_4su`,tcConversions_injection$`105min_4su` ,tcConversions_injection$`140min_4su`,tcConversions_injection$`170min_4su`)
colnames(tcConversions_injection) = c("0min","30min","65min","105min","140min","170min")

#### first let us plot control v/s injections

injectionMeanConversion = colMeans(tcConversions_injection)
controlMeanConversion = colMeans(tcConversions_control)
timePoints = names(injectionMeanConversion)
meanConversions = cbind.data.frame(injectionMeanConversion,controlMeanConversion)
colnames(meanConversions) = c("injection","control")
meanConversions = reshape::melt(meanConversions)
meanConversions$time= timePoints
meanConversions$timepoints = c(30,60,105,128,140,170)

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/backGroundTCvsInjection.pdf",height = 7,width=10)
ggplot(meanConversions,aes(x=timepoints,y=value,group=variable)) + geom_line(aes(col=variable),size=1.25) + theme_bw() + xlab("time (min)") + ylab("mean T>C conversions")+theme(axis.text=element_text(size=12),axis.title=element_text(size=14),legend.text=element_text(size=14))
dev.off()

tcConversions_injection  = as.matrix(tcConversions_injection) - as.matrix( tcConversions_control)
tcConversions_injection[which(tcConversions_injection<0)]<-0


##### creating summarized experiments of both these :



library(GenomicRanges)
metadata = data_counts_injection[[1]][,c(1:6)]
metadata$id = paste0(data_counts_injection[[1]]$Chromosome,data_counts_injection[[1]]$Start,data_counts_injection[[1]]$End,data_counts_injection[[1]]$Name,data_counts_injection[[1]]$Strand)
metadata = makeGRangesFromDataFrame(metadata,keep.extra.columns = T,ignore.strand = F,starts.in.df.are.0based = T,seqnames.field = "Chromosome",start.field = "Start",end.field = "End",strand.field = "Strand")

data_counts_injection_readsCPM  = SummarizedExperiment(assays = as.matrix(RPM_injection),rowData = metadata)
data_counts_injection_tcConversions  = SummarizedExperiment(assays = tcConversions_injection,rowData = metadata)

index_greaterthan5 = which(apply(assay(data_counts_injection_readsCPM),1,max)>5)
data_counts_injection_readsCPM = data_counts_injection_readsCPM[index_greaterthan5,]
data_counts_injection_tcConversions = data_counts_injection_tcConversions[index_greaterthan5,]


# test to check overlap  
# a = as.data.frame(rowRanges(data_counts_injection_readsCPM))
# b= as.data.frame(rowRanges(data_counts_injection_tcConversions))
# a$id = paste0(a[,1],a[,2],a[,3],a[,4],a[,5])
# b$id = paste0(b[,1],b[,2],b[,3],b[,4],b[,5])

#### clustering the expression data : 
### to do this, i first need to normalize. each gene will be normalised to the time point at which it is namimally exoressed.
### this will help maintain shape of expression. 


max_gene = apply(assay(data_counts_injection_readsCPM),1,max)
assay(data_counts_injection_readsCPM) = assay(data_counts_injection_readsCPM)/max_gene
### elbow mehtod to determine the number of clusters


#### by looking at the plot we can see there are ~ 6 clusters that define the data

nclus = 8


set.seed(123)
clusterRPMs = kmeans(assay(data_counts_injection_readsCPM),centers = nclus)
names_clus = paste0("clus",c(1:nclus))
diff_clusters= vector("list",nclus)
plot_clusters = vector("list",nclus)
names(plot_clusters) = names_clus
names(diff_clusters) = names_clus
TcConversion_clusters = vector("list",nclus)
plot_clusters_tc = vector("list",nclus)
### plotting these clusters now: 
totalPlot = c()

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/clusteringBasedOnQuantSeqReads.pdf")
set.seed(123)
k.max = 15
data = assay(data_counts_injection_readsCPM)
wss <- sapply(1:k.max, 
              function(k){kmeans(data, k, nstart=50,iter.max = 1000 ,algorithm="MacQueen")$tot.withinss})


plot(1:k.max, wss,
     type="b", pchw = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
for(i in 1:nclus){
  

  
  diff_clusters[[i]] <- data_counts_injection_readsCPM[which(clusterRPMs$cluster==i),]
  plot_clusters[[i]] = melt(assay(diff_clusters[[i]]))
  plot_clusters[[i]]$timepoint = rep(c(0,30,65,105,140,170),each=nrow(diff_clusters[[i]]))
  plot_clusters[[i]]$category  = paste0("clus",i) 
  q = ggplot(plot_clusters[[i]],aes(x=timepoint,y=value,group=X1)) + geom_line(aes(col=X1)) + scale_x_discrete(limits=c(0,30,65,105,140,170)) + theme(legend.position="none")
  q  = q+ ggtitle(paste0("clus",i," n=",nrow(assay(diff_clusters[[i]])))) + theme_bw() + xlab("time in minutes") + ylab("normalized CPM")
  print(q)
  p =  ggplot(plot_clusters[[i]],aes(x=timepoint,y=value,group=timepoint)) +geom_boxplot(outlier.shape = NA) + ggtitle(paste0("clus",i," n=",nrow(assay(diff_clusters[[i]])))) + theme_bw() + xlab("time in minutes") + ylab("normalized CPM")
  p = p + scale_x_discrete(limits=c(0,30,65,105,140,170)) 
  print(p)
  #TcConversion_clusters[[i]] = subsetByOverlaps(data_counts_injection_tcConversions,diff_clusters[[i]])
  TcConversion_clusters[[i]] = data_counts_injection_tcConversions[elementMetadata(data_counts_injection_tcConversions)[,3] %in% elementMetadata(diff_clusters[[i]])[,3]]
  plot_clusters_tc[[i]] = melt(assay(TcConversion_clusters[[i]]))
  plot_clusters_tc[[i]]$timepoint = rep(c(0,30,65,105,140,170),each=nrow(TcConversion_clusters[[i]]))
  r =  ggplot(plot_clusters_tc[[i]],aes(x=timepoint,y=value,group=timepoint)) +geom_boxplot(outlier.shape=NA) + ggtitle(paste0("clus",i," n=",nrow(assay(TcConversion_clusters[[i]])))) + theme_bw() + xlab("time in minutes") + ylab("raw TC conversion rate")
  r = r + scale_x_discrete(limits=c(0,30,65,105,140,170)) + ylim(c(0,0.05))
  print(r)

  
    }
dev.off()


### I want to now do a GO analysis on the different clusters by RPM
options(scipen = 999)

geneList = clusterRPMs$cluster
names(geneList) = rowRanges(data_counts_injection_readsCPM)$Name
allRes= vector("list",8)
for(i in 1:nclus){
  getGenes = function (allScore) 
  {
    return(allScore == i)
  }
  
  Odata <- new("topGOdata",ontology = "BP", allGenes = geneList, geneSel = getGenes,nodeSize = 5, annot = annFUN.org, mapping = "org.Dr.eg.db",ID = "symbol") 
  
  
  resultFisher <- runTest(Odata, algorithm = "classic", statistic = "fisher")
  allRes[[i]] <- GenTable(Odata, classicFisher = resultFisher)
}





significantResults = lapply(allRes,function(x) x[which(x$classicFisher<0.05),])

ggplotPvals = vector("list",nclus)
pdf("/Volumes//groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/goTermAnalysis.pdf")
for(i in 1:nclus){
  a = cbind.data.frame(significantResults[[i]]$Term,significantResults[[i]]$classicFisher,stringsAsFactors=F)
  colnames(a) = c("term","pVal")
  a$pVal = -log10(as.numeric(a$pVal))
  a$term = factor(a$term, levels = a$term)
  ggplotPvals[[i]] = ggplot(a,aes(y=pVal,x=term,fill=pVal)) + geom_bar(stat = "identity")+ coord_flip() + ggtitle(paste0("Cluster ",i))
  print(ggplotPvals[[i]])
}

dev.off()


##### by manual inspection it looks like cluster 6 is zygotically expressed. Now looking at this cluster in deail and reclustering it




zygoticTranscripts = TcConversion_clusters[[6]]

# assay(zygoticTranscripts) = assay(zygoticTranscripts)-assay(zygoticTranscripts)[,1]
# assay(zygoticTranscripts)[which(assay(zygoticTranscripts)<0)]<-0


assay(zygoticTranscripts) = assay(zygoticTranscripts)+ 0.00001
maxTC  = apply(assay(zygoticTranscripts),1,max)
assay(zygoticTranscripts) = assay(zygoticTranscripts)/maxTC



nclus=8

cluster_tc = kmeans(assay(zygoticTranscripts),centers = nclus)
diff_tcClusters = vector("list",nclus)

plot_clusters =  vector("list",nclus)

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/ZGAgene_cluster.pdf")
set.seed(123)
k.max = 15
data = assay(zygoticTranscripts)
wss <- sapply(1:k.max, 
              function(k){kmeans(data, k, nstart=50,iter.max = 1000 ,algorithm="MacQueen")$tot.withinss})


plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

for(i in 1:nclus){
  
  diff_tcClusters[[i]] <- zygoticTranscripts[which(cluster_tc$cluster==i),]
  plot_clusters[[i]] = melt(assay(diff_tcClusters[[i]]))
  plot_clusters[[i]]$timepoint = rep(c(0,30,65,105,140,170),each=nrow(diff_tcClusters[[i]]))
  plot_clusters[[i]]$category  = paste0("clus",i) 
  q = ggplot(plot_clusters[[i]],aes(x=timepoint,y=value,group=X1)) + geom_line(aes(col=X1)) + scale_x_discrete(limits=c(0,30,65,105,140,170)) + theme(legend.position="none")
  q  = q+ ggtitle(paste0("clus",i," n=",nrow(assay(diff_tcClusters[[i]])))) + theme_bw() + xlab("time in minutes") + ylab("normalized TC conversion rate")
  print(q)
  p =  ggplot(plot_clusters[[i]],aes(x=timepoint,y=value,group=timepoint)) +geom_boxplot() + ggtitle(paste0("clus",i," n=",nrow(assay(diff_tcClusters[[i]])))) + theme_bw() + xlab("time in minutes") + ylab("normalized TC conversion rate")
  p = p + scale_x_discrete(limits=c(0,30,65,105,140,170)) 
  print(p)
  
}


### trying out the elbow method : 


dev.off()

##### looking at the GO term analusis for this subset of ZGA genes 



options(scipen = 999)

geneList = cluster_tc$cluster
names(geneList) = rowRanges(zygoticTranscripts)$Name
allRes= vector("list",8)
for(i in 1:nclus){
  getGenes = function (allScore) 
  {
    return(allScore == i)
  }
  
  Odata <- new("topGOdata",ontology = "BP", allGenes = geneList, geneSel = getGenes,nodeSize = 5, annot = annFUN.org, mapping = "org.Dr.eg.db",ID = "symbol") 
  
  
  resultFisher <- runTest(Odata, algorithm = "classic", statistic = "fisher")
  allRes[[i]] <- GenTable(Odata, classicFisher = resultFisher)
}





significantResults = lapply(allRes,function(x) x[which(x$classicFisher<0.05),])

ggplotPvals = vector("list",nclus)
pdf("/Volumes//groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/goTermAnalysis_ZGA.pdf")
for(i in 1:nclus){
  a = cbind.data.frame(significantResults[[i]]$Term,significantResults[[i]]$classicFisher,stringsAsFactors=F)
  colnames(a) = c("term","pVal")
  a$pVal = -log10(as.numeric(a$pVal))
  a$term = factor(a$term, levels = a$term)
  ggplotPvals[[i]] = ggplot(a,aes(y=pVal,x=term,fill=pVal)) + geom_bar(stat = "identity")+ coord_flip() + ggtitle(paste0("Cluster ",i))
  print(ggplotPvals[[i]])
}

dev.off()




########



nclus = 8


set.seed(123)
clusterRPMs = kmeans(assay(data_counts_injection_tcConversions),centers = nclus,iter.max = 50)
names_clus = paste0("clus",c(1:nclus))
diff_clusters= vector("list",nclus)
plot_clusters = vector("list",nclus)
names(plot_clusters) = names_clus
names(diff_clusters) = names_clus
TcConversion_clusters = vector("list",nclus)
plot_clusters_tc = vector("list",nclus)

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/clusteringTCreads.pdf")

for(i in 1:nclus){
  
  
  
  diff_clusters[[i]] <- data_counts_injection_tcConversions[which(clusterRPMs$cluster==i),]
  plot_clusters[[i]] = melt(assay(diff_clusters[[i]]))
  plot_clusters[[i]]$timepoint = rep(c(0,30,65,105,140,170),each=nrow(diff_clusters[[i]]))
  plot_clusters[[i]]$category  = paste0("clus",i) 
  q = ggplot(plot_clusters[[i]],aes(x=timepoint,y=value,group=X1)) + geom_line(aes(col=X1)) + scale_x_discrete(limits=c(0,30,65,105,140,170)) + theme(legend.position="none")
  q  = q+ ggtitle(paste0("clus",i," n=",nrow(assay(diff_clusters[[i]])))) + theme_bw() + xlab("time in minutes") + ylab("normalized CPM")
  print(q)
  p =  ggplot(plot_clusters[[i]],aes(x=timepoint,y=value,group=timepoint)) +geom_boxplot(outlier.shape = NA) + ggtitle(paste0("clus",i," n=",nrow(assay(diff_clusters[[i]])))) + theme_bw() + xlab("time in minutes") + ylab("normalized CPM")
  p = p + scale_x_discrete(limits=c(0,30,65,105,140,170)) 
  print(p)
 
}
dev.off()

### doind the go term analysis for this : 


### I want to now do a GO analysis on the different clusters by RPM
options(scipen = 999)

geneList = clusterRPMs$cluster
names(geneList) = rowRanges(data_counts_injection_readsCPM)$Name
allRes= vector("list",8)
for(i in 1:nclus){
  getGenes = function (allScore) 
  {
    return(allScore == i)
  }
  
  Odata <- new("topGOdata",ontology = "BP", allGenes = geneList, geneSel = getGenes,nodeSize = 5, annot = annFUN.org, mapping = "org.Dr.eg.db",ID = "symbol") 
  
  
  resultFisher <- runTest(Odata, algorithm = "classic", statistic = "fisher")
  allRes[[i]] <- GenTable(Odata, classicFisher = resultFisher)
}





significantResults = lapply(allRes,function(x) x[which(x$classicFisher<0.05),])

ggplotPvals = vector("list",nclus)
pdf("/Volumes//groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/goTermAnalysis_clustersTCreads.pdf")
for(i in 1:nclus){
  a = cbind.data.frame(significantResults[[i]]$Term,significantResults[[i]]$classicFisher,stringsAsFactors=F)
  colnames(a) = c("term","pVal")
  a$pVal = -log10(as.numeric(a$pVal))
  a$term = factor(a$term, levels = a$term)
  ggplotPvals[[i]] = ggplot(a,aes(y=pVal,x=term,fill=pVal)) + geom_bar(stat = "identity")+ coord_flip() + ggtitle(paste0("Cluster ",i)) + xlab("-log10(pVal)")
  print(ggplotPvals[[i]])
}

dev.off()


