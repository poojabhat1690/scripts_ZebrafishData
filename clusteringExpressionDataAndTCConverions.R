## Including Plots
library(reshape)
library(ggplot2)
library(SummarizedExperiment)

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



library(GenomicRanges)
metadata = data_counts_injection[[1]][,c(1:6)]
metadata = makeGRangesFromDataFrame(metadata,keep.extra.columns = T,ignore.strand = F,starts.in.df.are.0based = T,seqnames.field = "Chromosome",start.field = "Start",end.field = "End",strand.field = "Strand")

data_counts_injection_readsCPM  = SummarizedExperiment(assays = RPM_injection,rowData = (metadata))
data_counts_injection_tcConversions  = SummarizedExperiment(assays = tcConversions_injection,rowData = (metadata))

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
nclus = 3
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

pdf("/Volumes/groups-1/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/clusteringBasedOnQuantSeqReads.pdf")

for(i in 1:nclus){
   
  diff_clusters[[i]] <- data_counts_injection_readsCPM[which(clusterRPMs$cluster==i),]
  plot_clusters[[i]] = melt(assay(diff_clusters[[i]]))
  plot_clusters[[i]]$timepoint = rep(c(0,30,65,105,140,170),each=nrow(diff_clusters[[i]]))
  plot_clusters[[i]]$category  = paste0("clus",i) 
  q = ggplot(plot_clusters[[i]],aes(x=timepoint,y=value,group=X1)) + geom_line(aes(col=X1)) + scale_x_discrete(limits=c(0,30,65,105,140,170)) + theme(legend.position="none")
  q  = q+ ggtitle(paste0("clus",i," n=",nrow(assay(diff_clusters[[i]])))) + theme_bw() + xlab("time in minutes") + ylab("normalized CPM")
  print(q)
  p =  ggplot(plot_clusters[[i]],aes(x=timepoint,y=value,group=timepoint)) +geom_boxplot() + ggtitle(paste0("clus",i," n=",nrow(assay(diff_clusters[[i]])))) + theme_bw() + xlab("time in minutes") + ylab("normalized CPM")
  p = p + scale_x_discrete(limits=c(0,30,65,105,140,170)) 
  print(p)
  TcConversion_clusters[[i]] = subsetByOverlaps(data_counts_injection_tcConversions,diff_clusters[[i]])
  plot_clusters_tc[[i]] = melt(assay(TcConversion_clusters[[i]]))
  plot_clusters_tc[[i]]$timepoint = rep(c(0,30,65,105,140,170),each=nrow(TcConversion_clusters[[i]]))
  r =  ggplot(plot_clusters_tc[[i]],aes(x=timepoint,y=value,group=timepoint)) +geom_boxplot() + ggtitle(paste0("clus",i," n=",nrow(assay(TcConversion_clusters[[i]])))) + theme_bw() + xlab("time in minutes") + ylab("raw TC conversion rate")
  r = r + scale_x_discrete(limits=c(0,30,65,105,140,170)) 
  print(r)
  }
dev.off()

#### it looks like there are 3 clusters based on thw quantSeq expression levels - 1. maternally provided transcripts that will devrease
#### transcripts that can be maternal or zygotic - so the expression should remain constant 3- increase in zygotically expressed transcripts. 

### it makes sense that there are no changes in TC conversions for maternal + zygotic transcripts, as they are probably those that are procused and degraded at some time
### for maternal transcripts, we would not be able to catch the production. 

### now looking a bit more in detail on the different 'waves of zygotic transctiption'. this would require clustering of the cluster 1. 


zygoticTranscripts = TcConversion_clusters[[1]]

assay(zygoticTranscripts) = assay(zygoticTranscripts)-assay(zygoticTranscripts)[,1]
assay(zygoticTranscripts)[which(assay(zygoticTranscripts)<0)]<-0


assay(zygoticTranscripts) = assay(zygoticTranscripts)+ 0.00001
maxTC  = apply(assay(zygoticTranscripts),1,max)
assay(zygoticTranscripts) = assay(zygoticTranscripts)/maxTC
nclus=5
cluster_tc = kmeans(assay(zygoticTranscripts),centers = 5)
diff_tcClusters = vector("list",nclus)

plot_clusters =  vector("list",nclus)

pdf("/Volumes/groups-1/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/ZGAgene_cluster.pdf")
for(i in 1:nclus){
  
  diff_tcClusters[[i]] <- zygoticTranscripts[which(cluster_tc$cluster==i),]
  plot_clusters[[i]] = melt(assay(diff_tcClusters[[i]]))
  plot_clusters[[i]]$timepoint = rep(c(0,30,65,105,140,170),each=nrow(diff_tcClusters[[i]]))
  plot_clusters[[i]]$category  = paste0("clus",i) 
  q = ggplot(plot_clusters[[i]],aes(x=timepoint,y=value,group=X1)) + geom_line(aes(col=X1)) + scale_x_discrete(limits=c(0,30,65,105,140,170)) + theme(legend.position="none")
  q  = q+ ggtitle(paste0("clus",i," n=",nrow(assay(diff_tcClusters[[i]])))) + theme_bw() + xlab("time in minutes") + ylab("normalized CPM")
  print(q)
  p =  ggplot(plot_clusters[[i]],aes(x=timepoint,y=value,group=timepoint)) +geom_boxplot() + ggtitle(paste0("clus",i," n=",nrow(assay(diff_tcClusters[[i]])))) + theme_bw() + xlab("time in minutes") + ylab("normalized CPM")
  p = p + scale_x_discrete(limits=c(0,30,65,105,140,170)) 
  print(p)
  
}

dev.off()