## Including Plots
library(reshape)
library(ggplot2)
library(SummarizedExperiment)
library(topGO)
library(BSgenome.Drerio.UCSC.danRer10)
library(GenomicRanges)
library(checkmate)

sampleInfoPath = read.delim("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10_data/quantseq/sampleInfoIndex.txt",sep="\t",header = F,stringsAsFactors = F)
sampleInfoPath = sampleInfoPath[11:24,]
sampleInfoPath$countFile = paste0(sampleInfoPath$V1,"_adapterTrimmed_slamdunk_mapped_filtered_tcount.tsv")
filePaths = file.path("/Volumes//groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/annotation_including3pSeq/",sampleInfoPath$countFile)

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
tcConversions_injection = cbind.data.frame(tcConversions_injection$`128_4su`,tcConversions_injection$`30min_4su`,tcConversions_injection$`65min_4su`,tcConversions_injection$`105min_4su` ,tcConversions_injection$`140min_4su`,tcConversions_injection$`170min_4su`)
colnames(tcConversions_injection) = c("0min","30min","65min","105min","140min","170min")



### get the CPMs and TC counts of control 

RPM_control = do.call(cbind.data.frame,lapply(data_counts_control,function(x) x$ReadsCPM))
names(RPM_control) = names(data_counts_control)
RPM_control = cbind.data.frame(RPM_control$`128_con`,RPM_control$`30min_con`,RPM_control$`65min_con`,RPM_control$`105min_con`,RPM_control$`140min_con`,RPM_control$`170min_con`)
colnames(RPM_control) = c("0min","30min","65min","105min","140min","170min")



tcConversions_control = do.call(cbind.data.frame,lapply(data_counts_control,function(x) x$ConversionRate))
names(tcConversions_control) = names(data_counts_control)
tcConversions_control = cbind.data.frame(tcConversions_control$`128_con`,tcConversions_control$`30min_con`,tcConversions_control$`65min_con`,tcConversions_control$`105min_con`,tcConversions_control$`140min_con`,tcConversions_control$`170min_con`)
colnames(tcConversions_control) = c("0min","30min","65min","105min","140min","170min")

#### checking control for injection for : 
    ### the TC conversion rate calcualated - i.e normalized per position
    ### the T>C conversion rate per read .. i,e #T>C/(#T * number of reads)

    ### firstly checking thr T>C conversion rate (normalized per posiiton)

        injectionMeanConversion = colMeans(tcConversions_injection)
        controlMeanConversion = colMeans(tcConversions_control)
        timePoints = names(injectionMeanConversion)
        meanConversions = cbind.data.frame(injectionMeanConversion,controlMeanConversion)
        colnames(meanConversions) = c("injection","control")
        meanConversions = reshape::melt(meanConversions)
        meanConversions$time= timePoints
        meanConversions$timepoints = c(30,60,105,128,140,170)


        assertSetEqual(nrow(tcConversions_injection),nrow(tcConversions_control))
        myColors <- c("#56ddc5", "#ff3db7")
        pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/backGroundTCvsInjection.pdf",height = 7,width=10)
        p = ggplot(meanConversions,aes(x=timepoints,y=value,group=variable,col=as.factor(variable))) + geom_line(size=1.25) + theme_bw() + xlab("time (min)") + ylab("mean T>C conversions")+theme(axis.text=element_text(size=12),axis.title=element_text(size=14),legend.text=element_text(size=14))
        p = p + ggtitle(paste0("number of transcritps = ",nrow(tcConversions_injection))) + scale_color_manual(values=myColors)
        print(p)
        dev.off()

    ##### now looking at the rate of T>C conversions, i.e #T>C reads/ #Ts 

        TCconversionRate = lapply(data_counts_injection,function(x) x$ConversionsOnTs/(x$Tcontent * x$ReadCount ))
        TCconversionRate_complete = lapply(TCconversionRate,function(x) x[complete.cases(x)])
        TCconversionRate_complete_mean = lapply(TCconversionRate_complete,mean)
        TCconversionRate_complete_mean_melt = melt(TCconversionRate_complete_mean)
        TCconversionRate_complete_mean_melt$L1 = c(0,170,105,65,140,30)
        TCconversionRate_complete_mean_melt$category = "injection"
        
        TCconversionRate_control = lapply(data_counts_control,function(x) x$ConversionsOnTs/(x$Tcontent * x$ReadCount ))
        TCconversionRate_control_complete = lapply(TCconversionRate_control,function(x) x[complete.cases(x)])
        TCconversionRate_control_complete_mean = lapply(TCconversionRate_control_complete,mean)
        TCconversionRate_control_complete_mean_melt = melt(TCconversionRate_control_complete_mean)
        TCconversionRate_control_complete_mean_melt$L1 = c(65,30,140,128,170,105)
        TCconversionRate_control_complete_mean_melt$category = "control"
        totalTCconversion = rbind(TCconversionRate_complete_mean_melt,TCconversionRate_control_complete_mean_melt)

        pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/backGroundTCvsInjection_conversionRatePerRead.pdf",height = 7,width=10)
        p = ggplot(totalTCconversion,aes(x=L1,y=value,group=category,col=as.factor(category))) + geom_line(size=1.25) 
        p =  p + scale_color_manual(values=myColors) + ggtitle(paste0("number of transcripts = ",nrow(data_counts_control[[1]])))
        p
        dev.off()
        
        
##################################################################################################################################################
        

################## clustering analysis ... this involves the following steps ... 
                    ##### removing background from T>C conversion reads
                    ##### write the foloowing data frame with the following 
                        ### chr
                        ### UTRstart
                        ### UTRend
                        ### genename
                        ### score
                        ### strand
                        ### ensembl gene id
                        ### RPM at each stage
                        ### TC conversion rate at each stage
                    ##### setting a cutoff on the RPM (mean RPM)
                    ##### clustering T>C conversion patterns. 
        
##################################################################################################################################################
        
        
  
### remove the background

      tcConversions_injection  = as.matrix(tcConversions_injection) - as.matrix( tcConversions_control)
      tcConversions_injection[which(tcConversions_injection<0)]<-0
      tcConversions_injection = as.data.frame(tcConversions_injection)
      tcConversions_injection$mean = apply(tcConversions_injection,1,mean)
      tcConversions_injection$max = apply(tcConversions_injection,1,max)
      RPM_control$mean = apply(RPM_control,1,mean)
    
        ## creating a master Table with T>C conversions and RPM (using only the control samples for this) at each stage
      
      metadata = data_counts_injection[[1]][,c(1:6)]
      colnames(tcConversions_injection) = paste0(colnames(tcConversions_injection),"_tc")
      colnames(RPM_control) = paste0(colnames(RPM_control),"_RPM")
      masterTable = cbind(metadata,tcConversions_injection,RPM_control)

## subset based on RPM ..mean RPM > 5
      library(dplyr)      
      #masterTable_5RPM = masterTable %>% filter(mean_RPM > 1)
      masterTable_5RPM = masterTable
####### now clustering the T>C conversions... 
        ### to do this I have to normalize T>C conversion rates to 1 (i.e max should be 1)
#######
      which_tc = masterTable_5RPM[,grep("min_tc",colnames(masterTable_5RPM))]
      maxTC = masterTable_5RPM[,"max_tc"]
      tcFraction =  which_tc /maxTC
      colnames(tcFraction) = paste0(colnames(tcFraction),"_fraction")
      masterTable_5RPM = cbind(masterTable_5RPM,tcFraction)  

### also normalizing the RPM to 1 
      which_rpm = masterTable_5RPM[,grep("min_RPM",colnames(masterTable_5RPM))]
      max_RPM = apply(which_rpm,1,max)
      rpm_fraction = which_rpm/max_RPM
      colnames(rpm_fraction) =  paste0(colnames(rpm_fraction),"_fraction")
      masterTable_5RPM = cbind(masterTable_5RPM,rpm_fraction)
      
      #### now I have the normalized T>C conversions and I can cluster this. there are however some transcripts that do not have any T>C conversions
            ### these could these be :
              ## maternal transcripts that do not incorportate any T>C conversions?
              ## maternal + zygotic transcripts that for some reason do not incorporate T>C conversions?
              ## zygotic transcripts that do not incorportate T>C conversions at all?
              ## could be highly stable maternal transcripts that are getting degraded very slowly?
            ### i can test this by looking at the expression profiles of these counting windows... 
      
      noTC = masterTable_5RPM[!complete.cases(masterTable_5RPM),]
      noTC_RPMs = noTC[,grep("_RPM_fraction",colnames(noTC))]
      noTC_RPMs_melt = melt(noTC_RPMs)
      
      pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/transcripts_0tcConversions.pdf",height = 7,width=10)
      ggplot(noTC_RPMs_melt,aes(x = variable,y=value,group=variable)) + geom_boxplot() +ggtitle(paste0("number of trasncripts=",nrow(noTC_RPMs)))
      dev.off()
      ### it looks like these are not maternal (or are very stable maternal transcripts) or are maternal + zygotic 
   
            #### looking at the GO terms associated with some of these example seems like they are structural components... 
      
      write.table(noTC,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/clusteringTCreads/noTCconversion.txt",sep="\t",quote = F,row.names = F)
    
      
    ### i will remove these 155 transcritps from the further analysis "
  
      
      
      masterTable_5RPM_TCpresent = masterTable_5RPM[complete.cases(masterTable_5RPM),]
      assertSetEqual(x = nrow(masterTable_5RPM),y=(nrow(masterTable_5RPM_TCpresent)+nrow(noTC))) ### quick check 
      nclus=6
      set.seed(123)
      clusterRPMs = kmeans(masterTable_5RPM_TCpresent[,grep("_RPM_fraction",colnames(masterTable_5RPM_TCpresent))],centers = nclus)
      masterTable_5RPM_TCpresent$cluster = clusterRPMs$cluster
      names_clus = paste0("clus",c(1:nclus))
      diff_clusters= vector("list",nclus)
      plot_clusters = vector("list",nclus)
      names(plot_clusters) = names_clus
      names(diff_clusters) = names_clus
      TcConversion_clusters = vector("list",nclus)
      plot_clusters_tc = vector("list",nclus)
      
      diff_clusters_TC = vector("list",nclus)
      plot_clusters_TC = vector("list",nclus)
      names(plot_clusters_TC) = names_clus
      names(diff_clusters_TC) = names_clus
      TcConversion_clusters_TC = vector("list",nclus)
      
      ### plotting these clusters now: 
      totalPlot = c()
    ### now I can cluster these transcripts based on RPM.. then I can clasify further based on T>C conversions
      
      pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/clusteringBasedOnQuantSeqReads.pdf")
      for(i in 1:nclus){
        
        TCfraction = grep("min_tc",colnames(masterTable_5RPM_TCpresent))
        TCfraction = TCfraction[c(1:6)]
        diff_clusters_TC[[i]] <- masterTable_5RPM_TCpresent[which(masterTable_5RPM_TCpresent$cluster==i),]
        tmp_table = masterTable_5RPM_TCpresent[,grep("min_tc",colnames(masterTable_5RPM_TCpresent))]
        tmp_table = tmp_table[,c(1:6)]
        tmp_table = tmp_table[which(masterTable_5RPM_TCpresent$cluster==i),]
        
        plot_clusters_TC[[i]] = melt(tmp_table)
        plot_clusters_TC[[i]]$timepoint = rep(c(0,30,65,105,140,170),each=nrow(diff_clusters_TC[[i]]))
        plot_clusters_TC[[i]]$category  = paste0("clus",i) 
        q = ggplot(plot_clusters_TC[[i]],aes(x=timepoint,y=value,group=variable)) + geom_boxplot() + scale_x_discrete(limits=c(0,30,65,105,140,170)) + theme(legend.position="none")
        q  = q+ ggtitle(paste0("clus",i," n=",nrow(diff_clusters_TC[[i]]))) + theme_bw() + xlab("time in minutes") + ylab(" T>C")
        print(q)
        
        ### plotting the RPMs of these corresponding clusters 
        RPMfraction = masterTable_5RPM_TCpresent[,grep("_RPM_fraction",colnames(masterTable_5RPM_TCpresent))]
        RPMfraction$name = masterTable_5RPM_TCpresent$Name
        diff_clusters[[i]] <- RPMfraction[which(clusterRPMs$cluster==i),]
        plot_clusters[[i]] = melt(diff_clusters[[i]])
        plot_clusters[[i]]$timepoint = rep(c(0,30,65,105,140,170),each=nrow(diff_clusters[[i]]))
        plot_clusters[[i]]$category  = paste0("clus",i) 
        p = ggplot(plot_clusters[[i]],aes(x=timepoint,y=value,group=variable)) + geom_boxplot() + scale_x_discrete(limits=c(0,30,65,105,140,170)) + theme(legend.position="none")
        p  = p+ ggtitle(paste0("clus",i," n=",nrow(diff_clusters[[i]]))) + theme_bw() + xlab("time in minutes") + ylab("normalized CPM")
        print(p)
      }
      dev.off()
      
    
      
      
            #### plotting all this together : 
     
      plot_clusters_TC_melt = do.call(rbind,plot_clusters_TC)
      plot_clusters_TC_melt$type = "TC"
      plot_clusters_RPM_melt = do.call(rbind,plot_clusters)
      plot_clusters_RPM_melt$type = "NormalizedRPM"
      
      
      pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/boxplots_clusteringBasedOnRPM.pdf")
      ggplot(plot_clusters_TC_melt,aes(x=timepoint,y=value,group=timepoint)) + geom_boxplot() + facet_wrap(~category) + ylim(c(0,0.02)) + ylab("T>C reads")
      ggplot(plot_clusters_RPM_melt,aes(x=timepoint,y=value,group=timepoint)) + geom_boxplot() + facet_wrap(~category) + ylab("normalized RPM")
      dev.off()
      
      ### clusters 1,3,4,5 look to be maternl + zygotic and cluster 6 looks to be zygotic and cluster 2 to be maternal, but cluster 2 still has some T>C conversions!!! why?
      diff_clusters_TC[[1]]$categoty = "MZ"  
      diff_clusters_TC[[3]]$categoty = "MZ"  
      diff_clusters_TC[[4]]$categoty = "MZ"  
      diff_clusters_TC[[5]]$categoty = "MZ"  
      diff_clusters_TC[[2]]$categoty = "M"  
      diff_clusters_TC[[6]]$categoty = "Z"  
      
      allData = do.call(rbind,diff_clusters_TC)
      write.table(allData, "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/clusteringTCreads/categories_slamSeq.txt",sep="\t",quote = F,row.names = F)
      # pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/boxplots_subclusteringZygoticTranscripts.pdf")
      # clusters = kmeans(x = diff_clusters_TC[[6]][,c(7:12)],3 )
      # test_zygoticCluster = diff_clusters_TC[[6]]
      # test_zygoticCluster$cluster = clusters$cluster
      # 
      # a = split(test_zygoticCluster,test_zygoticCluster$cluster,T)
      # a = lapply(a,function(x) x[,c(7:12)])
      # a_melt = melt(a)
      # timepointSplit = strsplit(as.character(a_melt$variable),"min",T)
      # timepointSplit = unlist(lapply(timepointSplit,function(x) x[1]))
      # timepointSplit = as.numeric(timepointSplit)
      # a_melt$timepoint = timepointSplit
      # ggplot(a_melt,aes(x=timepoint,y=value,group=variable)) + geom_boxplot() + facet_wrap(~L1) 
      # 
      # ### checking the RPMs of the same
      # test_zygoticCluster = diff_clusters[[6]]
      # test_zygoticCluster$cluster = clusters$cluster
      # 
      # a = split(test_zygoticCluster,test_zygoticCluster$cluster,T)
      # a = lapply(a,function(x) x[,c(1:6)])
      # a_melt = melt(a)
      # timepointSplit = strsplit(as.character(a_melt$variable),"min",T)
      # timepointSplit = unlist(lapply(timepointSplit,function(x) x[1]))
      # timepointSplit = as.numeric(timepointSplit)
      # a_melt$timepoint = timepointSplit
      # ggplot(a_melt,aes(x=timepoint,y=value,group=variable)) + geom_boxplot() + facet_wrap(~L1)
      # 
      # dev.off()
      # 
      # 
      
##### now I want to compare this to Lee et. al 2014 that did this classification based on intron signal
       
      leeEtal = read.delim("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/externalDatat/zygoticAndMaternalTranscripts_differentStudies/LeeEtAl_2014.txt",sep="\t")
      table(leeEtal$Maternal_contr)
      names_zygotic = allData %>% filter(categoty == "Z")
      intersect(leeEtal$Symbol,names_zygotic$name)
      setdiff(leeEtal$Symbol,names_zygotic$name)
   
      
      
      
      
      
      
      
      
      
      
      
      
      
      
         
      
          ##### creating summarized experiments for the RPM and the T>C conversion matrices :
       
  

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

test_data = cbind.data.frame(assay(data_counts_injection_readsCPM),assay(data_counts_injection_tcConversions))

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
names(geneList) = unique(rowRanges(zygoticTranscripts)$Name)
allRes= vector("list",8)
for(i in 1:nclus){
  getGenes = function (allScore) 
  {
    return(unique(allScore == i))
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
clusters_rpm = vector("list",nclus)
rpmClusters_melt = vector("list",nclus)
pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/clusteringTCreads.pdf")

for(i in 1:nclus){
  
  clusters_rpm[[i]] = data_counts_injection_readsCPM[which(clusterRPMs$cluster==i),]
  rpmClusters_melt[[i]]= melt(assay(clusters_rpm[[i]]))
  rpmClusters_melt[[i]]$timepoint = rep(c(0,30,65,105,140,170),each=nrow(clusters_rpm[[i]]))
  rpmClusters_melt[[i]]$category  = paste0("clus",i) 
  r =  ggplot(rpmClusters_melt[[i]],aes(x=timepoint,y=value,group=timepoint)) +geom_boxplot(outlier.shape = NA)  + ggtitle(paste0("clus",i," n=",nrow((rpmClusters_melt[[i]]))))+ theme_bw() + xlab("time in minutes") + ylab("normalized CPM")
  r = r + scale_x_discrete(limits=c(0,30,65,105,140,170)) 
  print(r)
  
  diff_clusters[[i]] <- data_counts_injection_tcConversions[which(clusterRPMs$cluster==i),]
  plot_clusters[[i]] = melt(assay(diff_clusters[[i]]))
  plot_clusters[[i]]$timepoint = rep(c(0,30,65,105,140,170),each=nrow(diff_clusters[[i]]))
  plot_clusters[[i]]$category  = paste0("clus",i) 
  q = ggplot(plot_clusters[[i]],aes(x=timepoint,y=value,group=X1)) + geom_line(aes(col=X1)) + scale_x_discrete(limits=c(0,30,65,105,140,170)) + theme(legend.position="none")
  q  = q+ ggtitle(paste0("clus",i," n=",nrow(assay(diff_clusters[[i]])))) + theme_bw() + xlab("time in minutes") + ylab("T>C conversion rate") + ylim(c(0,0.25))
  print(q)
  p =  ggplot(plot_clusters[[i]],aes(x=timepoint,y=value,group=timepoint)) +geom_boxplot(outlier.shape = NA) + ggtitle(paste0("clus",i," n=",nrow(assay(diff_clusters[[i]])))) + theme_bw() + xlab("time in minutes") + ylab("T>C conversion")
  p = p + scale_x_discrete(limits=c(0,30,65,105,140,170)) + ylim(c(0,0.01))
  print(p)
 
}
dev.off()

#### plottong the mean of the expression 



### doind the go term analysis for this : 

diff_clusters_numbers = unlist(lapply(diff_clusters,function(x) nrow(assay(x))))
names(diff_clusters) = paste(names(diff_clusters),"(n=",diff_clusters_numbers,")")
mean_clusters = lapply(diff_clusters,function(x) colMeans(assay(x)))
mean_clusters_melt = melt(mean_clusters)
mean_clusters_melt$time = c(0,30,65,105,140,170)

meanClusters_rpm = lapply(clusters_rpm,function(x) colMeans(assay(x)))
meanClusters_rpm_melt = melt(meanClusters_rpm)
meanClusters_rpm_melt$time =  c(0,30,65,105,140,170)

colMeans(assay(clusters_rpm[[1]]))
pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/Clusters_TC_meanConversions.pdf")
ggplot(mean_clusters_melt,aes(x=time,y=value,group=L1,col=L1))+geom_line(size=1) + xlab("time") + ylab("mean T>C conversions")
q = ggplot(meanClusters_rpm_melt,aes(x=time,y=value,group=L1,col=L1))+geom_line(size=1) + xlab("time") + ylab("normalized mean conversions")
print(q)
dev.off()
### I want to now do a GO analysis on the different clusters by RPM
options(scipen = 999)

geneList = clusterRPMs$cluster
names(geneList) = rowRanges(data_counts_injection_readsCPM)$Name
allRes= vector("list",nclus)
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


### the GO term analysis reveals some interesting patterns of gene activation.



