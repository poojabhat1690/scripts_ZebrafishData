  library(reshape)
  library(ggplot2)
  library(SummarizedExperiment)
  library(topGO)
  library(BSgenome.Drerio.UCSC.danRer10)
  library(GenomicRanges)
  library(checkmate)
  library(dplyr)      
  
  
  
  sampleInfoPath = read.delim("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10_data/quantseq/sampleInfoIndex.txt",sep="\t",header = F,stringsAsFactors = F)
  sampleInfoPath = sampleInfoPath[11:24,]
  sampleInfoPath$countFile = paste0(sampleInfoPath$V1,"_adapterTrimmed_slamdunk_mapped_filtered_tcount.tsv")
  # filePaths = file.path("/Volumes//groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/annotation_including3pseq_ensemblCountingWindows_18102017//",sampleInfoPath$countFile)
  filePaths = file.path("/Volumes//clustertmp/pooja/mapping_dr10_12052017/adapterTrimmed/allStages/includingEnsemblCountingWindows/TCreads_2/count/",sampleInfoPath$countFile)
  
  data_counts = lapply(filePaths,function(x) read.table(x,stringsAsFactors = F,sep="\t",header = T))
  names(data_counts) = sampleInfoPath$V3
  
  ### separating into injection and control
  
  data_counts_injection = data_counts[grep("4su",names(data_counts))]
  data_counts_control = data_counts[grep("con",names(data_counts))]
  
  ### get the CPMS and tc counts of the injection
  
  RPM_injection = do.call(cbind.data.frame,lapply(data_counts_injection,function(x) x$ReadsCPM))
  names(RPM_injection) = names(data_counts_injection)
  RPM_injection = cbind.data.frame(RPM_injection$`128_4su`,RPM_injection$`30min_4su`,RPM_injection$`65min_4su`,RPM_injection$`105min_4su`,RPM_injection$`140min_4su`,RPM_injection$`170min_4su`)
  colnames(RPM_injection) = c("2","2.5","3","3.75","4.3","48")
  
  
  tcConversions_injection = do.call(cbind.data.frame,lapply(data_counts_injection,function(x) x$ConversionRate))
  names(tcConversions_injection) = names(data_counts_injection)
  tcConversions_injection = cbind.data.frame(tcConversions_injection$`128_4su`,tcConversions_injection$`30min_4su`,tcConversions_injection$`65min_4su`,tcConversions_injection$`105min_4su` ,tcConversions_injection$`140min_4su`,tcConversions_injection$`170min_4su`)
  colnames(tcConversions_injection) =  c("2","2.5","3","3.75","4.3","48")
  
  
  #### getting the cpm and the conversion rate of the control
  
  
  ### get the CPMs and TC counts of control 
  
  RPM_control = do.call(cbind.data.frame,lapply(data_counts_control,function(x) x$ReadsCPM))
  names(RPM_control) = names(data_counts_control)
  RPM_control = cbind.data.frame(RPM_control$`128_con`,RPM_control$`30min_con`,RPM_control$`65min_con`,RPM_control$`105min_con`,RPM_control$`140min_con`,RPM_control$`170min_con`)
  colnames(RPM_control) =  c("2","2.5","3","3.75","4.3","48")
  
  
  
  tcConversions_control = do.call(cbind.data.frame,lapply(data_counts_control,function(x) x$ConversionRate))
  names(tcConversions_control) = names(data_counts_control)
  tcConversions_control = cbind.data.frame(tcConversions_control$`128_con`,tcConversions_control$`30min_con`,tcConversions_control$`65min_con`,tcConversions_control$`105min_con`,tcConversions_control$`140min_con`,tcConversions_control$`170min_con`)
  colnames(tcConversions_control) = c("2","2.5","3","3.75","4.3","48")
  
  
  
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
  meanConversions$timepoints = c(2,2.5,3,3.75,4.3,4.8)
  
  
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
  TCconversionRate_complete_mean_melt$L1 = c(2,4.8,3.75,3,4.3,2.5)
  TCconversionRate_complete_mean_melt$category = "injection"
  
  TCconversionRate_control = lapply(data_counts_control,function(x) x$ConversionsOnTs/(x$Tcontent * x$ReadCount ))
  TCconversionRate_control_complete = lapply(TCconversionRate_control,function(x) x[complete.cases(x)])
  TCconversionRate_control_complete_mean = lapply(TCconversionRate_control_complete,mean)
  TCconversionRate_control_complete_mean_melt = melt(TCconversionRate_control_complete_mean)
  TCconversionRate_control_complete_mean_melt$L1 = c(2,4.8,3.75,3,4.3,2.5)
  TCconversionRate_control_complete_mean_melt$category = "control"
  totalTCconversion = rbind(TCconversionRate_complete_mean_melt,TCconversionRate_control_complete_mean_melt)
  
  pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/backGroundTCvsInjection_conversionRatePerRead.pdf",height = 7,width=10)
  p = ggplot(totalTCconversion,aes(x=L1,y=value,group=category,col=as.factor(category))) + geom_line(size=1.25) 
  p =  p + scale_color_manual(values=myColors) + ggtitle(paste0("number of transcripts = ",nrow(data_counts_control[[1]])))
  p
  dev.off()
  
  #################################### k-means clusterg... i should have a function for this to check varying CPMS. 
  
  
  
  
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
  write.table(masterTable,'/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/masterTable_allTranscripta.txt',sep="\t",quote = F,row.names = F)
  masterTable$id = paste(masterTable$Chromosome,masterTable$Start,masterTable$End,masterTable$Name,masterTable$Strand,sep="_")
  ### so this master table is the main table from which any subset based on RPM can be made
  
  ### i should write a function here : 
    
    clusterBasedOnthreshold = function(masteTable,thresholdRPM){
      
      ###### looking at the RPMs
      rpms = masterTable[ ,grep("_RPM",colnames(masterTable))]
      rpms =  apply(rpms[,c(1:6)],1, function(x) length(which(x>=thresholdRPM))) #### the number of time points that have expression > threshold
      masterTable_5RPM = masterTable[which(rpms>1),] ### getting the genes that have at neast >1 time points greater than the threshold/
      masterTable_removed = masterTable[which(rpms<=1),]
      
      ###### working with the T>C reads
      which_tc = masterTable_5RPM[,grep("_tc",colnames(masterTable_5RPM))]
      maxTC = masterTable_5RPM[,"max_tc"]
      tcFraction =  which_tc /maxTC 
      colnames(tcFraction) = paste0(colnames(tcFraction),"_fraction")
      
      masterTable_5RPM = cbind(masterTable_5RPM,tcFraction)  
      which_rpm = masterTable_5RPM[,grep("_RPM",colnames(masterTable_5RPM))][,c(1:6)]
      max_RPM = apply(which_rpm,1,max)
      rpm_fraction = which_rpm/max_RPM
      colnames(rpm_fraction) =  paste0(colnames(rpm_fraction),"_fraction")
      masterTable_5RPM = cbind(masterTable_5RPM,rpm_fraction)
      
      noTC = masterTable_5RPM[!complete.cases(masterTable_5RPM),]
      noTC_RPMs = noTC[,grep("_RPM_fraction",colnames(noTC))]
      noTC_RPMs_melt = melt(noTC_RPMs)
      destination = "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/clusteringAnalysis/"
      fileName = paste0("transcripts_0tcConversions",thresholdRPM,".pdf")
      filePath = paste0(destination,fileName)
      pdf(filePath,height = 7,width=10)
      ggplot(noTC_RPMs_melt,aes(x = variable,y=value,group=variable)) + geom_boxplot() +ggtitle(paste0("number of trasncripts=",nrow(noTC_RPMs)))
      dev.off()
      
      
      masterTable_5RPM_TCpresent = masterTable_5RPM[complete.cases(masterTable_5RPM),]
      assertSetEqual(x = nrow(masterTable_5RPM),y=(nrow(masterTable_5RPM_TCpresent)+nrow(noTC))) ### quick check 
      nclus=8
      
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
      
      diff_clusters_RPM= vector("list",nclus)
      plot_clusters_RPM = vector("list",nclus)
      names(plot_clusters_RPM) = names_clus
      names(plot_clusters_RPM) = names_clus
      TcConversion_clusters_RPM = vector("list",nclus)
      
      
      ### plotting these clusters now: 
      totalPlot = c()
      fileNames = paste0("clusteringBasedOnQuantSeqReads_",thresholdRPM,".pdf")
      plotPath = paste0(destination,fileNames)
      ### now I can cluster these transcripts based on RPM.. then I can clasify further based on T>C conversions
      
      pdf(plotPath)
      for(i in 1:nclus){
        
        TCfraction = grep("_tc",colnames(masterTable_5RPM_TCpresent))
        TCfraction = TCfraction[c(1:6)]
        diff_clusters_TC[[i]] <- masterTable_5RPM_TCpresent[which(masterTable_5RPM_TCpresent$cluster==i),]
        tmp_table = masterTable_5RPM_TCpresent[,grep("_tc",colnames(masterTable_5RPM_TCpresent))]
        tmp_table = tmp_table[,c(1:6)]
        tmp_table = tmp_table[which(masterTable_5RPM_TCpresent$cluster==i),]
        
        plot_clusters_TC[[i]] = melt(tmp_table)
        plot_clusters_TC[[i]]$timepoint = rep(c(2,2.5,3,3.75,4.3,4.8),each=nrow(diff_clusters_TC[[i]]))
        plot_clusters_TC[[i]]$category  = paste0("clus",i) 
        q = ggplot(plot_clusters_TC[[i]],aes(x=timepoint,y=value,group=variable)) + geom_boxplot() + scale_x_discrete(limits=c(2,2.5,3,3.75,4.3,4.8)) + theme(legend.position="none")
        q  = q+ ggtitle(paste0("clus",i," n=",nrow(diff_clusters_TC[[i]]))) + theme_bw() + xlab("time in minutes") + ylab(" T>C")
        print(q)
        
        ### plotting the RPMs of these corresponding clusters 
        RPMfraction = masterTable_5RPM_TCpresent[,grep("_RPM_fraction",colnames(masterTable_5RPM_TCpresent))]
        RPMfraction$name = masterTable_5RPM_TCpresent$Name
        diff_clusters[[i]] <- RPMfraction[which(clusterRPMs$cluster==i),]
        plot_clusters[[i]] = melt(diff_clusters[[i]])
        plot_clusters[[i]]$timepoint = rep(c(2,2.5,3,3.75,4.3,4.8),each=nrow(diff_clusters[[i]]))
        plot_clusters[[i]]$category  = paste0("clus",i) 
        p = ggplot(plot_clusters[[i]],aes(x=timepoint,y=value,group=variable)) + geom_boxplot() + scale_x_discrete(limits=c(2,2.5,3,3.75,4.3,4.8)) + theme(legend.position="none")
        p  = p+ ggtitle(paste0("clus",i," n=",nrow(diff_clusters[[i]]))) + theme_bw() + xlab("time in minutes") + ylab("normalized CPM")
        print(p)
        
        
        ### total RPM 
        
        RPMfraction = masterTable_5RPM_TCpresent[,grep("_RPM",colnames(masterTable_5RPM_TCpresent))]
        RPMfraction = RPMfraction[,c(1:6)]
        RPMfraction$name = masterTable_5RPM_TCpresent$Name
        diff_clusters_RPM[[i]] <- RPMfraction[which(clusterRPMs$cluster==i),]
        plot_clusters_RPM[[i]] = melt(diff_clusters_RPM[[i]])
        plot_clusters_RPM[[i]]$timepoint = rep(c(2,2.5,3,3.75,4.3,4.8),each=nrow(diff_clusters_RPM[[i]]))
        plot_clusters_RPM[[i]]$category  = paste0("clus",i) 
        p = ggplot(plot_clusters_RPM[[i]],aes(x=timepoint,y=value,group=variable)) + geom_boxplot() + scale_x_discrete(limits=c(2,2.5,3,3.75,4.3,4.8)) + theme(legend.position="none")
        p  = p+ ggtitle(paste0("clus",i," n=",nrow(diff_clusters[[i]]))) + theme_bw() + xlab("time in minutes") + ylab("normalized CPM")
        print(p)
        
        
        
      }
      dev.off()
      
      #### plotting all this together : 
      
      plot_clusters_TC_melt = do.call(rbind,plot_clusters_TC)
      plot_clusters_TC_melt$type = "TC"
      plot_clusters_RPM_melt = do.call(rbind,plot_clusters)
      plot_clusters_RPM_melt$type = "NormalizedRPM"
      plot_clusters_RPMTotal_melt = do.call(rbind,plot_clusters_RPM)
      plot_clusters_RPMTotal_melt$type = "RPM"
      fileNames = paste0("boxplots_clusteringBasedOnRPM_",thresholdRPM,".pdf")
      plotPath = paste0(destination,fileNames)
      pdf(plotPath)
      ggplot(plot_clusters_TC_melt,aes(x=timepoint,y=value,group=timepoint)) + geom_boxplot() + facet_wrap(~category) + ylim(c(0,0.02)) + ylab("T>C reads")
      ggplot(plot_clusters_RPM_melt,aes(x=timepoint,y=value,group=timepoint)) + geom_boxplot() + facet_wrap(~category) + ylab("normalized RPM")
      ggplot(plot_clusters_RPMTotal_melt,aes(x=timepoint,y=value,group=timepoint)) + geom_boxplot() + facet_wrap(~category) + ylab("RPM")
      
      dev.off()
      
      ### plotting only the medians 
      
      split_clusters = lapply(plot_clusters_TC,function(x) split(x,x$variable,T))
      split_clusters_TC = lapply(split_clusters,function(x) lapply(x,function(y) y[,2]))
      split_clusters_TC = lapply(split_clusters_TC,function(x) lapply(x,function(y) mean(y)))
      split_clusters_TC_melt = melt(split_clusters_TC)
      split_clusters_TC_melt$tome =c(2,2.5,3,3.75,4.3,4.8)
      numberOfentriesPercluster = melt(lapply(plot_clusters_TC,function(x) nrow(x)/6))
      colnames(numberOfentriesPercluster) = c("value1","L1")
      library(plyr)
      split_clusters_TC_melt = join(split_clusters_TC_melt,numberOfentriesPercluster)
      split_clusters_TC_melt$group = paste(split_clusters_TC_melt$L1,split_clusters_TC_melt$value1,sep="_")
      fileNames= paste0("meanRPMandTC_clusters",thresholdRPM,".pdf")
      pdf(fileNames)
      
      
      
      ggplot(split_clusters_TC_melt,aes(x=tome,y=value,group=group)) + geom_line(aes(col = group),size=1 ) + ylab("Mean T>C conversion rate") + xlab("time")+scale_x_continuous(breaks=c(2,2.5,3,3.75,4.3,4.8)) 
      
      split_clusters = lapply(plot_clusters,function(x) split(x,x$variable,T))
      split_clusters_rpm = lapply(split_clusters,function(x) lapply(x,function(y) y[,3]))
      split_clusters_rpm = lapply(split_clusters_rpm,function(x) lapply(x,function(y) mean(y)))
      split_clusters_rpm_melt = melt(split_clusters_rpm)
      split_clusters_rpm_melt$tome =c(2,2.5,3,3.75,4.3,4.8)
      numberOfentriesPercluster = melt(lapply(plot_clusters,function(x) nrow(x)/6))
      colnames(numberOfentriesPercluster) = c("value1","L1")
      split_clusters_rpm_melt = join(split_clusters_rpm_melt,numberOfentriesPercluster)
      split_clusters_rpm_melt$group = paste(split_clusters_rpm_melt$L1,split_clusters_rpm_melt$value1,sep="_")
      
      ggplot(split_clusters_rpm_melt,aes(x=tome,y=value,group=group)) + geom_line(aes(col = group) ,size=1) + ylab("Mean normalized RPM") + xlab("time")+scale_x_continuous(breaks=c(2,2.5,3,3.75,4.3,4.8)) 
      
      
      #### plotting absolute RPM 
      
      
      split_clusters = lapply(plot_clusters_RPM,function(x) split(x,x$variable,T))
      split_clusters_RPM = lapply(split_clusters,function(x) lapply(x,function(y) y[,3]))
      split_clusters_RPM = lapply(split_clusters_RPM,function(x) lapply(x,function(y) mean(y)))
      split_clusters_RPM_melt = melt(split_clusters_RPM)
      split_clusters_RPM_melt$tome = c(2,2.5,3,3.75,4.3,4.8)
      numberOfentriesPercluster = melt(lapply(plot_clusters_RPM,function(x) nrow(x)/6))
      colnames(numberOfentriesPercluster) = c("value1","L1")
      library(plyr)
      split_clusters_RPM_melt = join(split_clusters_RPM_melt,numberOfentriesPercluster)
      split_clusters_RPM_melt$group = paste(split_clusters_RPM_melt$L1,split_clusters_RPM_melt$value1,sep="_")
      
      
      ggplot(split_clusters_RPM_melt,aes(x=tome,y=value,group=group)) + geom_line(aes(col = group) ,size=1) + ylab("Mean normalized RPM") + xlab("time")+scale_x_continuous(breaks=c(2,2.5,3,3.75,4.3,4.8)) 
      
      dev.off()
      
      allData = do.call(rbind,diff_clusters_TC)
      
      return(allData)
      
    }
    
    
    
  allData_1 = clusterBasedOnthreshold(masteTable = masterTable,thresholdRPM = 1)  
  allData_2 = clusterBasedOnthreshold(masteTable = masterTable,thresholdRPM = 2)
  allData_3  = allDataclusterBasedOnthreshold(masteTable = masterTable,thresholdRPM = 3)
  allData_4 = clusterBasedOnthreshold(masteTable = masterTable,thresholdRPM = 4)
  allData_5 = clusterBasedOnthreshold(masteTable = masterTable,thresholdRPM = 5)
  allData_0 = clusterBasedOnthreshold(masteTable = masterTable,thresholdRPM = 0)  
### based on the mESC dataset we can consider a RPM of 2 to be sure we have goos signal to noise ratio. 
  
  ##### so I would use a ratio of rpm 2 to define high confidence clusters. 
    
  
  allData_0_sub = allData_0[,c("id","cluster")]
  allData_2_sub = allData_2[,c("id","cluster")]
  
  ####
  
  
  