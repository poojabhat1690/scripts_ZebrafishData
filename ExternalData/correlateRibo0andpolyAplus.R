### checking how well polyA+ and ribo zero libraries correlate with each other. 
### for detailed description of the results look on one note for the note : polyA+ Ribo0 QuantSeq comparison in Zebrafish project notes

####installation

#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

###

### loading libraries 
library(corrplot)
library(plyr)
library(reshape)
library(ggplot2)
library(DESeq2)
library(biomaRt)
#############################

## Data details ####
#### reading in the gene counts that Thomas has created at some time from Andi. She has performed polyA+ RNAseq and Ribo0 seq from the same smaples. 
#### I will use this for the comparison....
#t1 = 4-cell embryos (~1 hour post fertilization (hpf))
#t2 = 64-128 cell (~2.5 hpf)
#t3 = 256-512 cell (~3.5 hpf)
#t4 = oblong/sphere (~4.5 hpf)




####################################################################### 
#### analysis - correlating polyA+ and riboCop data. 
  ### corrlate all genes - pearson and spearman
  ### correlate highly expressed genes - pearson and spearman
####################################################################### 

## loading the data. 

polyAAndRibo0 = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/externalDatat/fromAndi_ribo0VsPolyA/tpm_genes.txt",header=T)

### renaming samples based on the column names and the description Andi sent (see above data description)

sampleNames = strsplit(colnames(polyAAndRibo0),"_",T)
sampleNames = unlist(lapply(sampleNames,function(x) paste(x[2],x[3],x[4],sep = "_")))
sampleNames = sampleNames[-1] ### gene id was also a column name

times = rep(c("4cell","64-128cell","256-512cell","oblong/sphere"),each=3) ### based on data description
times = paste(times,sampleNames,sep="_")

rownames(polyAAndRibo0) = polyAAndRibo0$geneID
polyAAndRibo0$geneID = NULL
colnames(polyAAndRibo0) = times

##### do the replicates correlate amongst each other?
  

  ### with all genes 
  

  correlation_allGenes_pearson = cor(polyAAndRibo0)
  correlation_allGenes_spearman = cor(polyAAndRibo0,method = "spearman")
  
  jpeg("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/comparisonOfPolyARiboZeroQuantSeq/correlation_allGenes_pearson_polyA_ribo0.jpeg")
  corrplot(correlation_allGenes_pearson)
  dev.off()
  
  jpeg("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/comparisonOfPolyARiboZeroQuantSeq/correlation_allGenes_spearman_polyA_ribo0.jpeg")
  corrplot(correlation_allGenes_spearman)
  dev.off()
  
  ### from this it looks like one RiboCop sample (at t1) and one polyA + samples (at t1) do not correlate well with other replicates at the same time point
  ### when i look at why this is : from the mapping and counting data that thomas looked at, it seems that : 
    ### sample 
      ### 4cell_t1_polyA_49983 - has 2368295 reads mapping to exon in comparison to 7188554 and 10561980 reads in the other replicates. (not sequenced at enough depth- but still comparable), look at the fastqC file to see what is going 
          ## this sample also correlates well with the ribo0 seq, maybe the polyA selection did not work so well?
      ### 4cell_t1_RiboCop_50011 - has only 477322 reads mapping to exons in comparison to 3507766 and 4242330 reads of the same replicates - has enough reads aligned but 97% of these align to introns.
  
  ### the spearman correlation is good, this means that there is no proportional increase/decrease of transcripts. But between replicates I would think that the peasron correlation should be good. 
  
  
  
  ### looking if this trend also holds good for hightly expressed genes. Since we see that there are two samples that are outliers, accounting for this: 
      ### and eliminating them from the selection for highly expressed genes. 
  
    #### looking for genes that have at lease 2 replicates (any time point, excluding the outlier samples) expressed greater than 5TPM
  
  ### excluding the outlier samples
  
  polyAAndRibo0_subset = polyAAndRibo0
  polyAAndRibo0_subset$`4cell_t1_polyA_49983` = NULL 
  polyAAndRibo0_subset$`4cell_t1_RiboCop_50011` = NULL
  
  maxTpm = apply(polyAAndRibo0_subset,1,function(x) length(which(x>5)))
  polyAAndRibo0_5TPM = polyAAndRibo0[which(maxTpm>1),]
  
  correlation_highExpression_pearson = cor(polyAAndRibo0_5TPM)
  correlation_highExpression_spearman = cor(polyAAndRibo0_5TPM,method = "spearman")
  
  
  jpeg("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/comparisonOfPolyARiboZeroQuantSeq/correlation_5TPMcutoff_pearson_polyA_ribo0.jpeg")
  corrplot(correlation_highExpression_pearson)
  dev.off()
  
  jpeg("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/comparisonOfPolyARiboZeroQuantSeq/correlation_5TPMcutoff_spearman_polyA_ribo0.jpeg")
  corrplot(correlation_highExpression_spearman)
  dev.off()
  
  #### looking at fold changes between these  plots of these two libraries (only highly expressed genes) to get a better idea:
  
  
  #### inspecting sample 4cell_t1_polyA_49983
  
  foldChange_49983_49979 = (polyAAndRibo0_5TPM$`4cell_t1_polyA_49983`)/(polyAAndRibo0_5TPM$`4cell_t1_polyA_49979`)
  foldChange_49983_49987 = (polyAAndRibo0_5TPM$`4cell_t1_polyA_49983`)/(polyAAndRibo0_5TPM$`4cell_t1_polyA_49987`)
  foldChange_49979_49987 = (polyAAndRibo0_5TPM$`4cell_t1_polyA_49987`)/(polyAAndRibo0_5TPM$`4cell_t1_polyA_49979`)
  
  foldChange_49983_49979 = foldChange_49983_49979[-which(foldChange_49983_49979 == "NaN"| foldChange_49983_49979 == "Inf")]
  foldChange_49983_49987 = foldChange_49983_49987[-which(foldChange_49983_49987 == "NaN"| foldChange_49983_49987 == "Inf")]
  foldChange_49979_49987 = foldChange_49979_49987[-which(foldChange_49979_49987=="NaN"| foldChange_49979_49987 == "Inf")]
  
  
  allFoldChanges_polyA = list(foldChange_49983_49979,foldChange_49983_49987,foldChange_49979_49987)
  names(allFoldChanges_polyA) = c("foldChange_49983_49979","foldChange_49983_49987","foldChange_49979_49987")
  
  allFoldChanges_polyA_melt = melt(allFoldChanges_polyA)
  jpeg("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/comparisonOfPolyARiboZeroQuantSeq/densityOfFoldChanges_49983.jpg")
  ggplot(allFoldChanges_polyA_melt,aes(log2(value),group=L1 )) + geom_density(aes(col=L1),size=1) + xlab("log2(FoldChange)")
  dev.off()
  
  ### the sample 4cell_t1_polyA_49983 just reports lower values in comparison to the other samples 
  
  
  #### inspecting sample 4cell_t1_RiboCop_50011
  foldChange_50011_50007 = polyAAndRibo0_5TPM$`4cell_t1_RiboCop_50011`/polyAAndRibo0_5TPM$`4cell_t1_RiboCop_50007`
  foldChange_50011_50015 = polyAAndRibo0_5TPM$`4cell_t1_RiboCop_50011`/polyAAndRibo0_5TPM$`4cell_t1_RiboCop_50015`
  foldChange_50007_50015 = polyAAndRibo0_5TPM$`4cell_t1_RiboCop_50015`/polyAAndRibo0_5TPM$`4cell_t1_RiboCop_50007`
  
  foldChange_50011_50007 = foldChange_50011_50007[-which(foldChange_50011_50007 == "NaN" | foldChange_50011_50007 == "Inf")]
  foldChange_50011_50015 = foldChange_50011_50015[-which(foldChange_50011_50015 == "NaN"| foldChange_50011_50015 == "Inf")]
  foldChange_50007_50015 = foldChange_50007_50015[-which(foldChange_50007_50015 == "NaN"| foldChange_50007_50015 == "Inf")]
  
  allFoldChanges_riboCop = list(foldChange_50011_50007,foldChange_50011_50015,foldChange_50007_50015)
  names(allFoldChanges_riboCop) = c("foldChange_50011_50007","foldChange_50011_50015","foldChange_50007_50015")
  
  allFoldChanges_riboCop_melt = melt(allFoldChanges_riboCop)
  
  jpeg("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/comparisonOfPolyARiboZeroQuantSeq/densityOfFoldChanges_50011.jpg")
  ggplot(allFoldChanges_riboCop_melt,aes(log2(value),group=L1 )) + geom_density(aes(col=L1),size=1)+ xlab("log2(FoldChange)")
  dev.off()
  ### the difference is more pronounced in this case and the sample 4cell_t1_RiboCop_50011 reports much lower values than the other samples
  
  
  
### the replicates of polyA+ RNAseq and ribo0RNA seq correlate well amongst each other, but there are two outliers that do not correlate well. These are thechnically not very good samples
  
  
  
### but it also seems that there is no correlation between the same time point between polyA+ rnaSeq and Ribo0 based RNAseq. To test this I will do DESeq2 analysis 
  ## always considering polyA samples as untreated
  
  ### reading in the count data 
  
  
  countMatrix = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/externalDatat/fromAndi_ribo0VsPolyA/rawcounts.txt")
  
  runDeSeq2= function(countData,timePoint){
    countMatrix = countData 
    condition=c("polyA","polyA","polyA","ribo0","ribo0","ribo0")
    type=c("single-read","single-read","single-read","single-read","single-read","single-read")  
    coldata = cbind.data.frame(condition,type)
    row.names(coldata) = colnames(countMatrix)
    all(rownames(coldata) %in% colnames(countMatrix))
    
    dds <- DESeqDataSetFromMatrix(countData = countMatrix,colData = coldata,design = ~ condition)
    dds <- dds[ rowSums(counts(dds)) > 1, ]
    dds <- DESeq(dds)
    res <- results(dds)
    fileName = paste0("differnetialExpression_",timePoint) 
    destination= "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/comparisonOfPolyARiboZeroQuantSeq/"
    jpeg(paste0(destination,fileName,".jpg"))
    plotMA(res, ylim=c(-6,6))
    dev.off()
    res05 <- results(dds, alpha=0.05)
    res05_df = as.data.frame(res05)
    return(res05_df)
  }
  
  
  
  
#   condition=c("polyA","polyA","polyA","ribo0","ribo0","ribo0")
#   type=c("single-read","single-read","single-read","single-read","single-read","single-read")  
# coldata = cbind.data.frame(condition,type)
  # row.names(coldata) = colnames(countMatrix[1:6])
  # all(rownames(coldata) %in% colnames(countMatrix))
  # 
  # dds <- DESeqDataSetFromMatrix(countData = countMatrix[,c(1:6)],colData = coldata,design = ~ condition)
  # 
  # ### removing low count data  
  # dds <- dds[ rowSums(counts(dds)) > 1, ]
  # dds <- DESeq(dds)
  # res <- results(dds)
  #   
  
  
  DeSeq_t1 =  runDeSeq2(countData = countMatrix[,c(1:6)],timePoint = "t1")
  DeSeq_t2 = runDeSeq2(countData = countMatrix[,c(7:12)],timePoint = "t2")
  DeSeq_t3 = runDeSeq2(countData = countMatrix[,c(13:18)],timePoint = "t3")
  DeSeq_t4 = runDeSeq2(countData = countMatrix[,c(19:24)],timePoint = "t4")
 
  
  
  ####### now lookin at this data in a bit detail ... 
  
  
    ### first I want to see how many genes are upregulated and how many downregulated at each time point. This is a comparison between ribo0 and polyA+.
    ### I want to find out the following with this analysis : 
      ### 1. What are the genes that are upregulated (expressed higher in ribo0)
      ### 2. what are the genes that are downregulated (expressed higher in polyA+)
      ### 3. Are the genes that are downregulated common in all samples --> do these have high genomic A content?
  

  getUpregulatedAndDownregulated  = function(dataFrame){
    deGenes = dataFrame[which(dataFrame$padj<0.05),]
    upRegulated_deGenes = deGenes[which(deGenes$log2FoldChange>0),]
    downRegulated_deGenes = deGenes[which(deGenes$log2FoldChange<0),]
    upregulatedAndDownregulated = list(upRegulated_deGenes,downRegulated_deGenes)
    names(upregulatedAndDownregulated) = c("upRegulated","downRegulated")
    return(upregulatedAndDownregulated)
  }
    

  
  t1_upreg_downreg = getUpregulatedAndDownregulated(dataFrame = DeSeq_t1)  
  t2_upreg_downreg = getUpregulatedAndDownregulated(dataFrame = DeSeq_t2)  
  t3_upreg_downreg = getUpregulatedAndDownregulated(dataFrame = DeSeq_t3)  
  t4_upreg_downreg = getUpregulatedAndDownregulated(dataFrame = DeSeq_t4)  
  
  Num_t1  = lapply(t1_upreg_downreg ,nrow)
  Num_t2  = lapply(t2_upreg_downreg ,nrow)
  Num_t3  = lapply(t3_upreg_downreg ,nrow)
  Num_t4  = lapply(t4_upreg_downreg ,nrow)  
  allTimes_num =  list(Num_t1,Num_t2,Num_t3,Num_t4)
  names(allTimes_num) = c("t1","t2","t3","t4")
  allTimes_num_melt = melt(allTimes_num)
  jpeg("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/comparisonOfPolyARiboZeroQuantSeq/numberOfUpAndDownRegulatedGenes_timepoint.jpg")
  ggplot(allTimes_num_melt,aes(x=L1,y=value,fill=L2))+geom_bar(stat = "identity") + ylab("Number of genes") + xlab("timepoint")
  dev.off()
  ### so the number of genes that are differnetially expressed actually decreases at later time points. So at least partially the polyA tail lengths may make a difference.
  
  ### before we go into looking at the polyA tails of transcripts, I want to test the following : 
  
    ## The transcripts that are not picked up by polyA+ RNAseq have short polyA tails? - up regulated transcripts
    ## What are the transcripts that are upregulated in polyA+ RNA and not in RiboCop?- do these have long genomic A stretches that get primed? - downregulated transcripts

  
  
  ### to check if up regulated transcripts are enriched for short poly A tails, I will compare polyA tails of the up regulated transcripts with polyA tail lengths of all transcripts. 
  
    ## I have the polyA tail lengths from Subtenley et.al 
  
  samplesPolyAtail = list.files("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/externalDatat/polyAtailLengths_SutbtenlyEtal",pattern = "*.txt")
  samplesPolyAtail = samplesPolyAtail[-grep(".gz",samplesPolyAtail)]
  samplesPolyAtail = samplesPolyAtail[grep("mock",samplesPolyAtail)]
  path_samplesPolyAtail = paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/externalDatat/polyAtailLengths_SutbtenlyEtal/",samplesPolyAtail)
  
  ## make sure that the first two lines are deleted from these files
  
  polyAtailData = lapply(path_samplesPolyAtail,function(x) read.delim(x))
  names(polyAtailData) = c("2hpf","4hpf","6hpf")
  polyAtailData_meanTailLength = lapply(polyAtailData,function(x) cbind.data.frame(x$Mean.TL,x$Gene.name))
  colnames(polyAtailData_meanTailLength$`2hpf`) = c("meanTail_2hpf","ensembl_transcript_id")
  colnames(polyAtailData_meanTailLength$`4hpf`) = c("meanTail_4hpf","ensembl_transcript_id")
  colnames(polyAtailData_meanTailLength$`6hpf`) = c("meanTail_6hpf","ensembl_transcript_id")
  
  ###### now I have the tail length, first to check this data if there is actually an increase in the polyA tail length through MZT
  onlyPolyAtail = lapply(polyAtailData,function(x) x$Mean.TL)
  onlyPolyAtail_melt = melt(onlyPolyAtail)
  jpeg("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/comparisonOfPolyARiboZeroQuantSeq/polyAtailLengths_Subtenley.jpg")
  ggplot(onlyPolyAtail_melt,aes(x=L1,y=value)) + geom_boxplot() + xlab("timepoint") + ylab("polyA tail length (nt)") 
  dev.off()
  
  
  join_2hpf_4hpf = join(polyAtailData_meanTailLength$`2hpf`,polyAtailData_meanTailLength$`4hpf`)
  join_2hpf_4hpf_6hpf = join(join_2hpf_4hpf,polyAtailData_meanTailLength$`6hpf`)
  
  ## I need to get the ensembl gene ids to compare it to the differnetially expressed genes. 
  
  ensembl = useEnsembl(biomart="ensembl")
  ensembl = useEnsembl(biomart="ensembl", dataset="drerio_gene_ensembl")
  geneIdsAndNames = getBM(attributes = c("ensembl_gene_id","external_gene_name"),filters = "external_gene_name",values = join_2hpf_4hpf_6hpf$ensembl_transcript_id,mart = ensembl)
  
  colnames(join_2hpf_4hpf_6hpf) = c("meanTail_2hpf","external_gene_name","meanTail_4hpf","meanTail_6hpf")
  join_2hpf_4hpf_6hpf_geneIds = join(geneIdsAndNames,join_2hpf_4hpf_6hpf)
  
  #### now that i have this data formatted, I want to look at the mean tail lengths of genes enriched in each of the differnetially expressed sets
  

  plotTailLengths = function(upDown,timepoint,tailDf,time_de){
    t1_upregulated = upDown$upRegulated
    t1_upregulated$ensembl_gene_id = row.names(t1_upregulated)
    t1_upregulated = join(t1_upregulated,tailDf)
    t1_downregulated = upDown$downRegulated 
    t1_downregulated$ensembl_gene_id = row.names(t1_downregulated)
    t1_downregulated = join(t1_downregulated,tailDf)
    
    tail_downRegulated =  t1_downregulated[,grep(timepoint,colnames(t1_downregulated))]
    tail_upRegulated  = t1_upregulated[,grep(timepoint,colnames(t1_upregulated))]
    allTailLengths = tailDf[,grep(timepoint,colnames(tailDf))]
    allTailLengths = list(tail_upRegulated,tail_downRegulated,allTailLengths)
    names(allTailLengths)  = c("upregulated","downegulated","allTailLengths")
    allTailLengths_melt = melt(allTailLengths)
    fileName = paste0("tailLengthVariability_",time_de) 
    destination= "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/comparisonOfPolyARiboZeroQuantSeq/"
    jpeg(paste0(destination,fileName,".jpg"))
    p = ggplot(allTailLengths_melt,aes(y=value,x=L1,group=L1)) + geom_boxplot() + ggtitle(timepoint) + xlab("category") + ylab("tail length") + ylim(c(0,200))
    print(p)
    dev.off()
    fileName = paste0("tailLengthVariabilityDensity_",time_de) 
    
    jpeg(paste0(destination,fileName,".jpg"))
    p = ggplot(allTailLengths_melt,aes(value,group=L1)) + geom_density(aes(col=L1)) + ggtitle(timepoint) 
    print(p)
    dev.off()
    return(allTailLengths_melt)
  }
  
  tailDistributions_t1 = plotTailLengths(upDown =t1_upreg_downreg,timepoint = "2hpf",tailDf = join_2hpf_4hpf_6hpf_geneIds ,"t1")
  tailDistributions_t2 = plotTailLengths(upDown =t2_upreg_downreg,timepoint = "2hpf",tailDf = join_2hpf_4hpf_6hpf_geneIds,"t2" )
  tailDistributions_t3 = plotTailLengths(upDown =t3_upreg_downreg,timepoint = "4hpf",tailDf = join_2hpf_4hpf_6hpf_geneIds,"t3" )
  tailDistributions_t4 = plotTailLengths(upDown =t4_upreg_downreg,timepoint = "4hpf",tailDf = join_2hpf_4hpf_6hpf_geneIds ,"t4")
  
  
  
  ### just want to plot this together
  
  
  allTails = list(tailDistributions_t1,tailDistributions_t2,tailDistributions_t3,tailDistributions_t4)
  names(allTails) = c("t1","t2","t3","t4")
  allTails_melt = melt(allTails,level = "names")
  jpeg("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/comparisonOfPolyARiboZeroQuantSeq/tailLengths_allTimes.jpg")
  p = ggplot(allTails_melt,aes(value,group=L1))+ geom_density(aes(col=L1))+ facet_wrap(~Lnames) + xlab("tail length(nt)")
  print(p)
  dev.off()
  
  ##### so it does look like at all time points the upregulated genes (captured in ribo0 and not in polyA+ have shorter polyA tails)
  
  
###  on the other hand there are also genes upregulated in the polyA+ libraries in comparison to ribo0 libraries. 
first_intersect =   intersect(row.names(t1_upreg_downreg$downRegulated),row.names(t2_upreg_downreg$downRegulated))
second_intersect = intersect(first_intersect,row.names(t3_upreg_downreg$downRegulated)) 
third_intersect = intersect(second_intersect,row.names(t4_upreg_downreg$downRegulated))

downRegulatedInAll = third_intersect


  ### dont see any clear A enrichment anywhere


#### looking at quantSeq with respect to this analysis so far
### first I want to see how well the quantSeq samples correlate with each other



countDataSets_quantSeq = list.files(path = "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/analysis_ensemblAnnotation/countMatrices/",pattern = ".tsv")
countDataSets_quantSeq_paths = paste0( "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/analysis_ensemblAnnotation/countMatrices/",countDataSets_quantSeq)
countDataSets_quantSeq_data = lapply(countDataSets_quantSeq_paths,function(x) read.table(x,sep="\t",header=T))
names(countDataSets_quantSeq_data) = countDataSets_quantSeq

countDataSets_quantSeq_cpm = lapply(countDataSets_quantSeq_data,function(x) x$ReadsCPM)
countDataSets_quantSeq_cpm = do.call(cbind.data.frame,countDataSets_quantSeq_cpm)
countDataSets_quantSeq_cpm$names = countDataSets_quantSeq_data[[1]]$Name

#names(countDataSets_quantSeq_cpm) = countDataSets_quantSeq
a = split( countDataSets_quantSeq_cpm,countDataSets_quantSeq_cpm$names,T )
a_summed = lapply(a,function(x) colSums(x[1:24]))
a_summed = do.call(rbind,a_summed)
a_summed = as.data.frame(a_summed)
a_summed$ensembl_gene_id = row.names(a_summed)

samplesNames_quantSeq = read.delim("/Volumes//groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/sampleInfo.txt",header = F,stringsAsFactors = F)
colnames(a_summed) = c(samplesNames_quantSeq$V3,"GeneName")
a_summed$GeneName = NULL

jpeg("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/comparisonOfPolyARiboZeroQuantSeq/correlation_onlyQuantSeq.jpeg")
corrplot(cor(a_summed))
dev.off()


########## now I want to check how the polyA+ samples relate to the quantSeq samples : 
### cleaning up the columns (I dont want to have random samples)


relavantQuantSeq  = a_summed[,c("2cell","128_con" , "30min_con", "256cell",    "65min_con","sphere","140min_con","170min_con")]
### 2 cell = 0.5 hr
### 256 cell - 3.5 hr
### sphere - 4.5 hr
### dome - 6.5 hrs??

colnames(relavantQuantSeq) = c("0.5hpf",   "2.5hpf", "3.5hpf",    "3.5hpf" ,"4hpf","4.5hpf","5hpf","6hpf")
colnames(relavantQuantSeq) = paste0(colnames(relavantQuantSeq),"_quantSeq")

# geneNames_ensembl_gene_id = getBM(attributes = c("external_gene_name","ensembl_gene_id"),filters = "external_gene_name",values = row.names(relavantQuantSeq),mart = ensembl)
# relavantQuantSeq$external _gene_name = row.names(relavantQuantSeq)
# 
# relavantQuantSeq = join(relavantQuantSeq,geneNames_ensembl_gene_id)
# relavantQuantSeq = relavantQuantSeq[complete.cases(relavantQuantSeq),]
# row.names(relavantQuantSeq) = relavantQuantSeq$ensembl_gene_id

### get the mean of time points of polyA + and riboCop data 

###### not considering the two outlier smaples in calculating the means 

t1_mean_polyA = (polyAAndRibo0_5TPM$`4cell_t1_polyA_49979`+polyAAndRibo0_5TPM$`4cell_t1_polyA_49987`)/2
t2_mean_polyA =( polyAAndRibo0_5TPM$`64-128cell_t2_polyA_49980`+polyAAndRibo0_5TPM$`64-128cell_t2_polyA_49984`+polyAAndRibo0_5TPM$`64-128cell_t2_polyA_49988`  )/3
t3_mean_polyA = (polyAAndRibo0_5TPM$`256-512cell_t3_polyA_49981`+polyAAndRibo0_5TPM$`256-512cell_t3_polyA_49985`+polyAAndRibo0_5TPM$`256-512cell_t3_polyA_49989`)/3
t4_mean_polyA = (polyAAndRibo0_5TPM$`oblong/sphere_t4_polyA_49982`+polyAAndRibo0_5TPM$`oblong/sphere_t4_polyA_49986`+polyAAndRibo0_5TPM$`oblong/sphere_t4_polyA_49990`)/3



t1_mean_riboCop = (polyAAndRibo0_5TPM$`4cell_t1_RiboCop_50007`+polyAAndRibo0_5TPM$`4cell_t1_RiboCop_50015`)/2
t2_mean_riboCop = (polyAAndRibo0_5TPM$`64-128cell_t2_RiboCop_50008`+polyAAndRibo0_5TPM$`64-128cell_t2_RiboCop_50012`+polyAAndRibo0_5TPM$`4cell_t1_RiboCop_50015`)/3
t3_mean_riboCop = (polyAAndRibo0_5TPM$`256-512cell_t3_RiboCop_50009`+polyAAndRibo0_5TPM$`256-512cell_t3_RiboCop_50013`+polyAAndRibo0_5TPM$`256-512cell_t3_RiboCop_50017`)/3
t4_mean_riboCop = (polyAAndRibo0_5TPM$`oblong/sphere_t4_RiboCop_50010`+polyAAndRibo0_5TPM$`oblong/sphere_t4_RiboCop_50014` + polyAAndRibo0_5TPM$`oblong/sphere_t4_RiboCop_50018`)/3

polyA_ribo0_mean  = cbind.data.frame(t1_mean_polyA,t2_mean_polyA,t3_mean_polyA,t4_mean_polyA,t1_mean_riboCop,t2_mean_riboCop,t3_mean_riboCop,t4_mean_riboCop)

colnames(polyA_ribo0_mean) = c("1hpf_mean_polyA","2.5hpf_mean_polyA","3.5_mean_polyA","4.5_mean_polyA","1hpf_mean_riboCop","2.5hpf_mean_riboCop","3.5hpf_mean_riboCop","4.5hpf_mean_riboCop")
polyA_ribo0_mean$geneName = row.names(polyAAndRibo0_5TPM)
relavantQuantSeq$geneName = row.names(relavantQuantSeq)

allData_together = join(polyA_ribo0_mean,relavantQuantSeq)

allData_together = allData_together[complete.cases(allData_together),]

row.names(allData_together)  = allData_together$geneName
allData_together$geneName = NULL

jpeg("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/comparisonOfPolyARiboZeroQuantSeq/correlation_quantSeq_polyASeq_riboSeq.jpg")
corrplot(cor(allData_together))
dev.off()

jpeg("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/comparisonOfPolyARiboZeroQuantSeq/correlation_quantSeq_polyASeq_riboSeq_spearman.jpg")
corrplot(cor(allData_together,method = "spearman"))
dev.off()

#### the pearson correlation of quantSeq data with polyA+ rnaseq is better than that of quantSeq with RiboCop. 

### I exoect quantSeq to correlate very well with polyA +, but the correlation is not great, maybe scatter plots will show this? 

### choose 3 common time points to the 3 datasets - 2.5hpf, 3.5hpf, 4.5 hpf

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/comparisonOfPolyARiboZeroQuantSeq/scatterplots_quantSeq_ribo0_polyA.pdf")
plot(log2(allData_together$`2.5hpf_mean_polyA`),log2(allData_together$`2.5hpf_quantSeq`),xlim=c(-5,10),ylim = c(-5,10))
plot(log2(allData_together$`3.5_mean_polyA`),log2(allData_together$`3.5hpf_quantSeq`),xlim=c(-5,10),ylim = c(-5,10))
plot(log2(allData_together$`4.5_mean_polyA`),log2(allData_together$`4.5hpf_quantSeq`),xlim=c(-5,10),ylim = c(-5,10))

plot(log2(allData_together$`2.5hpf_mean_riboCop`),log2(allData_together$`2.5hpf_quantSeq`),xlim=c(-5,10),ylim = c(-5,10))
plot(log2(allData_together$`3.5hpf_mean_riboCop`),log2(allData_together$`3.5hpf_quantSeq`),xlim=c(-5,10),ylim = c(-5,10))
plot(log2(allData_together$`4.5hpf_mean_riboCop`),log2(allData_together$`4.5hpf_quantSeq`),xlim=c(-5,10),ylim = c(-5,10))

dev.off()




### from the above plots it looks like there are several genes that have very low quantSeq experssion but high polyA+ expression 

### checking the identiy of these genes,

allData_together$foldChange_2.5hpf_polyA_quantSeq = allData_together$`2.5hpf_mean_polyA`/allData_together$`2.5hpf_quantSeq`
allData_together$foldChange_3.5hpf_polyA_quantSeq = allData_together$`3.5_mean_polyA`/allData_together$`3.5hpf_quantSeq`
allData_together$foldChange_4.5hpf_polyA_quantSeq = allData_together$`4.5_mean_polyA`/allData_together$`4.5hpf_quantSeq`

## this could be due to incomplete bad annotation that does not allow quantSeq estimation correctly. (check screenshots). 

# allCombinations = expand.grid(c(1:ncol(averageOfTimepoints_5tpm)),c(1:ncol(averageOfTimepoints_5tpm)))
# pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/plots/externalDataComparison/scatterplots_roboZeroVsPolyAplur.pdf")
# for(i in 1:nrow(allCombinations)){
#   rowGet = allCombinations[i,1]
#   colGet = allCombinations[i,2]
#   cor_two = round(cor(averageOfTimepoints_5tpm[,rowGet],averageOfTimepoints_5tpm[,colGet],method = "spearman"),2)
#   plot(log2(averageOfTimepoints_5tpm[,rowGet]),log2(averageOfTimepoints_5tpm[,colGet]),xlab = colnames(averageOfTimepoints_5tpm)[rowGet],ylab=colnames(averageOfTimepoints_5tpm)[colGet],main=cor_two)
#   
# }
# dev.off()
# 
# 
### check fold changes of t1_polyA vs t1_ribo0

