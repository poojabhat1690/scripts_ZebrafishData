
#### reading in the gene counts that Thomas has created at some time from Andi. She has performed polyA+ RNAseq and Ribo0 seq from the same smaples. 
#### I will use this for the comparison....
#t1 = 4-cell embryos (~1 hour post fertilization (hpf))
#t2 = 64-128 cell (~2.5 hpf)
#t3 = 256-512 cell (~3.5 hpf)
#t4 = oblong/sphere (~4.5 hpf)



polyAAndRibo0 = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/externalDatat/fromAndi_ribo0VsPolyA/tpm_genes.txt",header=T)

sampleNames = strsplit(colnames(polyAAndRibo0),"_",T)
sampleNames = unlist(lapply(sampleNames,function(x) paste(x[2],x[3],x[4],sep = "_")))
sampleNames = sampleNames[-1]
times = rep(c("1hpf","2.5hpf","3.5hpf","4.5hpf"),each=3)
times = paste(times,sampleNames,sep="_")
rownames(polyAAndRibo0) = polyAAndRibo0$geneID
polyAAndRibo0$geneID = NULL
colnames(polyAAndRibo0) = times
cor_polyAandRibo0 = cor(polyAAndRibo0)
corrplot(cor_polyAandRibo0)


maxTPM = apply(polyAAndRibo0,1,max)
polyAAndRibo0_5tpm = polyAAndRibo0[which(maxTPM>5),]
cor_polyAandRibo0 = cor(polyAAndRibo0_5tpm)
corrplot(cor_polyAandRibo0)
pairs(polyAAndRibo0_5tpm)

#### now use quantSeq 


countDataSets_quantSeq = list.files(path = "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/analysis_ensemblAnnotation/countMatrices/",pattern = ".tsv")
countDataSets_quantSeq_paths = paste0( "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/analysis_ensemblAnnotation/countMatrices/",countDataSets_quantSeq)
countDataSets_quantSeq_data = lapply(countDataSets_quantSeq_paths,function(x) read.table(x,sep="\t",header=T))
names(countDataSets_quantSeq_data) = countDataSets_quantSeq
names(countDataSets_quantSeq_cpm) = countDataSets_quantSeq
countDataSets_quantSeq_cpm = lapply(countDataSets_quantSeq_data,function(x) x$ReadsCPM)
countDataSets_quantSeq_cpm = do.call(cbind.data.frame,countDataSets_quantSeq_cpm)
countDataSets_quantSeq_cpm$names = countDataSets_quantSeq_data[[1]]$Name


a = split( countDataSets_quantSeq_cpm,countDataSets_quantSeq_cpm$names,T )
a_summed = lapply(a,function(x) colSums(x[1:24]))
a_summed = do.call(rbind,a_summed)
a_summed = as.data.frame(a_summed)
a_summed$ensembl_gene_id = row.names(a_summed)

samplesNames_quantSeq = read.delim("/Volumes//groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/sampleInfo.txt",header = F,stringsAsFactors = F)
colnames(a_summed) = c(samplesNames_quantSeq$V3,"GeneName")


polyAAndRibo0$GeneName = row.names(polyAAndRibo0)
allData = join(polyAAndRibo0,a_summed)
allData = allData[complete.cases(allData),]
rownames(allData) = allData$GeneName 
allData$GeneName = NULL
allData_corr = cor(allData)
corrplot(allData_corr)
