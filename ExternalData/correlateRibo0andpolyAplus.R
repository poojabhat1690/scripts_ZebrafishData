#### I want to check if the maternal isoforms are influenced by the short polyA tails they have. 

## since QuantSeq is an approach based on oligoDT priming, it could be possible that sequences with short polyA sequences donot get primed well.


## to ceck this I want to comapre the expression of QUantSeq to riobo0 RNAseq. 

  ## the first data set i used for this was wild type samples from thw paper of mishimaAndTomari 2016. 


  #### first I  want to do a gene specific counting based on RNA seq data : 


### loading edge R to caluclate RPM data : 

source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")


library(edgeR)

### reading in RNAseq data from polyA+ RNA (from andi) and Ribo0 rnaSeq (from Mashima Tomari)

counts_polyAplus = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/rnaSeqReads_featureCounts_strandSpecific_allGenes.txt",header=T,stringsAsFactors = F)
counts_ribozero = read.table("/Volumes//groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/externalDatat/countsRNAseq_mashimaTomari.txt",header = T,stringsAsFactors = F)

sampleInfo_mashima = read.table("/Volumes//groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/externalDataResources/mashimaTomari/sampleInfo_emblEbi.txt",sep="\t",header=T)
countMatrices_polyAplus = counts_polyAplus[,7:ncol(counts_polyAplus)]
countMatrices_riboZero = counts_ribozero[,7:ncol(counts_ribozero)]

colnames(countMatrices_polyAplus) = c("1KCell_1","2.4Cell_1","2.4Cell_2","256Cell_1","28hpf_1","28hpf_2","2dpf_1","2dpf_2","5dpf_1","5dpf_2","Bud_1","Bud_2","Dome_1","Dome_2","1Kcell_2","Oblong","Shield_1","Shield_2","Shield_3","oocyte_1_5","oocyte_1_6","oocyte_1_7","oocyte_2_5", "oocyte_2_6","oocyte_2_7","ovary_1_5","ovary_1_6","ovary_1_7","sperm1_longRNA","sperm2_longRNA","testis_1_5","testis_1_6","testis_1_7","testis_2_5","testis_2_6","testis_2_7")
colnames(countMatrices_riboZero) = paste(sampleInfo_mashima$sampleDescription[1:10],c(1:10),sep="_")

### subset this to look only at the relavant stages 

#countMatrices_riboZero = countMatrices_riboZero[,1:6]
colnames(countMatrices_riboZero) = paste("R",colnames(countMatrices_riboZero),sep="_")
#countMatrices_polyAplus = countMatrices_polyAplus[,c(1:4,11:15)]
colnames(countMatrices_polyAplus) =paste("P",colnames(countMatrices_polyAplus),sep="_")
countMatrices_riboZero_rpkm = rpkm(x = countMatrices_riboZero,gene.length =counts_ribozero$Length )
countMatrices_polyAplus_rpkm = rpkm(x = countMatrices_polyAplus,gene.length = counts_polyAplus$Length)


total_RNAseq = cbind(countMatrices_riboZero_rpkm,countMatrices_polyAplus_rpkm)
library(corrplot)
corrplot(cor(total_RNAseq))

### now I want to correlate this with gene expression from quantSeq data : 


countDataSets_quantSeq = list.files(path = "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/analysis_ensemblAnnotation/countMatrices/",pattern = ".tsv")
countDataSets_quantSeq_paths = paste0( "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/analysis_ensemblAnnotation/countMatrices/",countDataSets_quantSeq)
countDataSets_quantSeq_data = lapply(countDataSets_quantSeq_paths,function(x) read.table(x,sep="\t",header=T))
names(countDataSets_quantSeq_data) = countDataSets_quantSeq
names(countDataSets_quantSeq_cpm) = countDataSets_quantSeq
countDataSets_quantSeq_cpm = lapply(countDataSets_quantSeq_data,function(x) x$ReadsCPM)
countDataSets_quantSeq_cpm = do.call(cbind.data.frame,countDataSets_quantSeq_cpm)
countDataSets_quantSeq_cpm$names = countDataSets_quantSeq_data[[1]]$Name
