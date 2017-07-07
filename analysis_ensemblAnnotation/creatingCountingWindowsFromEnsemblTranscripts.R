### creating counting windows from ENSEMBL zebrafish transctipt annotation and merging windows <250nts : 

library(dplyr)
ensemblTranscriptAnnotation = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10/ensembl_dr10_Ensembl_Genes_88/transcriptStartsAndEnds_all.polyABiotype.txt",sep="\t",stringsAsFactors = F,header = T)


ensemblTranscriptAnnotation = ensemblTranscriptAnnotation %>% mutate(end_CountingWindow = ifelse(strand == "+", transcript_end , transcript_start + 250)) %>% mutate(start_countingWindows = ifelse(strand == "+",transcript_end-250,transcript_start )) %>% mutate(score=0)



ensemblTranscriptAnnotation = ensemblTranscriptAnnotation[,c("chromosome_name","start_countingWindows","end_CountingWindow","external_gene_name","score","strand","ensembl_gene_id","ensembl_transcript_id")]

### converting to 1 based notation 

ensemblTranscriptAnnotation$start_countingWindows = ensemblTranscriptAnnotation$start_countingWindows + 1

ensemblTranscriptAnnotation_plus = ensemblTranscriptAnnotation %>% filter(strand=="+")
annotation_custom_positive_split = split(ensemblTranscriptAnnotation_plus,f = ensemblTranscriptAnnotation_plus$external_gene_name,drop = T )

ensemblTranscriptAnnotation_minus = ensemblTranscriptAnnotation %>% filter(strand=="-")
annotation_custom_negative_split = split(ensemblTranscriptAnnotation_minus,f = ensemblTranscriptAnnotation_minus$external_gene_name,drop = T )


library(GenomicRanges)
positive_ranges = lapply(annotation_custom_positive_split,function(x) with(x,GRanges(chromosome_name, IRanges(start = start_countingWindows,end = end_CountingWindow),strand = strand,score=0,names=external_gene_name,geneId=ensembl_gene_id,transcriptId=ensembl_transcript_id)))
negative_ranges = lapply(annotation_custom_negative_split,function(x) with(x,GRanges(chromosome_name, IRanges(start = start_countingWindows,end = end_CountingWindow),strand = strand,score=0,names=external_gene_name,geneId=ensembl_gene_id,transcriptId=ensembl_transcript_id)))

### reducing the ranges : 

allAnnotations_plus_ranges_reduced = lapply(positive_ranges,function(x) unlist(reduce(split(x, elementMetadata(x)$names))) )



reducedToDf = function(reduced){
  reduced <- data.frame(seqnames=seqnames(reduced),
                        starts=start(reduced),
                        ends=end(reduced),
                        names=c(names(reduced)),
                        scores=0,strand = strand(reduced))
  return(reduced)
}

allAnnotations_plus_ranges_reduced_df = lapply(allAnnotations_plus_ranges_reduced,function(x) reducedToDf(x))
allAnnotations_plus_ranges_reduced_df = do.call(rbind,allAnnotations_plus_ranges_reduced_df)




allAnnotations_minus_ranges_reduced = lapply(negative_ranges,function(x) unlist(reduce(split(x, elementMetadata(x)$names))) )


allAnnotations_minus_ranges_reduced_df = lapply(allAnnotations_minus_ranges_reduced,function(x) reducedToDf(x))
allAnnotations_minus_ranges_reduced_df = do.call(rbind,allAnnotations_minus_ranges_reduced_df)



## all annotations_merged 

countingWindows = rbind(allAnnotations_plus_ranges_reduced_df,allAnnotations_minus_ranges_reduced_df)
countingWindows$starts = countingWindows$starts -1
chromosomes_mm = paste0("chr",c("2L","2R","3L","3R","4","X","Y"))
write.table(countingWindows, "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/analysis_ensemblAnnotation/countingWindoes_ensemblAnnotation.bed",sep="\t",quote = F,row.names = F,col.names = F)
