### performing difffernetial exon analysis based on the annotations (counting windows)


### first I need to convert the .bed to GTF file :
### reading int the combined file of all counting windows  created in

allCountingWindows.bed = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/comparingAnnotationsBetweenStages/allCountingWindows_stagesCombined.bed",sep="\t",stringsAsFactors = F,header = T)

### based on the ensembl website this is the structure of the GTF file : 

#Fields must be tab-separated. Also, all but the final field in each feature line must contain a value; "empty" columns should be denoted with a '.'

  #seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
  #source - name of the program that generated this feature, or the data source (database or project name)
  #feature - feature type name, e.g. Gene, Variation, Similarity
  #start - Start position of the feature, with sequence numbering starting at 1.
  #end - End position of the feature, with sequence numbering starting at 1.
  #score - A floating point value.
  #strand - defined as + (forward) or - (reverse).
  #frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
  #attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.

seqname = allCountingWindows.bed$V1
source = "custom"
feature = "transcript"
start = allCountingWindows.bed$V2 + 1
end = allCountingWindows.bed$V3
score = "."
strand = allCountingWindows.bed$V6
frame="."

attribute = paste0('gene_id "',allCountingWindows.bed$V4,'";', ' gene_name "',allCountingWindows.bed$V4,'";', ' gene_source "',allCountingWindows.bed$stageCombine,'"')
completeGTF = cbind.data.frame(seqname,source,feature,start,end,score,strand,frame,attribute)
write.table(completeGTF,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/comparingAnnotationsBetweenStages/allStagesCountingWindows.gtf",sep="\t",quote = F)


GTFFile = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/comparingAnnotationsBetweenStages/allStagesCountingWindows.gtf",stringsAsFactors = F,sep="\t")
