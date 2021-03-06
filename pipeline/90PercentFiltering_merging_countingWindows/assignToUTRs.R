









## this script creates a file for the accepted refSeq and ensembl 3'UTR overlapping ends. 


library(checkmate)
library(GenomicRanges)
### the input should be info about the peaks that overlap with refSeq and ensembl 

refSeq_total = read.table(paste0(BOut, "/polyAmapping_allTimepoints/n_100_global_a0/refSeq_total.bed"),sep="\t",stringsAsFactors = F)

ensembl_total = read.table(paste0(BOut, "/polyAmapping_allTimepoints/n_100_global_a0/ensembl_total.bed"),sep="\t",stringsAsFactors = F)

assertDataFrame(x = refSeq_total,ncols = 21)
assertDataFrame(x = ensembl_total,ncols = 21)

refSeq_total$overlap = "refSeq"
ensembl_total$overlap = "ensembl"

refSeq_pas = refSeq_total[-grep("noPAS",refSeq_total$V21),]
refSeq_noPas = refSeq_total[grep("noPAS",refSeq_total$V21),]


ensembl_pas = ensembl_total[-grep("noPAS",ensembl_total$V21),]
ensembl_noPas = ensembl_total[grep("noPAS",ensembl_total$V21),]

total_pas = rbind(refSeq_pas,ensembl_pas)

total_nopas = rbind(refSeq_noPas, ensembl_noPas)

## filtering based on the A threshold . here V10 is contains the A thresholds


total_pas_accepted =  total_pas[which(total_pas$V10<0.36),]
total_nopas_accepted = total_nopas[which(total_nopas$V10<0.24),]


total_pas_notAccetped = total_pas[which(total_pas$V10>=0.36),]
total_noPas_notAccepted =   total_nopas[which(total_nopas$V10>=0.24),]

allRejected = rbind(total_noPas_notAccepted, total_pas_notAccetped)
allRejected_granges = makeGRangesFromDataFrame(allRejected,keep.extra.columns = T,ignore.strand = F,seqnames.field = "V1",start.field = "V2",end.field = "V3",strand.field = "V6",starts.in.df.are.0based = T )


#### I want to now overlap this with ends from 3pSeq to identify possibly missed ends... 

ulitskyData = read.table("//groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/externalDatat/3primeEnds_convertedTobed_dr10.bed",sep="\t",stringsAsFactors = F)

ulitskyData_ranges = makeGRangesFromDataFrame(df = ulitskyData,keep.extra.columns = T,ignore.strand = F,seqnames.field = "V1",start.field = "V2",end.field = "V3",strand.field = "V6",starts.in.df.are.0based = T)



allRejected_granges = makeGRangesFromDataFrame(allRejected,keep.extra.columns = T,ignore.strand = F,seqnames.field = "V1",start.field = "V2",end.field = "V3",strand.field = "V6",starts.in.df.are.0based = T )

allOverlaps_rejection = findOverlaps(allRejected_granges,ulitskyData_ranges)

rejectedAndOverlapping = allRejected[queryHits(allOverlaps_rejection),]
rejectedAndOverlapping = rejectedAndOverlapping[!duplicated(rejectedAndOverlapping),]


allPeaks = rbind(total_pas_accepted, total_nopas_accepted)

### now also including ends that are overlapping compared to the 3pseq data

allPeaks = rbind(allPeaks,rejectedAndOverlapping)


peaks_accepted_nonAccepted = rbind(total_pas, total_nopas)

acceptedPeaks_sum = sum(as.numeric(unlist(strsplit(allPeaks$V5,split = ",",fixed = T))))
peaks_accepted_nonAccepted_sum = sum(as.numeric(unlist(strsplit(peaks_accepted_nonAccepted$V5,split = ",",fixed = T))))


cat(paste0("The fraction of counts accepted:",acceptedPeaks_sum/peaks_accepted_nonAccepted_sum))



## now checking if the peaks are at the same position as the UTRs. 



# 
# allPeaks_positive = allPeaks[which(allPeaks$V6 == "+"),]
# allPeaks_negative = allPeaks[which(allPeaks$V6 == "-"),]
# 
# distanceFromUTR_minus = allPeaks_negative$V12 - allPeaks_negative$V2
# distanceFromUTR_minus[which(distanceFromUTR_minus!=0)] <-"short"
# distanceFromUTR_minus[which(distanceFromUTR_minus==0)] <-"equal"
# 
# distanceFromUTR_plus = allPeaks_positive$V13 - allPeaks_positive$V3
# distanceFromUTR_plus[which(distanceFromUTR_plus!=0)] <-"short"
# distanceFromUTR_plus[which(distanceFromUTR_plus==0)] <-"equal"
# 
# allPeaks_positive$distance = distanceFromUTR_plus
# allPeaks_negative$distance = distanceFromUTR_minus
# 
# allPeaks = rbind(allPeaks_positive,allPeaks_negative)



## writing these tables 


write.table(allPeaks,paste0(BOut, "/polyAmapping_allTimepoints/n_100_global_a0/ends_all_10threshold_n100.txt"),sep="\t",quote = F,row.names = F)

#write.table(allPeaks,"/groups/ameres/Pooja/Projects/wholePipeline/otherFiles/ends_all_10threshold_n100.txt",sep = "\t",quote = F,row.names = F)
#write.table(allPeaks,"/Users/pooja.bhat/Dropbox/UTRannotation/mESC/nucleotideProfiles_revised/OverlappingPrimingSitesWithAnnotations/data/utrOverlappin_AthresholdPassingEnds.txt",sep = "\t",quote = F,row.names = F)

### writng a part of this to a bed file : 

allPeaks$V4 = allPeaks$V14

allPeaks_short = cbind.data.frame(allPeaks[,c(1:6)],allPeaks$overlap)
write.table(allPeaks_short,paste0(BOut, "/polyAmapping_allTimepoints/n_100_global_a0/ends_all_10threshold_n100.bed"),sep = "\t",quote = F,row.names = F,col.names = F)



