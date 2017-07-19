#### checking custom annotations of different developmental stages :


### the stages are : 
  # 1dpf
  # 256 cell stage
  # 2 cell 
  # 2dpf
  # 4pdf
  # bud
  # dome
  # oocyte
  # sphere
  # testis

path_allStages = "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/2016_12_QuantSeq_Annotaton/raw_data/"
stages = c("1dpf","256cell","2cell","2dpf","4dpf","bud","dome","oocyte","sphere","testis")
completePath = paste0(path_allStages,stages,"/output/final90percent/")

##### comparing high confidence ends from the different stages

completePath_highConfidenceEnds = paste0(completePath,"ends_greater90percent_intergenic_n100.bed")

highConfidenceEnds = lapply(completePath_highConfidenceEnds,function(x) read.delim(x,stringsAsFactors = F,header = F))
names(highConfidenceEnds) = stages
