## 2022 Chaymae Ziyani @ UNIL
# Produces input files to run Hi-C analysis on enhancer-enhancer interactions

library(data.table)

minTotalCells = 100

##### Read the input file
# EnhancerData = fread( "/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/enh_enh_paper/data/enh_enh_correlation.tsv", stringsAsFactors = FALSE, header = T, sep="\t")
EnhancerData = fread( "/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/enh_enh_paper/data/enh_enh_correlation_all.tsv.gz", stringsAsFactors = FALSE, header = T, sep="\t")
EnhancerData$totcell = rowSums(EnhancerData[, c(4:7), with = FALSE], na.rm=TRUE)
EnhancerData = EnhancerData[totcell >= minTotalCells] 
EnhancerData = unique(EnhancerData) 

##### Enhancer 1 process
EnhancerData$enhancer1_chr = data.table(unlist(lapply(EnhancerData$enh1, function(x) unlist(strsplit(x,"[_]"))[1])))$V1
EnhancerData$enhancer1_start = as.numeric(data.table(unlist(lapply(EnhancerData$enh1, function(x) unlist(strsplit(x,"[_]"))[2])))$V1)
EnhancerData$enhancer1_end = as.numeric(data.table(unlist(lapply(EnhancerData$enh1, function(x) unlist(strsplit(x,"[_]"))[3])))$V1)

##### Enhancer 2 process
EnhancerData$enhancer2_chr = data.table(unlist(lapply(EnhancerData$enh2, function(x) unlist(strsplit(x,"[_]"))[1])))$V1
EnhancerData$enhancer2_start = as.numeric(data.table(unlist(lapply(EnhancerData$enh2, function(x) unlist(strsplit(x,"[_]"))[2])))$V1)
EnhancerData$enhancer2_end = as.numeric(data.table(unlist(lapply(EnhancerData$enh2, function(x) unlist(strsplit(x,"[_]"))[3])))$V1)

##### Merge data by gene
mergedData = EnhancerData

## Midpoint for enhancer 1 and 2
mergedData$cisStart = (mergedData$enhancer1_start + mergedData$enhancer1_end) / 2
mergedData$cisEnd = mergedData$cisStart + 1
mergedData$cisStrand = "+"

mergedData$centralStart = (mergedData$enhancer2_start + mergedData$enhancer2_end) / 2
mergedData$centralEnd = mergedData$centralStart + 1
mergedData$centralStrand = "+"

mergedData$`#chr` = data.table(unlist(lapply(mergedData$enhancer1_chr, function(x) unlist(strsplit(x,"chr"))[2])))$V1 
mergedData = unique(mergedData[,.(`#chr`,centralStart,centralEnd,centralStrand,cisStart,cisEnd,cisStrand,corr,pval )])
mergedData$real = 1
mergedData$match_id = seq(nrow(mergedData))

###### Control part
controlData = mergedData
controlData$tss = controlData$centralStart
controlData[centralStrand == "+"]$tss = controlData[centralStrand == "+"]$centralEnd
controlData$dist = controlData$tss - controlData$cisStart
controlData$cisStart = controlData$cisStart + controlData$dist * 2
controlData$cisEnd = controlData$cisStart + 1
controlData = unique(controlData[,.(`#chr`,centralStart,centralEnd,centralStrand,cisStart,cisEnd,cisStrand,corr,pval)])
controlData$real = 0
controlData$match_id = seq(nrow(controlData))

finalDT = rbind(mergedData,controlData)
table(finalDT$real)

# write.table(finalDT,"/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/enh_enh_paper/data/Hi-C/input_for_hic_analysis_significant.bed", row.names = F, sep = "\t", quote = F)
write.table(finalDT,"/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/enh_enh_paper/data/Hi-C/input_for_hic_analysis_all.bed", row.names = F, sep = "\t", quote = F)
