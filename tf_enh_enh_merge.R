#16-Aug-2022 Diogo Ribeiro @ UNIL

library(data.table)

minimumCells = 100

enhEnhData = fread("/home/dribeiro/EnhEnhPaper/data/enh_enh_correlation.tsv", header = T, sep = "\t")
enhEnhData$totalCells = enhEnhData$oneOne + enhEnhData$oneZero + enhEnhData$zeroOne + enhEnhData$zerozero
enhEnhData = enhEnhData[totalCells >= minimumCells]

## Exclude genes to simplify merge
enhEnhPairs = unique(enhEnhData[,.(enh1,enh2)])

tfData = fread("/home/dribeiro/EnhEnhPaper/data/TF/enhancers_sign_gene_enhancer_overlap_remap_F1.bed", header = F, sep= "\t")
tfData$tag = paste(tfData$V1,tfData$V2,tfData$V3,sep="_")
# Only care about the TFs, not the TFBS
tfData = unique(tfData[,.(tag,V7)])

# Merge, only care about cases where TF is the same between enh1 and enh2
m1 = merge(enhEnhPairs, tfData, by.x = "enh1", by.y = "tag", allow.cartesian = T)
m2 = merge(enhEnhPairs, tfData, by.x = "enh2", by.y = "tag", allow.cartesian = T)
mergedData = merge(m1,m2, by = c("V7","enh1","enh2"))
mergedData$tag = paste(mergedData$enh1,mergedData$enh2,sep="|")

# TF frequencies
freq = data.table(table(mergedData$tag))
# Account for zeroes
dt = unique(data.table(paste(enhEnhPairs$enh1,enhEnhPairs$enh2, sep="|")))
enhEnhFreq = merge(dt, freq, by = "V1", all.x = T)
enhEnhFreq[is.na(N)]$N = 0

## Remerge onto gene data
enhEnhData$tag = paste(enhEnhData$enh1,enhEnhData$enh2, sep="|")
finalDT = merge(enhEnhData, enhEnhFreq, by.x = "tag", by.y = "V1", all.x = T)

write.table(finalDT, "~/EnhEnhPaper/data/TF/enh_enh_correlation_remap.tsv", sep = "\t", quote=F, row.names=F)
