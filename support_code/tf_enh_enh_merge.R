#16-Aug-2022 Diogo Ribeiro @ UNIL
# Adding TF information of enh-enh data

library(data.table)

minimumCells = 100

enhEnhData = fread("/home/dribeiro/EnhEnhPaper/data/abc/enh_enh_correlation.tsv", header = T, sep = "\t")
enhEnhData$totalCells = enhEnhData$oneOne + enhEnhData$oneZero + enhEnhData$zeroOne + enhEnhData$zerozero
enhEnhData = enhEnhData[totalCells >= minimumCells]

## Exclude genes to simplify merge
enhEnhPairs = unique(enhEnhData[,.(enh1,enh2)])

# tfData = fread("/home/dribeiro/EnhEnhPaper/data/TF/enhancers_sign_gene_enhancer_overlap_remap_F1.bed", header = F, sep= "\t")
# tfData = fread("/home/dribeiro/EnhEnhPaper/data/TF/enhancers_sign_gene_enhancer_overlap_motifmap_F1.bed", header = F, sep= "\t")
# tfData = fread("/home/dribeiro/EnhEnhPaper/data/TF/abc/enhancers_sign_gene_enhancer_overlap_remap_F1.bed", header = F, sep= "\t")
tfData = fread("/home/dribeiro/EnhEnhPaper/data/TF/abc/enhancers_sign_gene_enhancer_overlap_motifmap_F1.bed", header = F, sep= "\t")
tfData$tag = paste(tfData$V1,tfData$V2,tfData$V3,sep="_")
# Only care about the TFs, not the TFBS
tfData = unique(tfData[,.(tag,V7)])
totDT = data.table(table(tfData$tag))

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

colnames(totDT) = c("enh1","enh1_totN")
m = merge(finalDT,totDT,by = "enh1",all.x=T)
colnames(totDT) = c("enh2","enh2_totN")
finalDT = merge(m,totDT,by="enh2",all.x=T)
finalDT[is.na(enh1_totN)]$enh1_totN = 0
finalDT[is.na(enh2_totN)]$enh2_totN = 0
finalDT$union = finalDT$N + (finalDT$enh1_totN - finalDT$N) + (finalDT$enh2_totN - finalDT$N)
finalDT$jaccard_index = finalDT$N / finalDT$union

# write.table(finalDT, "~/EnhEnhPaper/data/TF/enh_enh_correlation_remap.tsv", sep = "\t", quote=F, row.names=F)
# write.table(finalDT, "~/EnhEnhPaper/data/TF/enh_enh_correlation_motifmap.tsv", sep = "\t", quote=F, row.names=F)
# write.table(finalDT, "~/EnhEnhPaper/data/TF/abc/enh_enh_correlation_remap.tsv", sep = "\t", quote=F, row.names=F)
write.table(finalDT, "~/EnhEnhPaper/data/TF/abc/enh_enh_correlation_motifmap.tsv", sep = "\t", quote=F, row.names=F)
