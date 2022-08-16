
library(data.table)

inFile = "/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/enh_enh_paper/data/coex_peak_F0.5_coding_permuteR1000_fdr.tsv.gz"
peakAssociationFDRCutoff = 0.05
peakCorrCutoff = 0.05

data = fread( inFile, stringsAsFactors = FALSE, header = T, sep="\t")

data$gene = data.table(unlist(lapply(data$pairID, function(x) unlist(strsplit(x,"[|]"))[1])))$V1
data$peak = data.table(unlist(lapply(data$pairID, function(x) unlist(strsplit(x,"[|]"))[2])))$V1

length(unique(data$peak))
length(unique(data$gene))
table(data$corrSign)

data$fdr = p.adjust(data$corrPval, method = "BH")
data = data[,.(gene,peak,corr,corrSign,corrPval,fdr)]

#data[is.na(corr)]

sign = data[corr > peakCorrCutoff][fdr < peakAssociationFDRCutoff]

length(unique(sign$peak))
length(unique(sign$gene))
table(sign$corrSign)

write.table(sign[,.(gene,peak,corr,corrSign,corrPval)],"/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/enh_enh_paper/data/gene_enhancer_significant.tsv",quote=F,sep="\t",row.names=F)

dataFreq = data.table(table(data$gene))
signFreq = data.table(table(sign$gene))

mergedData = merge(dataFreq,signFreq,by="V1")
colnames(mergedData) = c("gene","total_enh","sign_enh")

mergedData[sign_enh > total_enh] # check

mergedData[sign_enh > 1]

write.table(mergedData,"/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/enh_enh_paper/data/gene_enhancer_frequency.tsv",quote=F,sep="\t",row.names=F)

