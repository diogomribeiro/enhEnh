#15-Aug-2022 Diogo Ribeiro @ UNIL
# Script to measure enhancer-enhancer associations across genes

library(data.table)
library(ggplot2)

enhEnhData = fread("../../source_data/pbmc_dataset/pbmc_enh_enh_correlation_positive.tsv.gz", header = T, sep = "\t") 
signData = fread("../../source_data/pbmc_dataset/pbmc_significant_enh_enh_positive.tsv.gz", header = T, sep = "\t") 

length(unique(enhEnhData$tag))

enhEnhData$tag = paste(enhEnhData$enh1,enhEnhData$enh2)
possible = data.table(table(enhEnhData$tag))

signData$tag = paste(signData$enh1,signData$enh2)
significant = data.table(table(signData$tag))

mergedData = merge(possible, significant, by = "V1", all.x = T)
colnames(mergedData) = c("enh_pair","possible","significant")
mergedData[is.na(significant)]$significant = 0

mergedData$ratio = mergedData$significant * 100 / mergedData$possible
mergedData$n_genes = as.factor(mergedData$possible)
# mergedData[possible >= 8][possible < 10]$n_genes = "8-10"
# mergedData[possible >= 10][possible < 15]$n_genes = "10-15"
# mergedData[possible >= 15]$n_genes = "=>15"
mergedData[possible >= 10]$n_genes = ">=10"


ggplot(mergedData, aes(x=n_genes, y=ratio)) + 
  geom_boxplot( fill="#bebada", size = 1.5, width = 0.7, alpha = 0.7, outlier.shape = 1) +
  labs(x="Number of genes", y = "% genes significant")+
  ylim (c(-5,100))+ 
  geom_text(aes(label=paste(..count.., sep = "") ), y=-7.5, stat='count', colour="black", size=10)+
  stat_summary(fun=mean, geom="point", size=7, color="#969696", shape = 19)+
  stat_summary(fun=mean, geom="text", size=10, color="black",vjust = 1.5, aes(label= paste( round(..y.., digits = 1))), fontface = "bold")+
  theme_linedraw() + 
  theme(text = element_text(size=44), # for big screen
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),  aspect.ratio = 1)  

t = cor.test(mergedData$possible,mergedData$ratio, method = "spearman")
t
t$p.value
