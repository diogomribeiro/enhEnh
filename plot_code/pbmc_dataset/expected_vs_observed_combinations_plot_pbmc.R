# 09-Aug-2022 Diogo Ribeiro, Chaymae Ziyani @ UNIL
# Script to calculate observed vs expected associations

library(data.table)
library(ggplot2)

enhPerGene = fread("../../source_data/pbmc_dataset/pbmc_gene_enhancer_frequency_positive.tsv.gz", header = T, sep = "\t") 

summary(enhPerGene[sign_enh > 1]$sign_enh)

enhPerGene$expected = 2 ** (enhPerGene$sign_enh) - 1
summary(enhPerGene$expected)

### Enhancer-enhancer pairs
enhEnhData = fread("../../source_data/pbmc_dataset/pbmc_gene_enhancer_combinations_positive.out.gz", header = T, sep = "\t") 

observed = data.table(table(enhEnhData$gene))
colnames(observed) = c("gene","observed")

mergedData = merge(enhPerGene,observed,by="gene", all.x = T)
mergedData[is.na(observed)]$observed = 0

mergedData$observed = as.numeric(mergedData$observed)
mergedData$expected = as.numeric(mergedData$expected)

# sanity check
mergedData[observed > expected]

mergedData$ratio = mergedData$observed / mergedData$expected * 100
mergedData$NEnh = as.factor(mergedData$sign_enh)
mergedData[sign_enh >= 8][sign_enh <= 10]$NEnh = "8-10"
mergedData[sign_enh >= 10][sign_enh <= 15]$NEnh = "10-15"
mergedData[sign_enh > 15]$NEnh = ">15"

#### boxplot
ggplot(mergedData[!is.na(ratio)][NEnh != 1], aes(x=NEnh, y=ratio)) + 
  geom_boxplot( fill="#69b3a2", size = 1.5, width = 0.7, alpha = 0.7, outlier.shape = 1) +
  labs(x="Number of enhancers", y = "% combinations observed")+
  ylim (c(-5,100))+ 
  geom_text(aes(label=paste(..count.., sep="")), y=-7.5, stat='count', colour="black", size=10)+
  stat_summary(fun=mean, geom="point", size=7, color="#969696", shape = 19)+
  stat_summary(fun=mean, geom="text", size=10, color="black",vjust = 1.5, aes(label= paste( round(..y.., digits = 1))), fontface = "bold")+
  theme_linedraw() + 
  theme(text = element_text(size=44),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2), aspect.ratio = 1)  

#### density plot
meltedData = melt(mergedData[,.(expected,observed)])
ggplot(meltedData, aes(x=log10(value), fill = variable)) + 
  geom_density( adjust = 5, alpha = 0.5) +
  scale_fill_brewer(palette = "Set2")+
  xlab("log10(combinations)") +
  theme_linedraw() + 
  theme(text = element_text(size=24),
        legend.text=element_text(size=20), legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio = 1)  

summary(mergedData$ratio)
