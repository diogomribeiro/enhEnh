# 09-Aug-2022 Diogo Ribeiro, Chaymae Ziyani @ UNIL
# Script to calculate observed vs expected associations

library(data.table)
library(ggplot2)

enhPerGene = fread("../../source_data/200kb_window/200kb_gene_enhancer_frequency.tsv.gz", header = T, sep = "\t") 

enhPerGene$expectedPairs = enhPerGene$sign_enh * (enhPerGene$sign_enh - 1) / 2

summary(enhPerGene$expectedPairs)

### Enhancer-enhancer pairs
enhEnhData = fread("../../source_data/200kb_window/200kb_significant_enh_enh_cutoff.tsv.gz", header = T, sep = "\t")

observed = data.table(table(enhEnhData$gene))
colnames(observed) = c("gene","obs")

mergedData = merge(enhPerGene,observed,by="gene", all.x = T)
mergedData[is.na(obs)]$obs = 0

ggplot( mergedData, aes(x = expectedPairs, y = obs ))  +
  geom_point( size = 0.8, alpha = 0.5) +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed") +
  ylab("Observed pairs") +
  xlab("Maximum possible pairs") + 
  xlim(c(0,max(mergedData$expectedPairs))) +
  ylim(c(0,max(mergedData$expectedPairs))) +
  theme_linedraw() + 
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24), 
        legend.text=element_text(size=20), legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio = 1)

mergedData$ratio = mergedData$obs / mergedData$expectedPairs * 100
mergedData$NEnh = as.factor(mergedData$sign_enh)
mergedData[sign_enh >= 8][sign_enh <= 10]$NEnh = "8-10"
mergedData[sign_enh >= 10][sign_enh <= 15]$NEnh = "10-15"
mergedData[sign_enh > 15]$NEnh = ">15"

#### plot
ggplot(mergedData[!is.na(ratio)], aes(x=NEnh, y=ratio)) + 
  geom_boxplot( fill="#69b3a2", size = 1.5,  width = 0.7, alpha = 0.7, outlier.shape = 1) +
  labs(x="Number of enhancers", y = "% enhancer pairs significant")+
  ylim (c(-5,100))+ 
  geom_text(aes(label=paste(..count.., sep = "") ), y=-7.5, stat='count', colour="black", size=10)+
  stat_summary(fun=mean, geom="point", size=7, color="#969696", shape = 19)+
  stat_summary(fun=mean, geom="text", size=10, color="black",vjust = 1.5, aes(label= paste( round(..y.., digits = 1))), fontface = "bold")+
  theme_linedraw() + 
  theme(text = element_text(size=44),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2), aspect.ratio = 1)  

t = cor.test(mergedData$ratio,mergedData$sign_enh,method = "spearman")
t
t$p.value

