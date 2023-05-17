#17-Aug-2022 Diogo Ribeiro @ UNIL
# Script to plot TF sharing in enh-enh associations

library(data.table)
library(ggplot2)

options(scipen = 1)

enhEnhData = fread("~/git/enhEnh/source_data/pbmc_dataset/enh_enh_correlation_motifmap_positive.tsv.gz", header = T, sep = "\t")

enhEnhData$fdr = p.adjust(enhEnhData$pval, method = "fdr")

enhEnhData$significant = "no"
enhEnhData[corr > 0.05][fdr < 0.05]$significant = "yes"

table(enhEnhData$significant)

mergedData = enhEnhData

## Number of TF shared boxplot
ggplot(mergedData, aes(x=significant, y=N, fill = significant)) + 
  geom_boxplot( width = 0.5, size = 2) +
  labs(x="Enhancer-enhancer significant", y = "# shared TFs")+
  ylim (c(-1,max(mergedData$N)))+
  geom_text(aes(label=paste("N=",..count.., sep = "") ), y=-1, stat='count', colour="black", size=10)+
  stat_summary(fun=mean, geom="text", size=10, color="black", vjust = 0, aes(label= paste( round(..y.., digits = 2))), fontface = "bold")+
  scale_fill_manual(values = c("#B3AF8F","#66999B")) +
  theme_linedraw() + 
  theme(text = element_text(size=45), # for big screen
        legend.text=element_text(size=20), legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2), legend.position = "None", aspect.ratio = 1)  

t = wilcox.test(mergedData[significant == "yes"]$N, mergedData[significant == "no"]$N)
t
t$p.value


## Distance-Matched 
match = fread("../../source_data/pbmc_dataset/pbmc_sign_nonsign_enh_enh_match_positive.txt.gz")
simpleData = unique(mergedData[,.(tag,N)])
signMatch = simpleData[tag %in% match$sig]
signMatch$significant = "yes"
nonsignMatch = simpleData[tag %in% match$nonsig]
nonsignMatch$significant = "no"
matchDT = rbind(signMatch,nonsignMatch)
ggplot(matchDT, aes(x=significant, y=N, fill = significant)) + 
  geom_boxplot( width = 0.5, size = 1, outlier.alpha = 0.1) +
  labs(x="Enhancer-enhancer significant", y = "# shared TFs")+
  ylim (c(-1,max(mergedData$N)))+
  geom_text(aes(label=paste("N=",..count.., sep = "") ), y=-1, stat='count', colour="black", size=5)+
  stat_summary(fun=mean, geom="text", size=5, color="black", vjust = 0, aes(label= paste( round(..y.., digits = 2))), fontface = "bold")+
  scale_fill_manual(values = c("#B3AF8F","#66999B")) +
  theme_linedraw() + 
  theme(text = element_text(size=22),
        legend.text=element_text(size=20), legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position = "None", aspect.ratio = 1)  

t = wilcox.test(matchDT[significant == "yes"]$N, matchDT[significant == "no"]$N)
t
t$p.value


## Bin2D plot TF shared vs correlation
t = cor.test(mergedData$N, mergedData$corr, method = "spearman")
text = paste("Spearman R = ",round(t$estimate,2), " P-value ",format.pval(t$p.value), sep = "")
ggplot(mergedData, aes(x=corr, y=N)) + 
  geom_bin2d( bins = 15) +
  geom_smooth(method = "lm", size = 1.5) +
  labs(x="Enhancer-enhancer correlation", y = "# shared TFs")+
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.05, vjust = 1.5, size = 10, fontface = "bold"  ) +
  scale_fill_gradient(high = "black", low = "#f0f0f0") +
  theme_linedraw() + 
  theme(text = element_text(size=45), # for big screen
        legend.text=element_text(size=20), #legend.position = "None",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2), aspect.ratio = 1)

