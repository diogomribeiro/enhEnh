#17-Aug-2022 Diogo Ribeiro @ UNIL
# Script to plot TF sharing in enh-enh associations

library(data.table)
library(ggplot2)

options(scipen = 1)

enhEnhData = fread("~/EnhEnhPaper/data/TF/enh_enh_correlation_remap.tsv", header = T, sep = "\t")
# enhEnhData = fread("~/EnhEnhPaper/data/TF/enh_enh_correlation_motifmap.tsv", header = T, sep = "\t")
# enhEnhData = fread("~/EnhEnhPaper/data/TF/abc/enh_enh_correlation_remap.tsv", header = T, sep = "\t")
# enhEnhData = fread("~/EnhEnhPaper/data/TF/abc/enh_enh_correlation_motifmap.tsv", header = T, sep = "\t")

signData = fread("~/EnhEnhPaper/data/significant_enh_enh.tsv", header = T, sep = "\t")
# signData = fread("~/EnhEnhPaper/data/abc/significant_enh_enh.tsv", header = T, sep = "\t")

enhEnhData$triplet = paste(enhEnhData$gene,enhEnhData$enh1,enhEnhData$enh2,sep="|")
signData$triplet = paste(signData$gene,signData$enh1,signData$enh2,sep="|")

signTF = enhEnhData[triplet %in% signData$triplet]
nonsignTF = enhEnhData[!triplet %in% signData$triplet]

signTF$significant = "yes"
nonsignTF$significant = "no"
mergedData = rbind(signTF,nonsignTF)

## Number of TF shared boxplot
ggplot(mergedData, aes(x=significant, y=N, fill = significant)) + 
  geom_boxplot( alpha = 0.7, width = 0.5) +
  labs(x="Enhancer-enhancer significant", y = "# shared TFs")+
  ylim (c(-1,max(mergedData$N)))+
  # ylim (c(-5,max(mergedData$N)))+
  geom_text(aes(label=paste("N=",..count.., sep = "") ), y=-1, stat='count', colour="black", size=4.5)+
  # stat_summary(fun=mean, geom="point", size=2, color="red")+
  stat_summary(fun=mean, geom="text", size=5, color="black",
               vjust = 0, aes(label= paste( round(..y.., digits = 2))))+
  scale_fill_brewer(palette = "Set2") +
  theme_linedraw() + 
  theme(text = element_text(size=24),
        legend.text=element_text(size=20), legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position = "None", aspect.ratio = 1)  

t = wilcox.test(mergedData[significant == "yes"]$N, mergedData[significant == "no"]$N)
t
t$p.value


## Distance-Matched 
match = fread("~/EnhEnhPaper/data/sign_nonsign_enh_enh_match.txt")
# match = fread("~/EnhEnhPaper/data/abc/sign_nonsign_enh_enh_match.txt")
simpleData = unique(mergedData[,.(tag,N)])
signMatch = simpleData[tag %in% match$sig]
signMatch$significant = "yes"
nonsignMatch = simpleData[tag %in% match$nonsig]
nonsignMatch$significant = "no"
matchDT = rbind(signMatch,nonsignMatch)
t = wilcox.test(signMatch$N,nonsignMatch$N)
t$p.value
ggplot(matchDT, aes(x=significant, y=N, fill = significant)) + 
  geom_boxplot( alpha = 0.7, width = 0.5) +
  labs(x="Enhancer-enhancer significant", y = "# shared TFs")+
  ylim (c(-1,max(mergedData$N)))+
  geom_text(aes(label=paste("N=",..count.., sep = "") ), y=-1, stat='count', colour="black", size=4.5)+
  # stat_summary(fun=mean, geom="point", size=2, color="red")+
  stat_summary(fun=mean, geom="text", size=5, color="black",
               vjust = 0, aes(label= paste( round(..y.., digits = 2))))+
  scale_fill_brewer(palette = "Set2") +
  theme_linedraw() + 
  theme(text = element_text(size=24),
        legend.text=element_text(size=20), legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position = "None", aspect.ratio = 1)  


## Bin2D plot TF shared vs correlation
t = cor.test(mergedData$N, mergedData$corr, method = "spearman")
text = paste("Spearman R = ",round(t$estimate,2), " P-value ",format.pval(t$p.value), sep = "")
ggplot(mergedData, aes(x=corr, y=N)) + 
  # geom_bin2d( bins = 50) +
  geom_bin2d( bins = 15) +
  geom_smooth(method = "lm") +
  labs(x="Enhancer-enhancer correlation", y = "# shared TFs")+
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.05, vjust = 1.5, size = 5, fontface = "bold"  ) +
  scale_fill_gradient(high = "black", low = "#f0f0f0") +
  theme_linedraw() + 
  theme(text = element_text(size=24),
        legend.text=element_text(size=20), #legend.position = "None",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio = 1)


## Bin2D plot TF shared JI vs correlation
t = cor.test(mergedData$jaccard_index, mergedData$corr, method = "spearman")
text = paste("Spearman R = ",round(t$estimate,2), " P-value ",format.pval(t$p.value), sep = "")
ggplot(mergedData, aes(x=corr, y=jaccard_index)) + 
  # geom_bin2d( bins = 50) +
  geom_bin2d( bins = 30) +
  geom_smooth(method = "lm") +
  labs(x="Enhancer-enhancer correlation", y = "TF sharing jaccard index")+
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.05, vjust = 1.5, size = 5, fontface = "bold"  ) +
  scale_fill_gradient(high = "black", low = "#f0f0f0") +
  theme_linedraw() + 
  theme(text = element_text(size=24),
        legend.text=element_text(size=20), #legend.position = "None",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio = 1)

