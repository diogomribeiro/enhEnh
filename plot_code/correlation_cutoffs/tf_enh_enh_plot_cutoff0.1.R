#17-Aug-2022 Diogo Ribeiro @ UNIL
# Script to plot TF sharing in enh-enh associations

library(data.table)
library(ggplot2)

options(scipen = 1)

# enhEnhData = fread("../../source_data/enh_enh_correlation_remap.tsv.gz", header = T, sep = "\t")
enhEnhData = fread("../../source_data/enh_enh_correlation_motifmap.tsv.gz", header = T, sep = "\t")

# signData = fread("../../source_data/correlation_cutoff/significant_enh_enh_cutoff0.1.tsv.gz", header = T, sep = "\t")
signData = fread("../../source_data/correlation_cutoff/significant_enh_enh_bonferroni0.05.tsv.gz", header = T, sep = "\t")

enhEnhData$triplet = paste(enhEnhData$gene,enhEnhData$enh1,enhEnhData$enh2,sep="|")
signData$triplet = paste(signData$gene,signData$enh1,signData$enh2,sep="|")

signTF = enhEnhData[triplet %in% signData$triplet]
nonsignTF = enhEnhData[!triplet %in% signData$triplet]

signTF$significant = "yes"
nonsignTF$significant = "no"
mergedData = rbind(signTF,nonsignTF)

table(mergedData$significant)

## Number of TF shared boxplot
ggplot(mergedData, aes(x=significant, y=N, fill = significant)) + 
  geom_boxplot( width = 0.5, size = 2) +
  labs(x="Enhancer-enhancer significant", y = "# shared TFs")+
  ylim (c(-1,max(mergedData$N)))+
  geom_text(aes(label=paste("N=",..count.., sep = "") ), y=-2, stat='count', colour="black", size=10)+
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
# match = fread("../../source_data/correlation_cutoff/sign_nonsign_enh_enh_match_cutoff0.1.txt.gz")
match = fread("../../source_data/correlation_cutoff/sign_nonsign_enh_enh_match_bonferroni.txt.gz")
simpleData = unique(mergedData[,.(tag,N)])
signMatch = simpleData[tag %in% match$sig]
signMatch$significant = "yes"
nonsignMatch = simpleData[tag %in% match$nonsig]
nonsignMatch$significant = "no"
matchDT = rbind(signMatch,nonsignMatch)
t = wilcox.test(signMatch$N,nonsignMatch$N)
t$p.value
ggplot(matchDT, aes(x=significant, y=N, fill = significant)) + 
  geom_boxplot( width = 0.5, size = 1) +
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

