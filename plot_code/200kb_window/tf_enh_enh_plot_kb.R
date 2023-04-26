#17-Aug-2022 Diogo Ribeiro @ UNIL
# Script to plot TF sharing in enh-enh associations

library(data.table)
library(ggplot2)

options(scipen = 1)

enhEnhData = fread("../source_data/enh_enh_correlation_remap.tsv.gz", header = T, sep = "\t")

filterData = fread("../../source_data/200kb_window/200kb_enh_enh_correlation.tsv", header = T, sep = "\t")

enhEnhData$p = paste0(enhEnhData$gene,"_",enhEnhData$enh1,"_",enhEnhData$enh2)
filterData$p = paste0(filterData$gene,"_",filterData$enh1,"_",filterData$enh2)
enhEnhData = enhEnhData[p %in% filterData$p]

enhEnhData$significant = "no"
enhEnhData[p %in% filterData[corr > 0.05][fdr < 0.05]$p]$significant = "yes"

mergedData = enhEnhData

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
match = fread("../../source_data/200kb_window/200kb_sign_nonsign_enh_enh_match.txt.gz")
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

