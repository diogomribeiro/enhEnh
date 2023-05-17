#16-May-2023 Diogo Ribeiro @ UNIL

library(data.table)
library(ggplot2)
library(patchwork)
library(ggsignif)

pdf("~/Figure3.pdf",18,18)

options(scipen = 1)

###################
# Panel A
##################
mergedData = fread("~/git/enhEnh/source_data/data_figure_3a.tsv.gz", header = T, sep="\t")

## Number of TF shared boxplot
g1 = ggplot(mergedData, aes(x=significant, y=N, fill = significant)) + 
  geom_boxplot( width = 0.5, size = 1.5) +
  geom_text(aes(label=paste("N=",..count.., sep = "") ), y=-2, stat='count', colour="black", size=8)+
  geom_signif(y_position = 90, xmin = 1.1, xmax = 1.9, textsize = 8, annotation = c("p<2e-308"), size = 1.5) + 
  stat_summary(fun=mean, geom="text", size=8, color="black", vjust = 0, aes(label= paste( round(..y.., digits = 2))), fontface = "bold")+
  scale_fill_manual(values = c("#B3AF8F","#66999B")) +
  labs(x="Enhancer-enhancer significant", y = "# shared TFs")+
  labs(tag = "a") +
  ylim (c(-1,max(mergedData$N)))+
  theme_linedraw() + 
  theme(text = element_text(size=32), # for big screen
        legend.text=element_text(size=17), legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2), legend.position = "None", aspect.ratio = 1)  

t = wilcox.test(mergedData[significant == "yes"]$N, mergedData[significant == "no"]$N)
t
t$p.value


##################
# Panel B
##################

## Bin2D plot TF shared vs correlation
t = cor.test(mergedData$N, mergedData$corr, method = "spearman")
pval = t$p.value 
if(t$p.value == 0){ pval = 2.2e-308 }
text = paste("Spearman R = ",round(t$estimate,2), " P-value <",pval, sep = "")
g2 = ggplot(mergedData, aes(x=corr, y=N)) + 
  geom_bin2d( bins = 50) +
  geom_smooth(method = "lm", size = 1.5) +
  labs(tag = "b") +
  labs(x="Enhancer-enhancer correlation", y = "# shared TFs")+
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.05, vjust = 1.5, size = 8, fontface = "bold"  ) +
  scale_fill_gradient(high = "black", low = "#f0f0f0") +
  theme_linedraw() + 
  theme(text = element_text(size=32), # for big screen
        legend.text=element_text(size=17), legend.title=element_blank(), legend.position = c(0.85,0.81),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2), aspect.ratio = 1)


##################
# Panel C
##################

mergedData = fread("~/git/enhEnh/source_data/data_figure_3c.tsv.gz", header=T, sep="\t")
# write.table(mergedData, "~/git/enhEnh/source_data/data_figure_3c.tsv", quote=F, sep="\t", row.names=F)

t = cor.test(mergedData$oe_lof_upper,mergedData$observed,method = "spearman")
text = paste("Spearman R = ",round(t$estimate,2),"P-value =",signif(t$p.value, digits = 3))
g3 = ggplot(mergedData, aes(x=oe_lof_upper, y=observed) ) + 
  geom_bin2d( bins = 50) +
  geom_smooth(method = "lm") +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.05, vjust = 1.5, size = 8, fontface = "bold"  ) +
  labs(tag = "c") +
  labs(x="LOEUF score", y = "Observed combinations")+
  scale_fill_gradient(trans = "log10", high = "#3f007d", low = "#dadaeb") +
  theme_linedraw() + 
  theme(text = element_text(size=32), # for big screen
        legend.text=element_text(size=17), legend.title=element_blank(), legend.position = c(0.85,0.83),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2), aspect.ratio = 1)  

##################
# Panel D
##################

t = cor.test(mergedData$oe_lof_upper,mergedData$total_enh,method = "spearman")
text = paste("Spearman R = ",round(t$estimate,2),"P-value =",signif(t$p.value, digits = 3))
g4 = ggplot(mergedData, aes(x=oe_lof_upper, y=total_enh) ) + 
  geom_bin2d( bins = 50) +
  labs(tag = "d") +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.05, vjust = 1.5, size = 8, fontface = "bold"  ) +
  geom_smooth(method = "lm") +
  scale_fill_gradient(trans = "log10", high = "#3f007d", low = "#dadaeb") +
  labs(x="LOEUF score", y = "# nearby enhancers")+
  theme_linedraw() + 
  theme(text = element_text(size=32), # for big screen
        legend.text=element_text(size=17), legend.title=element_blank(), legend.position = c(0.85,0.83), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2), aspect.ratio = 1)  


g1 + g2 + g3 + g4

dev.off()
