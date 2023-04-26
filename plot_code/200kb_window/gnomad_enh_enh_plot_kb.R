#17-Aug-2022 Diogo Ribeiro @ UNIL
# Script to plot gnomad/essentiality related plots

library(data.table)
library(ggplot2)

gnomadData = fread("../source_data/gnomad_gene_level_relevant_columns.txt.gz")
geneBackground = fread("../../source_data/200kb_window/200kb_gene_enhancer_frequency.tsv.gz")
geneBackground = geneBackground[sign_enh > 0]
gnomadData = gnomadData[gene_id %in% geneBackground$gene][!is.na(oe_lof_upper)]

## Significant enh-enh
enhEnhData = fread("../../source_data/200kb_window/200kb_significant_enh_enh_cutoff.tsv.gz", header = T, sep = "\t")

observed = data.table(table(enhEnhData$gene))
colnames(observed) = c("gene","observed")
observed = merge(geneBackground, observed, by = "gene", all.x = T)
observed[is.na(observed)]$observed = 0
mergedData = merge(observed,gnomadData, by.x="gene",by.y = "gene_id")

t = cor.test(mergedData$oe_lof_upper,mergedData$observed,method = "spearman")
text = paste("Spearman R = ",round(t$estimate,2),"P-value = ",signif(t$p.value, digits = 3))
ggplot(mergedData, aes(x=oe_lof_upper, y=observed) ) + 
  geom_bin2d( bins = 50) +
  # geom_bin2d( bins = 15) +
  # geom_point( fill="#69b3a2", size = 0.5, alpha = 0.5) +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.05, vjust = 1.5, size = 5, fontface = "bold"  ) +
  geom_smooth(method = "lm") +
  labs(x="LOEUF score", y = "Significant enhancer pairs")+
  scale_fill_gradient(trans = "log10", high = "#3f007d", low = "#dadaeb") +
  theme_linedraw() + 
  theme(text = element_text(size=24),
        legend.text=element_text(size=20), legend.title=element_blank(), #legend.position = "None",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio = 1)  
