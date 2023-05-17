#11-Aug-2022 Diogo Ribeiro @ UNIL
# Comparison of rep2 and rep3 enhancer-enhancer correlations

library(data.table)
library(ggplot2)

minimumCells = 100

rep3Sign = fread("~/git/enhEnh/source_data/enh_enh_correlation.tsv.gz", header = T, sep = "\t")
rep3Sign$totalCells = rep3Sign$oneOne + rep3Sign$oneZero + rep3Sign$zeroOne + rep3Sign$zerozero
rep3Sign = rep3Sign[totalCells >= minimumCells]

rep3Sign = rep3Sign[,.(gene,enh1,enh2,corr,pval)]

rep2Data = fread("~/git/enhEnh/source_data/rep2_enh_enh_correlation.tsv.gz") 

rep2Data$totalCells = rep2Data$oneOne + rep2Data$oneZero + rep2Data$zeroOne + rep2Data$zerozero
rep2Data = rep2Data[totalCells >= minimumCells]

rep2Data$fdr = p.adjust(rep2Data$pval, method = "BH")
rep2Data = unique(rep2Data[,.(gene,enh1,enh2,corr,pval)])

mergedData = merge(rep3Sign, rep2Data, by = c("gene","enh1","enh2"))

options(scipen = 5)

t = cor.test(mergedData$corr.x,mergedData$corr.y, method = "spearman")
if(t$p.value == 0){ pval = 2.2e-308 }
text = paste("Spearman R = ",round(t$estimate,2), " P-value <",pval, sep = "")
ggplot(mergedData, aes(x = corr.x, y = corr.y) ) +
  geom_bin2d(bins = 50) +
  scale_fill_gradient(high = "#3f007d", low = "#efedf5") +
  geom_smooth(method = "lm") +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.05, vjust = 1.5, size = 5, fontface = "bold"  ) +
  xlab("Main rep correlation (24844 cells)") +
  ylab("Other rep correlation (2788 cells)") +
  theme_minimal() +
  theme(text = element_text(size = 18), plot.title = element_text(hjust = 0.5), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1  )

