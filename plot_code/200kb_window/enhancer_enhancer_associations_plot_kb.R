# 09-Aug-2022 Diogo Ribeiro, Chaymae Ziyani @ UNIL
# Script to parse enh-enh correlations

library(data.table)
library(ggplot2)

data = fread("../../source_data/200kb_window/200kb_enh_enh_correlation.tsv.gz", header = T, sep = "\t")
corrCutoff = 0.05
fdrCutoff = 0.05
minTotalCells = 100

initialGenes = length(unique(data$gene))

# Apply filter for minimum number of cells expressing the gene
data$totalCells = data$oneOne + data$oneZero + data$zeroOne + data$zerozero
data = data[totalCells >= minTotalCells]
length(unique(data$gene))
data[corr < 0]

data$fdr = p.adjust(data$pval, method = "BH")
paste(nrow(data),"total tested")
paste(nrow(data[fdr < fdrCutoff]),"kept after FDR filter")
paste(nrow(data[abs(corr) > corrCutoff]),"kept after Corr filter")

data$significant = "no"
data[abs(corr) > corrCutoff][fdr < fdrCutoff]$significant = "yes"
sign = data[abs(corr) > corrCutoff][fdr < fdrCutoff]
paste(nrow(sign),"kept after Corr & FDR filter")
paste(nrow(sign[corr < 0]),"negatively correlated")

paste(round(nrow(sign)*100/nrow(data),2),"% associations significant")

paste(initialGenes,"initial genes")
paste(length(unique(data$gene)),"total genes")
paste(length(unique(sign$gene)),"genes with significant enh-enh associations")

geneCells = unique(data[,.(gene,totalCells)])
summary(geneCells$totalCells)

## Correlation distribution
ggplot(data, aes(x=corr, fill = significant) ) + 
  geom_histogram( color = "black", binwidth = 0.01, position = "stack", size = 0.5) +
  scale_fill_manual(values = c("#B3AF8F","#66999B")) +
  geom_vline(xintercept = 0.05, color = "#525252", linetype = "dashed", size = 1)+
  xlab("Enhancer-enhancer correlation") +
  ylab("Frequency") +
  theme_linedraw() + 
  theme(text = element_text(size=45), # big screen
        legend.text=element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2), aspect.ratio = 1)  

dt = data.table(significant = c("1:yes","2:no"), values = c(nrow(data[significant == "yes"]), nrow(data[significant == "no"])))

dataCorr0.1 = data
dataCorr0.1$significant = "no"
dataCorr0.1[abs(corr) > 0.1][fdr < fdrCutoff]$significant = "yes"
dataBH = data
dataBH$significant = "no"
dataBH[fdr < fdrCutoff]$significant = "yes"
dataBon = data
dataBon$fdr = p.adjust(data$pval, method = "bonferroni")
dataBon$significant = "no"
dataBon[fdr < fdrCutoff]$significant = "yes"

table(data$significant)


dt$perc = dt$values / nrow(data) * 100
dt$y = dt$values
dt[significant == "no"]$y = nrow(data)

ggplot(dt, aes(x= "", y = values, fill = significant ) ) + 
  geom_bar( stat = "identity", position = "stack", color = "black", size = 1) +
  scale_fill_manual(values = c("#66999B","#B3AF8F")) +
  geom_vline(xintercept = 0.05, color = "grey", linetype = "dashed")+
  geom_text(aes(y = y, label = paste0("N=",values)), size = 8, color = "black", vjust = 2, fontface = "bold") + 
  geom_text(aes(y = y, label = paste0(round(perc,1),"%")), size = 8, color = "black", vjust = 3.5, fontface = "bold") +
  xlab("") +
  ylab("Enhancer pairs") +
  theme_linedraw() + 
  theme(text = element_text(size=30), legend.position = "none",
        legend.text=element_text(size=20), 
        axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0))  
