# 09-Aug-2022 Diogo Ribeiro, Chaymae Ziyani @ UNIL
# Script to parse enh-enh correlations

library(data.table)
library(ggplot2)

data = fread("/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/enh_enh_paper/data/enh_enh_correlation.tsv", header = T, sep = "\t") 
corrCutoff = 0.05
fdrCutoff = 0.05
minTotalCells = 100

initialGenes = length(unique(data$gene))

# Apply filter for minimum number of cells expressing the gene
data$totalCells = data$oneOne + data$oneZero + data$zeroOne + data$zerozero
data = data[totalCells >= minTotalCells]

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
  geom_histogram( color = "black", size = 0.1, binwidth = 0.01, position = "stack") +
  scale_fill_brewer(palette = "Set2")+
  geom_vline(xintercept = 0.05, color = "grey", linetype = "dashed")+
  xlab("Enhancer-enhancer correlation") +
  ylab("Frequency") +
  theme_linedraw() + 
  theme(text = element_text(size=24),
        legend.text=element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio = 1)  

## Distance comparison
data$enh1_start = as.numeric(data.table(unlist(lapply(data$enh1, function(x) unlist(strsplit(x,"[_]"))[2])))$V1)
data$enh1_end = as.numeric(data.table(unlist(lapply(data$enh1, function(x) unlist(strsplit(x,"[_]"))[3])))$V1)
data$enh2_start = as.numeric(data.table(unlist(lapply(data$enh2, function(x) unlist(strsplit(x,"[_]"))[2])))$V1)
data$enh2_end = as.numeric(data.table(unlist(lapply(data$enh2, function(x) unlist(strsplit(x,"[_]"))[3])))$V1)

data$enh1_midpoint = (data$enh1_start + data$enh1_end) / 2
data$enh2_midpoint = (data$enh2_start + data$enh2_end) / 2

data$distance = abs(data$enh1_midpoint - data$enh2_midpoint)

summary(data[significant == "yes"]$distance)
summary(data[significant == "no"]$distance)
wilcox.test(data[significant == "yes"]$distance,data[significant == "no"]$distance)

ggplot(data, aes(x=distance/1000, fill = significant) ) + 
  geom_density(alpha = 0.5)+
  scale_fill_brewer(palette = "Set2")+
  # geom_vline(xintercept = 1000, color = "grey", linetype = "dashed")+
  xlab("Absolute distance (kb)") +
  theme_linedraw() + 
  theme(text = element_text(size=24),
        legend.text=element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), aspect.ratio = 1)  

nrow(data[distance > 1000000][significant == "yes"])/nrow(sign)
nrow(data[distance > 1000000][significant == "no"])/nrow(data[significant == "no"])

## Logistic regression
# data$sign = 0
# data[significant == "yes"]$sign = 1
# g = glm(sign ~ distance, data = data, family = binomial)
# summary(g)


write.table(sign,"/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/enh_enh_paper/data/significant_enh_enh.tsv",quote=F,sep="\t",row.names=F)
