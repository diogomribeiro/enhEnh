# 09-Aug-2022 Diogo Ribeiro, Chaymae Ziyani @ UNIL
# Script to parse enh-enh correlations

library(data.table)
library(ggplot2)

# data = fread("~/EnhEnhPaper/data/enh_enh_correlation.tsv", header = T, sep = "\t")
data = fread("~/EnhEnhPaper/data/abc/enh_enh_correlation.tsv", header = T, sep = "\t")
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

# write.table(data,"/home/dribeiro/EnhEnhPaper/data/enh_enh_supp_table.tsv",row.names=F,quote=F,sep="\t")

# nrow(data[fdr < 0.05])/nrow(data)
# nrow(data[significant == "yes"])/nrow(data)
# nrow(data[fdr < 0.05][corr > 0.05])/nrow(data)
# nrow(data[fdr < 0.05][corr > 0.1])/nrow(data)
# data$bon = p.adjust(data$pval, method = "bonferroni")
# data$significant = "no"
# data[abs(corr) > corrCutoff][bon < fdrCutoff]$significant = "yes"
# sign = data[abs(corr) > corrCutoff][bon < fdrCutoff]
# nrow(data[bon < 0.05])/nrow(data)

## Correlation distribution
ggplot(data, aes(x=corr, fill = significant) ) + 
  geom_histogram( color = "black", binwidth = 0.01, position = "stack", size = 0.5) +
  scale_fill_manual(values = c("#B3AF8F","#66999B")) +
  geom_vline(xintercept = 0.05, color = "#525252", linetype = "dashed", size = 1)+
  xlab("Enhancer-enhancer correlation") +
  ylab("Frequency") +
  theme_linedraw() + 
#  theme(text = element_text(size=20), # small screen
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

dt = data.table(cutoff = c("BH 5% &\n Corr > 0.05","BH 5% &\n Corr > 0.05","BH 5%","BH 5%","Bonf. 0.05","Bonf. 0.05","BH 5% &\n Corr > 0.1", "BH 5% &\n Corr > 0.1"),
                # significant = c("1:yes","2:no","1:yes","2:no","1:yes","2:no","1:yes","2:no"),
                significant = c("yes","no","yes","no","yes","no","yes","no"),
                values = c(nrow(data[significant == "yes"]), nrow(data[significant == "no"]), 
                           nrow(dataBH[significant == "yes"]), nrow(dataBH[significant == "no"]),
                           nrow(dataBon[significant == "yes"]), nrow(dataBon[significant == "no"]),
                           nrow(dataCorr0.1[significant == "yes"]), nrow(dataCorr0.1[significant == "no"]))
                )

dt$perc = dt$values / nrow(data) * 100
dt$y = dt$values
dt[significant == "no"]$y = nrow(data)

ggplot(dt, aes(x= cutoff, y = values, fill = significant ) ) + 
  geom_bar( stat = "identity", position = "stack", color = "black", size = 2) +
  # scale_fill_brewer(palette = "Set2")+
  scale_fill_manual(values = c("#B3AF8F","#66999B")) +
  geom_vline(xintercept = 0.05, color = "grey", linetype = "dashed")+
  geom_text(aes(y = y, label = paste("N =",values)), size = 10, color = "black", vjust = 1.5, fontface = "bold") + 
  geom_text(aes(y = y, label = paste0(round(perc,1),"%")), size = 10, color = "black", vjust = 3, fontface = "bold") +
  # stat_summary(fun=mean, geom="text", size=6, color="black", vjust = 2, aes(label= paste("N =",..y..)) )+
  xlab("Cutoff") +
  ylab("Enhancer pairs") +
  theme_linedraw() + 
  # theme(text = element_text(size=20),
  theme(text = element_text(size=42), legend.position = "none",
        legend.text=element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2), aspect.ratio = 1)  

ggplot(dt[cutoff == "BH 5% &\n Corr > 0.05"], aes(x= cutoff, y = values, fill = significant ) ) + 
  # geom_bar( stat = "identity", position = "dodge", color = "black", alpha = 0.9) +
  geom_bar( stat = "identity", position = "stack", color = "black", size = 1) +
  scale_fill_manual(values = c("#B3AF8F","#66999B")) +
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
  scale_fill_manual(values = c("#B3AF8F","#66999B")) +
  # geom_vline(xintercept = 1000, color = "grey", linetype = "dashed")+
  xlab("Absolute distance (kb)") +
  theme_linedraw() + 
  theme(text = element_text(size=20),
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

write.table(sign,"~/EnhEnhPaper/data/significant_enh_enh.tsv",quote=F,sep="\t",row.names=F)
