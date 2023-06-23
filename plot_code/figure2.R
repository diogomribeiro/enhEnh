#17-May-2023 Diogo Ribeiro @ UNIL

library(data.table)
library(ggplot2)
library(patchwork)

pdf("~/Figure2_new.pdf",18,16)

##############
# Panel A
##############

exampleGene = "ENSG00000100439"
data = fread("~/git/enhEnh/source_data/gene_enhancer_combinations.out.gz", header = T, sep = "\t") 

example = data[gene == exampleGene]
example = example[order(NBcell)]
example$combID = seq(1,nrow(example))

exampleDT = data.table()
for (i in seq(1,nrow(example)) ){
  row = example[i,]
  enhs = data.table(unlist(lapply(example[i,]$enhancerList, function(x) unlist(strsplit(x,"[,]")))))  
  for (j in seq(1,nrow(enhs)) ){
    enh = enhs[j,]$V1
    enh = gsub(" ","", enh)
    d = data.table(gene = row$gene, enh = enh, cells = row$NBcell, combID = row$combID)
    exampleDT = rbind(exampleDT, d)
  }
}

exampleDT

enhs = data.table(unique(exampleDT$enh))
enhs = enhs[order(V1)]
enhs$enhID = seq(1,nrow(enhs))

mergedData = merge(exampleDT, enhs, by.x = "enh", by.y = "V1")
mergedData$active = 1
mergedData$combID = as.factor(mergedData$combID)

order = data.table(table(mergedData$combID))
order = order[order(-N)]
order$order = seq(1,nrow(order))

mergedData = merge(mergedData,order, by.x = "combID", by.y = "V1")
mergedData = mergedData[order(order)]

## Add missing manually
mergedData = rbind(mergedData, data.table(combID = "2",enh = "", gene = "", cells = NA, enhID = 2, active = 0, N = 2, order = 2))
mergedData = rbind(mergedData, data.table(combID = "4",enh = "", gene = "", cells = NA, enhID = 3, active = 0, N = 2, order = 3))
mergedData = rbind(mergedData, data.table(combID = "5",enh = "", gene = "", cells = NA, enhID = 1, active = 0, N = 2, order = 4))
mergedData = rbind(mergedData, data.table(combID = "1",enh = "", gene = "", cells = NA, enhID = 2, active = 0, N = 2, order = 5))
mergedData = rbind(mergedData, data.table(combID = "1",enh = "", gene = "", cells = NA, enhID = 3, active = 0, N = 2, order = 5))
mergedData = rbind(mergedData, data.table(combID = "6",enh = "", gene = "", cells = NA, enhID = 1, active = 0, N = 2, order = 6))
mergedData = rbind(mergedData, data.table(combID = "6",enh = "", gene = "", cells = NA, enhID = 2, active = 0, N = 2, order = 6))
mergedData = rbind(mergedData, data.table(combID = "7",enh = "", gene = "", cells = NA, enhID = 1, active = 0, N = 2, order = 7))
mergedData = rbind(mergedData, data.table(combID = "7",enh = "", gene = "", cells = NA, enhID = 3, active = 0, N = 2, order = 7))

mergedData$active = as.factor(mergedData$active)

cellDT = unique(mergedData[,.(cells,order)][!is.na(cells)])

g1 = ggplot(mergedData, aes(x = enhID, y = as.factor(-order), fill = cells, color = active)) + 
  geom_tile( size = 1, width = 0.7, height = 0.5) +
  geom_hline(yintercept = 1.5, color = "grey") +
  geom_hline(yintercept = 2.5, color = "grey") +
  geom_hline(yintercept = 3.5, color = "grey") +
  geom_hline(yintercept = 4.5, color = "grey") +
  geom_hline(yintercept = 5.5, color = "grey") +
  geom_hline(yintercept = 6.5, color = "grey") +
  geom_text(data = cellDT, aes(y = as.factor(-order), x = 3.5, label = cells, color = "1"), size = 5.5 ) +
  annotate(geom = "text", label = "cells", x = 3.5, y = 7.4, size = 5 ) +
  scale_x_discrete(limits=c("chr14:22851-22854","chr14:23016-23026","chr14:23388-23388")) +
  scale_color_manual(values = c("grey","black")) +
  scale_fill_gradient(low = "#dadaeb", high = "#6a51a3", na.value = "white", trans = "log10") +
  labs(tag = "a") +
  xlab("Enhancers") +
  ylab("Enhancer combinations") +
  theme_linedraw() + 
  theme(text = element_text(size=25), axis.text.x = element_text(size = 16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2), aspect.ratio = 1)  


# ### Real coordinates
# 
# pdf("~/Figure2_detail.pdf",16,1.8)
# 
# geneModel = fread("/home/dribeiro/EnhEnhPaper/data/genes.bed")
# geneModel$gene = data.table(unlist(lapply(geneModel$V4, function(x) unlist(strsplit(x,"[.]"))[1])))$V1
# geneModel = geneModel[gene == exampleGene]
# enhs$chr = data.table(unlist(lapply(enhs$V1, function(x) unlist(strsplit(x,"[_]"))[1])))$V1
# enhs$enh_start = as.numeric(data.table(unlist(lapply(enhs$V1, function(x) unlist(strsplit(x,"[_]"))[2])))$V1)
# enhs$enh_end = as.numeric(data.table(unlist(lapply(enhs$V1, function(x) unlist(strsplit(x,"[_]"))[3])))$V1)
# 
# enhs$st = enhs$enh_start/1000
# enhs$en = enhs$enh_end/1000
# geneModel$st = geneModel$V2/1000
# geneModel$en = geneModel$V3/1000
# 
# ggplot() +
#   geom_rect( data = enhs, aes(xmin = st-1, xmax = en+1, ymin = 1, ymax = 1.1), fill = "grey", color = "black", size = 0.5) +
#   geom_rect( data = geneModel, aes(xmin = st-1, xmax = en+1, ymin = 1, ymax = 1.1), fill = "black", size = 0.5) +
#   ylim(c(0.98,1.12)) +
#   xlab("Genomic coordinates (Kb)") +
#   theme_minimal() +
#   theme(text = element_text(size = 42), plot.title = element_text(hjust = 0.5), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
#         panel.background = element_rect(colour = "black", fill = "white", size = 1), axis.text.y = element_blank())
# 
# dev.off()

##############
# Panel B
##############

mergedData = fread("~/git/enhEnh/source_data/data_figure_2b.tsv.gz", header = T, sep = "\t") 

t = cor.test(mergedData$sign_enh, mergedData$ratio, method = "spearman")
t$estimate
t$p.value

mergedData$NEnh = factor(mergedData$NEnh, levels=c("1", "2", "3", "4", "5", "6", "7","8-10", "10-15", ">15"))

g2 = ggplot(mergedData[!is.na(ratio)][NEnh != 1], aes(x=NEnh, y=ratio)) + 
  geom_boxplot( fill="#69b3a2", size = 1, width = 0.7, alpha = 0.7, outlier.shape = 1) +
  geom_text(aes(label=paste(..count.., sep="")), y=-7.5, stat='count', colour="black", size=6)+
  stat_summary(fun=mean, geom="point", size=3, color="#969696", shape = 19)+
  stat_summary(fun=mean, geom="text", size=6, color="black",vjust = 1.5, aes(label= paste( round(..y.., digits = 1))), fontface = "bold")+
  ylim (c(-5,100))+ 
  labs(x="Number of enhancers", y = "% combinations observed")+
  labs(tag = "b") +
  theme_linedraw() + 
  theme(text = element_text(size=25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2), aspect.ratio = 1)  

##############
# Panel C
##############

mergedData = fread("~/git/enhEnh/source_data/data_figure_2c.tsv.gz", header = T, sep = "\t") 

mergedData$NEnh = factor(mergedData$NEnh, levels=c("1", "2", "3", "4", "5", "6", "7","8-10", "10-15", ">15"))

g3 = ggplot(mergedData[!is.na(ratio)], aes(x=NEnh, y=ratio)) + 
  geom_boxplot( fill="#69b3a2", size = 1,  width = 0.7, alpha = 0.7, outlier.shape = 1) +
  geom_text(aes(label=paste(..count.., sep = "") ), y=-7.5, stat='count', colour="black", size=6)+
  stat_summary(fun=mean, geom="point", size=3, color="#969696", shape = 19)+
  stat_summary(fun=mean, geom="text", size=6, color="black",vjust = 1.5, aes(label= paste( round(..y.., digits = 1))), fontface = "bold")+
  ylim (c(-5,100))+ 
  labs(x="Number of enhancers", y = "% enhancer pairs significant")+
  labs(tag = "c") +
  theme_linedraw() + 
  theme(text = element_text(size=25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2), aspect.ratio = 1)  

t = cor.test(mergedData$ratio,mergedData$sign_enh,method = "spearman")
t
t$p.value

##############
# Panel D
##############

mergedData = fread("~/git/enhEnh/source_data/data_figure_2d.tsv.gz", header = T, sep = "\t")

mergedData$n_genes = factor(mergedData$n_genes, levels=c("1", "2", "3", "4", "5", "6", "7","8-10", "10-15", ">15"))

g4 = ggplot(mergedData, aes(x=n_genes, y=ratio)) + 
  geom_boxplot( fill="#bebada", size = 1, width = 0.7, alpha = 0.7, outlier.shape = 1) +
  geom_text(aes(label=paste(..count.., sep = "") ), y=-7.5, stat='count', colour="black", size=6)+
  stat_summary(fun=mean, geom="point", size=3, color="#969696", shape = 19)+
  stat_summary(fun=mean, geom="text", size=6, color="black",vjust = 1.5, aes(label= paste( round(..y.., digits = 1))), fontface = "bold")+
  ylim (c(-5,100))+ 
  labs(x="Number of genes", y = "% genes significant")+
  labs(tag = "d") +
  theme_linedraw() + 
  theme(text = element_text(size=24), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2),  aspect.ratio = 1)  

t = cor.test(mergedData$possible,mergedData$ratio, method = "spearman")
t
t$p.value

######################

g1 + g2 + g3 + g4
#  inset_element(g2, left = 0, bottom = 0.9, right = 1, top = 1.1, on_top = TRUE, align_to = "full")

dev.off()

