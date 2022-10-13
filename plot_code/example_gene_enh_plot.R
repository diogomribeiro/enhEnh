#12-Aug-2022 Diogo Ribeiro @ UNIL
# Script to plot combinations of enhancers and cell number

library(data.table)
library(ggplot2)

exampleGene = "ENSG00000100439"
data = fread("../source_data/gene_enhancer_combinations.out.gz", header = T, sep = "\t") 

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
# mergedData$enhID = as.factor(mergedData$enhID)

cellDT = unique(mergedData[,.(cells,order)][!is.na(cells)])

# ggplot(mergedData, aes(xmin = enhID-0.5, xmax = enhID+0.5, ymin = combID, ymax = combID+1, fill = log10(cells)) ) + 
ggplot(mergedData, aes(x = enhID, y = as.factor(-order), fill = cells, color = active)) + 
  geom_tile( size = 1, width = 0.7, height = 0.5) +
  geom_hline(yintercept = 1.5, color = "grey") +
  geom_hline(yintercept = 2.5, color = "grey") +
  geom_hline(yintercept = 3.5, color = "grey") +
  geom_hline(yintercept = 4.5, color = "grey") +
  geom_hline(yintercept = 5.5, color = "grey") +
  geom_hline(yintercept = 6.5, color = "grey") +
  geom_text(data = cellDT, aes(y = as.factor(-order), x = 3.6, label = cells, color = "1"), size = 5 ) +
  theme_minimal() +
  xlim(c(0.5,3.6)) +
  xlab("Enhancers") +
  ylab("Enhancer combinations") +
  scale_color_manual(values = c("grey","black")) +
  scale_fill_gradient(low = "#dadaeb", high = "#6a51a3", na.value = "white", trans = "log10") +
  theme(text = element_text(size = 18), plot.title = element_text(hjust = 0.5), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), axis.text.y = element_blank()  )
#aspect.ratio = 1, 


### Real coordinates
geneModel = fread("../source_data/genes.bed.gz")
geneModel$gene = data.table(unlist(lapply(geneModel$V4, function(x) unlist(strsplit(x,"[.]"))[1])))$V1
geneModel = geneModel[gene == exampleGene]
enhs$chr = data.table(unlist(lapply(enhs$V1, function(x) unlist(strsplit(x,"[_]"))[1])))$V1
enhs$enh_start = as.numeric(data.table(unlist(lapply(enhs$V1, function(x) unlist(strsplit(x,"[_]"))[2])))$V1)
enhs$enh_end = as.numeric(data.table(unlist(lapply(enhs$V1, function(x) unlist(strsplit(x,"[_]"))[3])))$V1)

enhs$st = enhs$enh_start/1000
enhs$en = enhs$enh_end/1000
geneModel$st = geneModel$V2/1000
geneModel$en = geneModel$V3/1000

ggplot() + 
  geom_rect( data = enhs, aes(xmin = st-1, xmax = en+1, ymin = 1, ymax = 1.1), fill = "grey", color = "black", size = 0.5) +
  geom_rect( data = geneModel, aes(xmin = st-1, xmax = en+1, ymin = 1, ymax = 1.1), fill = "black", size = 0.5) +
  ylim(c(0.98,1.12)) +
  xlab("Chr14 genomic coordiantes (Kb)") +
  theme_minimal() +
  theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), axis.text.y = element_blank())

