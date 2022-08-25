#15-Aug-2022 Diogo Ribeiro @ UNIL
# Script to measure enhancer-enhancer associations across genes

enhEnhData = fread("~/EnhEnhPaper/data/enh_enh_correlation.tsv", header = T, sep = "\t") 
signData = fread("~/EnhEnhPaper/data//significant_enh_enh.tsv", header = T, sep = "\t") 

length(unique(enhEnhData$tag))

enhEnhData$tag = paste(enhEnhData$enh1,enhEnhData$enh2)
possible = data.table(table(enhEnhData$tag))

signData$tag = paste(signData$enh1,signData$enh2)
significant = data.table(table(signData$tag))

mergedData = merge(possible, significant, by = "V1", all.x = T)
colnames(mergedData) = c("enh_pair","possible","significant")
mergedData[is.na(significant)]$significant = 0

#plot(mergedData$possible,mergedData$significant)

mergedData$ratio = mergedData$significant * 100 / mergedData$possible
mergedData$n_genes = as.factor(mergedData$possible)
mergedData[possible >= 8][possible <= 10]$n_genes = "8-10"
mergedData[possible >= 10][possible <= 15]$n_genes = "10-15"
mergedData[possible > 15]$n_genes = ">15"

ggplot(mergedData, aes(x=n_genes, y=ratio)) + 
  geom_boxplot( fill="#bebada", alpha = 0.7) +
  labs(x="Number of genes", y = "% genes significant")+
  ylim (c(-5,100))+ 
  geom_text(aes(label=paste("N=",..count.., sep = "") ), y=-5, stat='count', colour="black", size=4)+
  stat_summary(fun=mean, geom="point", size=2, color="#969696")+
  stat_summary(fun=mean, geom="text", size=5, color="black",vjust = 1.5, aes(label= paste( round(..y.., digits = 1))))+
  theme_linedraw() + 
  theme(text = element_text(size=22),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),  aspect.ratio = 1)  

t = cor.test(mergedData$possible,mergedData$ratio, method = "spearman")
t
t$p.value


mergedData[significant > 10]
mergedData[enh_pair == "chr10_103660600_103665400 chr10_103760000_103764400"]
enhEnhData[tag == "chr10_103660600_103665400 chr10_103760000_103764400"]


enhEnhData[gene == "ENSG00000187837"]
enhEnhData[tag == "chr6_26104800_26105400 chr6_26189200_26191000"]
mergedData[enh_pair == "chr6_26104800_26105400 chr6_26189200_26191000"]
