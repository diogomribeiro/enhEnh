# 10-Aug-2022 Diogo Ribeiro, Chaymae Ziyani @ UNIL
# Script to plot the gene-enhancer significance in the x axis and the enh-enh significant in the y axis

library(data.table)
library(ggplot2)
options(digits = 3)

corrCutoff = 0.05
fdrCutoff = 0.05
minTotalCells = 100

signGeneEnh = fread("/home/dribeiro/EnhEnhPaper/data/gene_enhancer_significant.tsv", header = T, sep = "\t") 
signGeneEnh$tag = paste(signGeneEnh$gene,"|chr",signGeneEnh$peak,sep="")

allEnhEnh = fread("/home/dribeiro/EnhEnhPaper/data/enh_enh_correlation_all.tsv.gz", header = T, sep = "\t") 
allEnhEnh$totalCells = allEnhEnh$oneOne + allEnhEnh$oneZero + allEnhEnh$zeroOne + allEnhEnh$zerozero
allEnhEnh = allEnhEnh[totalCells >= minTotalCells]
allEnhEnh$fdr = p.adjust(allEnhEnh$pval,method = "BH")

# allEnhEnh$significant = "no"
# allEnhEnh[fdr < 0.05][corr > 0.05]$significant = "yes"
# write.table(allEnhEnh,"~/EnhEnhPaper/data/enh_enh_all_supp_table.tsv.gz",row.names=F,quote=F,sep="\t")

allEnhEnh$tag1 = paste(allEnhEnh$gene,"|",allEnhEnh$enh1,sep="")
allEnhEnh$tag2 = paste(allEnhEnh$gene,"|",allEnhEnh$enh2,sep="")

paste(round(nrow(allEnhEnh[corr > corrCutoff][fdr < fdrCutoff])/nrow(allEnhEnh)*100,digits = 3),"% of all enh-enh are significant",sep="")

signEnhEnh = allEnhEnh[tag1 %in% signGeneEnh$tag & tag2 %in% signGeneEnh$tag]
nonsignEnhEnh = allEnhEnh[!tag1 %in% signGeneEnh$tag & !tag2 %in% signGeneEnh$tag]


##### calculate the propotion
v1 = nrow(signEnhEnh[corr > corrCutoff][fdr < fdrCutoff]) /nrow(signEnhEnh)
v2 = 1 - v1
v3 = nrow(nonsignEnhEnh[corr > corrCutoff][fdr < fdrCutoff]) / nrow(nonsignEnhEnh)
v4 = 1 - v3

##### table containing the 4 values
dt = data.table(prop = c(v2,v1,v4,v3), grp = c("significant", "significant", "non-significant", "non-significant"), enh_enh_grp = c("non-significant", "significant", "non-significant", "significant"))

m = matrix(c(nrow(signEnhEnh[corr > corrCutoff][fdr < fdrCutoff]), nrow(signEnhEnh), nrow(nonsignEnhEnh[corr > corrCutoff][fdr < fdrCutoff]), nrow(nonsignEnhEnh)),nrow=2)
t = fisher.test(m)
t
t$p.value

##### Plot
`enh-enh assoc.`=c("1","2","3","4")
ggplot(dt, aes(x=grp,y=prop*100, fill = `enh-enh assoc.`)) + 
  geom_bar( stat = "identity", width = .8, color = "black", size = 1)+
  ylab("% enh-enh significant") +
  xlab("Gene-enhancer significance") +
  #guides(fill="none")+
  scale_fill_manual( values =  c("#ccece6","#66999B","#d9d9d9","#969696"), labels = c("non-significant","significant", "non-significant","significant"  ) ) +
  #  scale_fill_brewer(palette = "Paired" )+
  geom_text(aes(label =paste(format(prop*100, nsmall=0), "%")), position = position_stack(vjust= 0.5),
            colour = "black", size = 6)+
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=22),
        legend.text=element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))  

