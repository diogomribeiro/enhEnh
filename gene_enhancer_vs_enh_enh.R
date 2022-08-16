# 10-Aug-2022 Diogo Ribeiro, Chaymae Ziyani @ UNIL
# Script to plot the gene-enhancer significance in the x axis and the enh-enh significant in the y axis

library(data.table)
library(ggplot2)
options(digits = 3)

corrCutoff = 0.05
fdrCutoff = 0.05
minTotalCells = 100

signGeneEnh = fread("/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/enh_enh_paper/data/gene_enhancer_significant.tsv", header = T, sep = "\t") 
signGeneEnh$tag = paste(signGeneEnh$gene,"|chr",signGeneEnh$peak,sep="")

allEnhEnh = fread("/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/enh_enh_paper/data/enh_enh_correlation_all.tsv.gz", header = T, sep = "\t") 
allEnhEnh$totalCells = allEnhEnh$oneOne + allEnhEnh$oneZero + allEnhEnh$zeroOne + allEnhEnh$zerozero
allEnhEnh = allEnhEnh[totalCells >= minTotalCells]
allEnhEnh$fdr = p.adjust(allEnhEnh$pval,method = "BH")

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
dt = data.table(prop = c(v2,v1,v4,v3), grp = c("sign", "sign", "nonsign", "nonsign"), enh_enh_grp = c("nonsign", "sign", "nonsign", "sign"))

##### Plot
group=c("1","2","3","4")
ggplot(dt, aes(x=grp,y=prop*100, fill = group)) + 
  geom_bar( stat = "identity", width = .8, color = "black")+
  ylab("% enh-enh significant") +
  xlab("Gene-enhancer significance") +
  #guides(fill="none")+
  scale_fill_manual( values =  c("#ccece6","#66C2A5","#d9d9d9","#969696"), labels = c("significant","non-signicant", "significant","non-signicant"  ) ) +
  #  scale_fill_brewer(palette = "Paired" )+
  geom_text(aes(label =paste(format(prop*100, nsmall=0), "%")), position = position_stack(vjust= 0.5),
            colour = "black", size = 5)+
  theme_linedraw() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=24),
        legend.text=element_text(size=20), legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))  
