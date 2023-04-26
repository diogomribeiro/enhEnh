#25-Jan-2023 Diogo Ribeiro @ UNIL
# Exploring partial correlation results

options(digits=3)

library(ggplot2)
library(data.table)

dt = fread("../../source_data/partial_correlation/share_seq_partial_corr_all.out.gz")

#################
# partial vs simple correlation (24k cells)
#################

test = cor.test(dt$enh1_enh2_pcor,dt$enh1_enh2_scor, method = "spearman")
pval <- ifelse(test$p.value == 0, 2.2e-308, test$p.value)
text = paste("Spearman R = ",round(test$estimate,3),"p-value <", format(pval, scientific = TRUE))
ggplot(dt, aes(x = enh1_enh2_pcor, y = enh1_enh2_scor)) + 
  geom_point(alpha = 0.5) + 
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = "grey", linetype = "dashed") +
  annotate(geom = "text", x = Inf, y = Inf, label = text, size = 8, hjust = 1.2, vjust = 2) +
  xlab("Partial correlation") +
  ylab("Normal correlation") +
  theme_minimal() + 
  theme(text = element_text(size=30), aspect.ratio = 1,
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2.5))

test = cor.test(-log10(dt$enh1_enh2_scor_pval),-log10(dt$enh1_enh2_pcor_pval), method = "spearman")
pval <- ifelse(test$p.value == 0, 2.2e-308, test$p.value)
text = paste("Spearman R = ",round(test$estimate,3),"p-value <", format(pval, scientific = TRUE))
ggplot(dt, aes(x = -log10(enh1_enh2_scor_pval), y = -log10(enh1_enh2_pcor_pval)) ) + 
  geom_point(alpha = 0.5) + 
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = "grey", linetype = "dashed") +
  annotate(geom = "text", x = Inf, y = Inf, label = text, size = 8, hjust = 1.2, vjust = 2) +
  xlab("Partial corr. p-values (-log10)") +
  ylab("Normal corr. p-values (-log10)") +
  theme_minimal() + 
  theme(text = element_text(size=30), aspect.ratio = 1,
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2.5))

##############
# partial vs previous approach
##############

enhEnhData = fread("../source_data/enh_enh_correlation.tsv")
minTotalCells = 100
enhEnhData$totalCells = enhEnhData$oneOne + enhEnhData$oneZero + enhEnhData$zeroOne + enhEnhData$zerozero
enhEnhData = enhEnhData[totalCells >= minTotalCells]
enhEnhData$tag = paste(enhEnhData$gene, enhEnhData$enh1, enhEnhData$enh2, sep = "|")
signEnhEnhData = fread("/home/dribeiro/EnhEnhPaper/data/significant_enh_enh.tsv")
signEnhEnhData$tag = paste(signEnhEnhData$gene, signEnhEnhData$enh1, signEnhEnhData$enh2, sep = "|")
enhEnhData$significant = "no"
enhEnhData[tag %in% signEnhEnhData$tag]$significant = "yes"

dt$fdr_partial = p.adjust(dt$enh1_enh2_pcor_pval, method = "BH")
dt$fdr_simple = p.adjust(dt$enh1_enh2_scor_pval, method = "BH")
dt$tag = paste(dt$gene, dt$enhancer1, dt$enhancer2, sep = "|")

mergedData = merge(dt, enhEnhData, by = "tag", all.x = TRUE)
mergedData = mergedData[!is.na(significant)]

m = mergedData[,.(tag,corr, pval, enh1_enh2_pcor, enh1_enh2_pcor_pval, enh1_enh2_scor, enh1_enh2_scor_pval)]

test = cor.test(m$corr,m$enh1_enh2_pcor, method = "spearman")
pval <- ifelse(test$p.value == 0, 2.2e-308, test$p.value)
text = paste("Spearman R = ",round(test$estimate,3),"p-value <", format(pval, scientific = TRUE))
ggplot(m, aes(x = enh1_enh2_pcor, y = corr) ) + 
  geom_point(alpha = 0.5) + 
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = "grey", linetype = "dashed") +
  annotate(geom = "text", x = Inf, y = Inf, label = text, size = 8, hjust = 1.2, vjust = 2) +
  xlab("Partial correlation (24k cells)") +
  ylab("Normal correlation (subset cells)") +
  theme_minimal() + 
  theme(text = element_text(size=30), aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2.5))

test = cor.test(-log10(m$pval),-log10(m$enh1_enh2_pcor_pval), method = "spearman")
pval <- ifelse(test$p.value == 0, 2.2e-308, test$p.value)
text = paste("Spearman R = ",round(test$estimate,3),"p-value <", format(pval, scientific = TRUE))
ggplot(m, aes(x = -log10(enh1_enh2_pcor_pval), y = -log10(pval)) ) + 
  geom_point(alpha = 0.5) + 
  geom_abline(slope = 1, intercept = 0, size = 1.5, color = "grey", linetype = "dashed") +
  annotate(geom = "text", x = Inf, y = Inf, label = text, size = 8, hjust = 1.2, vjust = 2) +
  xlab("Partial corr. p-value (-log10)") +
  ylab("Normal corr. p-value (-log10)") +
  theme_minimal() + 
  theme(text = element_text(size=30), aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2.5))


nrow(mergedData[significant == "yes"][fdr_partial > 0.05])
nrow(mergedData[significant == "yes"][fdr_partial < 0.05])

summary(mergedData[significant == "yes"]$enh1_enh2_pcor)
summary(mergedData[significant == "no"]$enh1_enh2_pcor)
ggplot(mergedData, aes(enh1_enh2_pcor, fill = significant)) + 
  geom_density(alpha = 0.5) + 
  scale_fill_brewer(palette = "Set2") +
  xlab("Partial correlation") +
  # xlim(c(-0.05,0.2)) +
  theme_minimal() + 
  theme(text = element_text(size=30), legend.position = c(0.8, 0.8),  aspect.ratio = 1,
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2.5))
w = wilcox.test(mergedData[significant == "yes"]$enh1_enh2_pcor,mergedData[significant == "no"]$enh1_enh2_pcor)
w$p.value

#########
# Number significant

bk = data.table(table(mergedData$significant))
d1 = data.table(table(mergedData[enh1_enh2_pcor >= 0.05]$significant))
d2 = data.table(d1$N * 100 / bk$N)
d2$dataset = c("Corr<0.05","Corr>=0.05")

m = matrix(c(d1$N[2],bk$N[2],d1$N[1],bk$N[1]),nrow=2)
f = fisher.test(m)
pval <- ifelse(f$p.value == 0, 2.2e-308, f$p.value)
text = paste("Fisher's exact test\n  OR = ",round(f$estimate,2),"\np-value <", format(pval, scientific = TRUE))

ggplot(d2,aes(x = dataset, y = V1, fill = dataset)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.7, alpha = 0.5, size = 1) + 
  annotate(geom = "text", x = -Inf, y = Inf, label = text, size = 8, hjust = -0.3, vjust = 2) +
  geom_text(label = paste(round(d2$V1,1),"%"), vjust = 1.5, size = 8 ) + 
  scale_fill_brewer(palette = "Set2") +
  xlab("Partial correlation") +
  ylab("% enh-enh significant") +
  guides(fill = FALSE) +
  theme_minimal() + 
  theme(text = element_text(size=30),  aspect.ratio = 1,
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=2.5))

