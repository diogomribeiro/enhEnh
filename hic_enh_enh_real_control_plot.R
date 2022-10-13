#11-Aug-2022 Diogo Ribeiro @ UNIL

library(data.table)
library(ggplot2)

##########Load/process data
data = fread("~/EnhEnhPaper/data/Hi-C/hic_contacts_significant.bed", header = T, sep = "\t")
options("scipen"=100, "digits"=2)

##### reaplce NA by 0
data[is.na(normalised_contact)]$normalised_contact = 0

##### log the HiC_cotact values
data$normalised_contact = log2(data$normalised_contact+1)

# accounting for distance 
data$dist = abs(data$centralStart - data$cisStart)
# correlation test between distance and Hi-C values
cor.test(data$dist,data$normalised_contact) ### cor : -0.8086667  / p-value < 2.2e-16
data$res = residuals(lm(data$normalised_contact ~ data$dist ))

######## Real regions  
p = cor.test(data[real == 1]$corr,data[real == 1]$res, method = "spearman")
text = paste("Spearman R:", round(p$estimate,2), "P-val:", signif(p$p.value, digits = 3))
ggplot(data[real == 1], aes(x = corr, y = res) ) +
  # geom_point(alpha = 0.5, size = 0.1) +
  geom_bin_2d(bins = 80) +
  geom_smooth(method = "lm") +
  # scale_fill_gradient(high = "#cb181d", low = "#fcbba1") +
  # annotate("text", x = Inf, y = Inf, label = text, hjust = 1.05, vjust = 1.5, size = 6, fontface = "bold"  ) +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.05, vjust = 1.5, size = 11, fontface = "bold"  ) +
  xlab("Enhancer-enhancer correlation") +
  ylab("Hi-C contact (log distance-scaled)") +
  theme_linedraw() +
  # theme(text = element_text(size = 24), # small screen
  theme(text = element_text(size = 50), # big screen
              plot.title = element_text(hjust = 0.5), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 2), aspect.ratio = 1  )

###### Control regions
data$group = "real"
data[real == 0]$group = "control"
ggplot(data, aes(x = group, y = normalised_contact, fill = group) ) +
  geom_boxplot(size = 1) +
  stat_summary(fun=mean, geom="text", size=6, color="black",
               vjust = 0, aes(label= paste( round(..y.., digits = 2))))+
  # geom_text(aes(label=..count..), y=0, stat='count', colour="black", size=5)+
  scale_fill_manual(values = c("#B3AF8F","#66999B")) +
  # scale_fill_brewer(palette = "Set2") +
  xlab("enhancer-enhancer regions") +
  ylab("Hi-C contact (log-scaled)") +
  theme_minimal() +
  theme(text = element_text(size = 22), plot.title = element_text(hjust = 0.5), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1  )
t = wilcox.test(data[group == "real"]$normalised_contact, data[group == "control"]$normalised_contact)
t
t$p.value

## Direct comparison
real = data[real == 1][,.(centralStart,cisStart,corr,match_id,normalised_contact)]
control = data[real == 0][,.(centralStart,cisStart,corr,match_id,normalised_contact)]
mergedData = merge(real,control,by = "match_id")

mergedData$diff = mergedData$normalised_contact.x - mergedData$normalised_contact.y
options(digits = 3)
text = paste(round(nrow(mergedData[diff > 0])/nrow(mergedData) * 100,2),"% cases real is higher control",sep="")
ggplot(mergedData, aes(x = diff) ) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.05, vjust = 1.5, size = 5, fontface = "bold"  ) +
  xlab("Contact difference") +
  ylab("Count") +
  theme_minimal() +
  theme(text = element_text(size = 18), plot.title = element_text(hjust = 0.5), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1  )


# ### Significant vs non-significant
# # Note: Hi-C processing means enh coordinates are midpoints
# data = data[real == 1]
# data$tag=paste("chr",data$`#chr`,"_",data$centralStart,"|","chr",data$`#chr`,"_",data$cisStart,sep="")
# simpleData = unique(data[,.(normalised_contact,res,tag)])
# 
# match = fread("~/EnhEnhPaper/data/sign_nonsign_enh_enh_match.txt")
# match$sig1 = data.table(unlist(lapply(match$sig, function(x) unlist(strsplit(x,"[|]"))[1])))$V1
# match$sig2 = data.table(unlist(lapply(match$sig, function(x) unlist(strsplit(x,"[|]"))[2])))$V1
# match$nonsig1 = data.table(unlist(lapply(match$nonsig, function(x) unlist(strsplit(x,"[|]"))[1])))$V1
# match$nonsig2 = data.table(unlist(lapply(match$nonsig, function(x) unlist(strsplit(x,"[|]"))[2])))$V1
# 
# match$sig1_chr = data.table(unlist(lapply(match$sig1, function(x) unlist(strsplit(x,"[_]"))[1])))$V1
# match$sig1_st = as.numeric(data.table(unlist(lapply(match$sig1, function(x) unlist(strsplit(x,"[_]"))[2])))$V1)
# match$sig1_end = as.numeric(data.table(unlist(lapply(match$sig1, function(x) unlist(strsplit(x,"[_]"))[3])))$V1)
# match$sig1_mid = (match$sig1_st + match$sig1_end) / 2
# 
# match$sig2_chr = data.table(unlist(lapply(match$sig2, function(x) unlist(strsplit(x,"[_]"))[1])))$V1
# match$sig2_st = as.numeric(data.table(unlist(lapply(match$sig2, function(x) unlist(strsplit(x,"[_]"))[2])))$V1)
# match$sig2_end = as.numeric(data.table(unlist(lapply(match$sig2, function(x) unlist(strsplit(x,"[_]"))[3])))$V1)
# match$sig2_mid = (match$sig2_st + match$sig2_end) / 2
# 
# match$sig1_coord = paste(match$sig1_chr,match$sig1_mid,sep="_")
# match$sig2_coord = paste(match$sig2_chr,match$sig2_mid,sep="_")
# match$sig_tag = paste(match$sig2_coord,match$sig1_coord,sep="|")
# 
# match$nonsig1_chr = data.table(unlist(lapply(match$nonsig1, function(x) unlist(strsplit(x,"[_]"))[1])))$V1
# match$nonsig1_st = as.numeric(data.table(unlist(lapply(match$nonsig1, function(x) unlist(strsplit(x,"[_]"))[2])))$V1)
# match$nonsig1_end = as.numeric(data.table(unlist(lapply(match$nonsig1, function(x) unlist(strsplit(x,"[_]"))[3])))$V1)
# match$nonsig1_mid = (match$nonsig1_st + match$nonsig1_end) / 2
# 
# match$nonsig2_chr = data.table(unlist(lapply(match$nonsig2, function(x) unlist(strsplit(x,"[_]"))[1])))$V1
# match$nonsig2_st = as.numeric(data.table(unlist(lapply(match$nonsig2, function(x) unlist(strsplit(x,"[_]"))[2])))$V1)
# match$nonsig2_end = as.numeric(data.table(unlist(lapply(match$nonsig2, function(x) unlist(strsplit(x,"[_]"))[3])))$V1)
# match$nonsig2_mid = (match$nonsig2_st + match$nonsig2_end) / 2
# 
# match$nonsig1_coord = paste(match$nonsig1_chr,match$nonsig1_mid,sep="_")
# match$nonsig2_coord = paste(match$nonsig2_chr,match$nonsig2_mid,sep="_")
# match$nonsig_tag = paste(match$nonsig2_coord,match$nonsig1_coord,sep="|")
# 
# signMatch = simpleData[tag %in% match$sig_tag]
# signMatch$significant = "yes"
# nonsignMatch = simpleData[tag %in% match$nonsig_tag]
# nonsignMatch$significant = "no"
# matchDT = rbind(signMatch,nonsignMatch)
# 
# ggplot(matchDT, aes(x = significant, y = normalised_contact, fill = significant) ) +
#   geom_boxplot() +
#   stat_summary(fun=mean, geom="text", size=5, color="black",
#                vjust = 0, aes(label= paste( round(..y.., digits = 2))))+
#   # geom_text(aes(label=..count..), y=0, stat='count', colour="black", size=5)+
#   scale_fill_brewer(palette = "Set2") +
#   xlab("enhancer-enhancer regions") +
#   ylab("Hi-C contact (log-scaled)") +
#   theme_minimal() +
#   theme(text = element_text(size = 18), plot.title = element_text(hjust = 0.5), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
#         panel.background = element_rect(colour = "black", fill = "white", size = 1), aspect.ratio = 1  )
# wilcox.test(signMatch$normalised_contact,nonsignMatch$normalised_contact)
# wilcox.test(signMatch$res,nonsignMatch$res)


