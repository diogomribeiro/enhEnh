#11-Aug-2022 Diogo Ribeiro @ UNIL

library(data.table)
library(ggplot2)

##########Load/process data
data = fread("/home/dribeiro/EnhEnhPaper/revision/pbmc/hic/pbmc_hic_contacts_significant.bed.gz", header = T, sep = "\t")
options("scipen"=100, "digits"=2)

##### reaplce NA by 0
data[is.na(normalised_contact)]$normalised_contact = 0

##### log the HiC_cotact values
data$normalised_contact = log2(data$normalised_contact+1)

# accounting for distance 
data$dist = abs(data$centralStart - data$cisStart)
# correlation test between distance and Hi-C values
cor.test(data$dist,data$normalised_contact)
data$res = residuals(lm(data$normalised_contact ~ data$dist ))

######## Real regions  
p = cor.test(data[real == 1]$corr,data[real == 1]$res, method = "spearman")
text = paste("Spearman R:", round(p$estimate,2), "P-val:", signif(p$p.value, digits = 3))
ggplot(data[real == 1], aes(x = corr, y = res) ) +
  geom_bin_2d(bins = 80) +
  geom_smooth(method = "lm") +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.05, vjust = 1.5, size = 10, fontface = "bold"  ) +
  xlab("Enhancer-enhancer correlation") +
  ylab("Hi-C contact (log distance-scaled)") +
  theme_linedraw() +
  theme(text = element_text(size = 40), # big screen
              plot.title = element_text(hjust = 0.5), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 2), aspect.ratio = 1  )

###### Control regions
data$group = "real"
data[real == 0]$group = "control"
ggplot(data, aes(x = group, y = normalised_contact, fill = group) ) +
  geom_boxplot(size = 1) +
  stat_summary(fun=mean, geom="text", size=6, color="black",
               vjust = 0, aes(label= paste( round(..y.., digits = 2))))+
  scale_fill_manual(values = c("#B3AF8F","#66999B")) +
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

