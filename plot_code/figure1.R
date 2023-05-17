# 16-May-2023 Diogo Ribeiro

library(data.table)
library(ggplot2)
library(patchwork)

# pdf("~/Figure1_panel_b_c.pdf",11,8)
png("~/Figure1_panel_b_c.png",2000,1000)

###############
# Panel B
###############
data = fread("~/git/enhEnh/source_data/data_figure_1b.tsv.gz",sep="\t",header=T)

## Correlation distribution
g1 = ggplot(data, aes(x=corr, fill = significant) ) + 
  geom_histogram( color = "black", binwidth = 0.01, position = "stack", size = 1.2) +
  scale_fill_manual(values = c("#B3AF8F","#66999B")) +
  geom_vline(xintercept = 0.05, color = "#525252", linetype = "dashed", size = 2.5)+
  labs(tag = "b") +
  xlab("Enhancer-enhancer correlation") +
  ylab("Frequency") +
  theme_linedraw() + 
  theme(text = element_text(size=50), legend.position = c(0.5,0.8),
        legend.text=element_text(size=45),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3), aspect.ratio = 1)  


dt = data.table(cutoff = c("BH 5% &\n Corr > 0.05","BH 5% &\n Corr > 0.05"),  significant = c("yes","no"), values = c(nrow(data[significant == "yes"]), nrow(data[significant == "no"])) )
dt$perc = dt$values / nrow(data) * 100
dt$y = dt$values
dt[significant == "no"]$y = nrow(data)
g2 = ggplot(dt, aes(x= cutoff, y = values, fill = significant ) ) + 
  geom_bar( stat = "identity", position = "stack", color = "black", size = 1) +
  scale_fill_manual(values = c("#B3AF8F","#66999B")) +
  geom_vline(xintercept = 0.05, color = "grey", linetype = "dashed") +
  annotate("text", x = 1, y = 0.5, size = 11, vjust = 4,  label = paste0("N=", dt$values[1],"\n", paste0(round(dt$perc[1],1),"%")), fontface = "bold") + 
  annotate("text", x = 1, y = 0.5, size = 11, vjust = 7.5,  label = paste0("N=", dt$values[2],"\n", paste0(round(dt$perc[2],1),"%")), fontface = "bold") + 
  xlab("") +
  ylab("Enhancer pairs") +
  theme_linedraw() + 
  theme(text = element_text(size=40), legend.position = "none",
        axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0))  


###############
# Panel C
###############

data = fread("~/git/enhEnh/source_data/data_figure_1c.tsv.gz",sep="\t",header=T)

p = cor.test(data$corr,data$res, method = "spearman")
text = paste("Spearman R:", round(p$estimate,2), "P-val =", signif(p$p.value, digits = 3))
g3 = ggplot(data, aes(x = corr, y = res) ) +
  geom_bin_2d(bins = 80) +
  geom_smooth(method = "lm") +
  annotate("text", x = Inf, y = Inf, label = text, hjust = 1.05, vjust = 1.5, size = 13, fontface = "bold"  ) +
  xlab("Enhancer-enhancer correlation") +
  ylab("Hi-C contact (log distance-scaled)") +
  labs(tag = "c") +
  theme_linedraw() +
  theme(text = element_text(size = 50), legend.position = c(0.85,0.20), legend.text=element_text(size=24), legend.title=element_text(size=30),
        plot.title = element_text(hjust = 0.5), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background = element_rect(colour = "black", fill = "white", size = 3), aspect.ratio = 1  )

g1 + inset_element(g2, left = 0.65, bottom = 0.8, right = 0.95, top = 0.1, on_top = TRUE, align_to = "panel") + g3


dev.off()
