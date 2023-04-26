# 25-Jan-2023 Diogo Ribeiro @ UNIL

library(data.table)
library(ppcor)

nCells = 24844

atacData = fread("rep3_atac_peak_epimap_enhancer_merge_0.5.bed.gz")
atacData$tag = paste(atacData$V1,atacData$V2,atacData$V3,sep = "_")
# for SHARE-seq
atacData$V7 = gsub(",P1.02","",as.character(atacData$V7))
atacData$V7 = gsub(",P1.03","",as.character(atacData$V7))

rnaData = fread("GSM4156603_GM12878.rep3.rna.counts.ensg.final.binary.coding.nofilt.sparse.bed.gz", header = F, sep = "\t")
# for SHARE-seq
colnames(rnaData) = c("gene","cell","value")
rnaData$cell = gsub(",P1.51","",as.character(rnaData$cell))
rnaData$cell = gsub(",P1.50","",as.character(rnaData$cell))

enhEnhData = fread("~/EnhEnhPaper/data/enh_enh_correlation.tsv")


# for Share-SEQ
wantedCells = fread("rep3_rna_atac_cells.txt",header = F)
wantedCells$tag = paste(wantedCells$V1,wantedCells$V2,wantedCells$V3, sep = ",")
atacData = atacData[V7 %in% wantedCells$tag]
rnaData = rnaData[cell %in% wantedCells$tag]

paste("Processing",nrow(enhEnhData),"entries..")

count = 0
finalDT = data.table()
for (i in seq(nrow(enhEnhData))){
  count = count + 1
  if (count %% 100 == 0){print(paste("Processed",count,"entries"))}
  
  g = enhEnhData[i,]$gene
  e1 = enhEnhData[i,]$enh1
  e2 = enhEnhData[i,]$enh2
  
  gene = rnaData[gene == g]$cell
  enh1 = atacData[tag == e1]$V7
  enh2 = atacData[tag == e2]$V7
  
  common3 = intersect(gene,enh1)
  common3 = intersect(common3,enh2)
  
  common2_1 = intersect(gene,enh1)
  common2_1 = setdiff(common2_1, common3)
  
  common2_2 = intersect(gene,enh2)
  common2_2 = setdiff(common2_2, common3)
  
  common2_3 = intersect(enh1,enh2)
  common2_3 = setdiff(common2_3, common3)
  
  uniq1 = setdiff(gene, common3)
  uniq1 = setdiff(uniq1, common2_1)
  uniq1 = setdiff(uniq1, common2_2)
  
  uniq2 = setdiff(enh1, common3)
  uniq2 = setdiff(uniq2, common2_1)
  uniq2 = setdiff(uniq2, common2_3)
  
  uniq3 = setdiff(enh2, common3)
  uniq3 = setdiff(uniq3, common2_2)
  uniq3 = setdiff(uniq3, common2_3)
  
  rest = nCells - (length(common3) + length(common2_1) + length(common2_2) + length(common2_3) + length(uniq1) + length(uniq2) + length(uniq3))
  
  g = c(rep(1,length(common3)), rep(1, length(common2_1)),  rep(1, length(common2_2)), rep(0, length(common2_3)), rep(1, length(uniq1)),  rep(0, length(uniq2)), rep(0, length(uniq3)), rep(0,rest)  )
  e1 = c(rep(1,length(common3)), rep(1, length(common2_1)),  rep(0, length(common2_2)), rep(1, length(common2_3)), rep(0, length(uniq1)),  rep(1, length(uniq2)), rep(0, length(uniq3)), rep(0,rest)  )
  e2 = c(rep(1,length(common3)), rep(0, length(common2_1)),  rep(1, length(common2_2)), rep(1, length(common2_3)), rep(0, length(uniq1)),  rep(0, length(uniq2)), rep(1, length(uniq3)), rep(0,rest)  )
  
  sum(g) == length(gene)
  sum(e1) == length(enh1)
  sum(e2) == length(enh2)
  
  data <- data.table(gene = g, enh1 = e1, enh2 = e2)
  
  d = pcor(data)
  ge1 = d$estimate[1,2]
  ge2 = d$estimate[1,3]
  e1e2 = d$estimate[2,3]
  
  ge1p = d$p.value[1,2]
  ge2p = d$p.value[1,3]
  e1e2p = d$p.value[2,3]
  
  sge1 = cor.test(data$gene, data$enh1)  
  sge2 = cor.test(data$gene, data$enh2)  
  se1e2 = cor.test(data$enh1, data$enh2)  
  
  sge1e = sge1$estimate
  sge2e = sge2$estimate
  se1e2e = se1e2$estimate
  
  sge1p = sge1$p.value
  sge2p = sge2$p.value
  se1e2p = se1e2$p.value

  finalDT = rbind(finalDT, data.table( gene = enhEnhData[i,]$gene, enhancer1 = enhEnhData[i,]$enh1, enhancer2 = enhEnhData[i,]$enh2,
                                       gene_enh1_pcor = ge1, gene_enh2_pcor = ge2, enh1_enh2_pcor = e1e2,  
                                       gene_enh1_pcor_pval = ge1p, gene_enh2_pcor_pval = ge2p, enh1_enh2_pcor_pval = e1e2p,
                                       gene_enh1_scor = sge1e, gene_enh2_scor = sge2e, enh1_enh2_scor = se1e2e,
                                       gene_enh1_scor_pval = sge1p, gene_enh2_scor_pval = sge2p, enh1_enh2_scor_pval = se1e2p
  ) )
  
  
}

finalDT

write.table(finalDT, "share_seq_partial_corr_all.out",sep="\t", row.names = F, quote=F)
