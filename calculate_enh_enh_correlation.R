# 09-Aug-2022 Diogo Ribeiro, Chaymae Ziyani @ UNIL
# Script to create the enhancer-enhancer correlation file 
# Requires file with all cells per gene, file with all cells with RNA-seq and ATAC-seq, and a file with gene-enhancer-cell

library(data.table)

args = commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
  stop("All arguments need to be provided", call.=FALSE)
}

# Background of cells with ATAC-seq and RNA-seq
cellBackground = fread(args[1], header = F, sep = "\t") # simple list of cell IDs to use, one cell per line
# Gene-enhancer-cell file
geneEnhCell= fread(args[2], header = T, sep = "\t") # TSV file with header with "gene", "tag" and "cell" columns (tag refers to chr_start_end of enhancers) 
# Gene expression
geneDT = fread(args[3], header = T, sep = "\t") # TSV file with header with "gene" and "cell" columns
# Filter for cells with ATAC-seq data
geneDT = geneDT[cell %in% cellBackground$V1]

paste("Number of genes in gene-enhancer-cell file:",length(unique(geneEnhCell$gene)))
paste("Number of genes in gene expression data:", length(unique(geneDT$gene)))
paste("Number of cells in gene expression data:", length(unique(geneDT$cell)))

###  Enh-enh correlation analysis
df = data.table()
geneEnhEnhCells = data.table()
count=0
### Loop through the geneID: calculate the cor between enh across genes and cells
for (geneID in unique(geneEnhCell$gene)){
  data = geneEnhCell[gene == geneID]
  geneEx = geneDT[geneID == gene]$cell
  listEnh = unique(data$tag)
  done = data.table()
  
  count=count+1
  if (count %% 100 == 0 ) print(paste("Processed",count,"genes"))
  listEnh = sort(listEnh)
  for (enh1 in listEnh){
    for (enh2 in listEnh){
      
      done = rbind(done, data.table(enh1 = enh1, enh2=enh2))
      if (enh1 %in% done$enh2 & enh2 %in% done$enh1) next
      if (enh1 == enh2) next
      
      L1 = data[data$tag == enh1]$cell
      L2 = data[data$tag == enh2]$cell
      
      overlap = intersect(L1 , L2)
      oneOne = length(overlap)
      oneZero = length(L1) - oneOne
      zeroOne = length(L2) - oneOne
      zerozero = length(unique(geneEx)) - (oneOne + oneZero + zeroOne)
      
      C1= c(rep(1,oneOne) , rep(1, oneZero), rep(0, zeroOne) , rep(0, zerozero))
      C2= c(rep(1,oneOne) , rep(0,oneZero), rep(1, zeroOne) , rep(0, zerozero))
      
      if (length(C1) <= 2 | length(C2) <=2 ) next
      
      ### real
      c = cor.test(C1, C2, method = "pearson")
      pval = c$p.value
      corr = c$estimate
      
      geneEnhEnhCells = rbind(geneEnhEnhCells,data.table(geneID, enh1, enh2, cell = overlap))
      
      df = rbind(df,data.table(gene = geneID, enh1=enh1, enh2=enh2, oneOne=oneOne, oneZero=oneZero,
                               zeroOne=zeroOne, zerozero=zerozero, corr = corr,pval = pval))

    }
  }
}

# Write enh-enh correlation
write.table(df, "enh_enh_correlation.tsv",col.names = T , row.names = F , quote = F, sep = "\t" )

# Write gene-enh-enh-cell file
write.table(geneEnhEnhCells, "gene_enh_enh_cell.tsv",col.names = T , row.names = F , quote = F, sep = "\t" )


