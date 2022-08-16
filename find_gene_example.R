

geneEnh = fread("/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/enh_enh_paper/data/gene_enhancer_significant.tsv")
geneFreq = data.table(table(geneEnh$gene))
wantedGenes = geneFreq[N == 3]

data = fread("/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/enh_enh_paper/data/significant_enh_enh.tsv")
data = data[gene %in% wantedGenes$V1]

enhEnhFreq = data.table(table(data$gene))

data = data[gene %in% enhEnhFreq[N == 3]$V1]

data[order(corr)]

data[gene == "ENSG00000100439"]
