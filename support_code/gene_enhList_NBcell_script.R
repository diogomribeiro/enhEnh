########### Script : Create a list of active enhancers in gene and
####   count in how many cells we have these combinations of gene + list of active enhancers

library(data.table)
library(plyr)
library(dplyr)
library(splitstackshape)
library(tidyverse)


######l Read the gene-enh-cell file (with filter)
final_data= fread("/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/enh_enh_paper/data/gene_enhancer_cell_corr_0.05_FDR_5.out.gz", header = T, sep = "\t") 
#### 6,487,801

#### Create a list of active enhancer per gene and cell
final_data1 = final_data %>% group_by(gene, cell) %>%mutate(enhancerList = toString(tag)) %>%as.data.frame()

### remove the enhancer-ID
final_data1$tag=NULL

### data.table
final_data1 = data.table(unique(final_data1)) ### 4,379,755

##### Count in how many cells(same cells) we have the same combination of gene-list of active enhancers
s = final_data1 %>% group_by(gene , enhancerList) %>% summarise("NBcell" = paste(length(cell) ))
s$NBcell = as.numeric(s$NBcell)
s = data.table(s)  #### 458,841 associations gene- list enh - NBcell

###save the file : gene / enhancerList / NBcell
write.table(s, "/work/FAC/FBM/DBC/odelanea/glcoex/dribeiro/single_cell/enh_enh_paper/data/gene_enhancer_combinations.out",col.names = T , row.names = F , quote = F, sep = "\t" )
