#17-Aug-2022 Diogo Ribeiro @ UNIL
# Script to match enh-enh pairs by distance

library(data.table)

maxDistancePerc = 0.05

minimumCells = 100

enhEnhData = fread("~/EnhEnhPaper/data/enh_enh_correlation.tsv", header = T, sep = "\t")
# enhEnhData = fread("~/EnhEnhPaper/data/abc/enh_enh_correlation.tsv", header = T, sep = "\t")
enhEnhData$totalCells = enhEnhData$oneOne + enhEnhData$oneZero + enhEnhData$zeroOne + enhEnhData$zerozero
enhEnhData = enhEnhData[totalCells >= minimumCells]

signData = fread("~/EnhEnhPaper/data/significant_enh_enh.tsv", header = T, sep = "\t")
# signData = fread("~/EnhEnhPaper/data/abc/significant_enh_enh.tsv", header = T, sep = "\t")

enhPairs = unique(enhEnhData[,.(enh1,enh2)])
sigEnhPairs = unique(signData[,.(enh1,enh2)])

enhPairs$enhancer1_start = as.numeric(data.table(unlist(lapply(enhPairs$enh1, function(x) unlist(strsplit(x,"[_]"))[2])))$V1)
enhPairs$enhancer1_end = as.numeric(data.table(unlist(lapply(enhPairs$enh1, function(x) unlist(strsplit(x,"[_]"))[3])))$V1)
enhPairs$enhancer2_start = as.numeric(data.table(unlist(lapply(enhPairs$enh2, function(x) unlist(strsplit(x,"[_]"))[2])))$V1)
enhPairs$enhancer2_end = as.numeric(data.table(unlist(lapply(enhPairs$enh2, function(x) unlist(strsplit(x,"[_]"))[3])))$V1)
enhPairs$mid1 = (enhPairs$enhancer1_start + enhPairs$enhancer1_end) / 2
enhPairs$mid2 = (enhPairs$enhancer2_start + enhPairs$enhancer2_end) / 2
enhPairs$dist = abs(enhPairs$mid1 - enhPairs$mid2)


enhPairs$tag = paste(enhPairs$enh1, enhPairs$enh2,sep="|")
sigEnhPairs$tag = paste(sigEnhPairs$enh1, sigEnhPairs$enh2,sep="|")

neverSignPairs = enhPairs[!tag %in% sigEnhPairs$tag]
signPairs = enhPairs[tag %in% sigEnhPairs$tag]

paste("Processing",nrow(neverSignPairs),"enh-enh pairs")
finalDT = data.table()
for (i in seq(1,nrow(neverSignPairs))){
  
  if (i %% 1000 == 0 ) print(paste("Processed",i,"pairs"))
  
  wantedDist = neverSignPairs[i,]$dist
  lowerDist = wantedDist - wantedDist * maxDistancePerc
  upperDist = wantedDist + wantedDist * maxDistancePerc
  
  possiblePairs = signPairs[dist > lowerDist & dist < upperDist]
  
  if (nrow(possiblePairs) > 0){
    
  rand = sample(nrow(possiblePairs),1)
  chosen = possiblePairs[rand]
  # remove chosen  
  signPairs = signPairs[tag != chosen$tag]

  finalDT = rbind(finalDT, data.table(nonsig = neverSignPairs[i,]$tag, sig = chosen$tag))
  } else{
    print(paste("No match for",neverSignPairs[i,]$tag))
  }
  
}

finalDT

write.table(finalDT,"~/EnhEnhPaper/data/sign_nonsign_enh_enh_match.txt",quote=F,sep="\t",row.names=F)
# write.table(finalDT,"~/EnhEnhPaper/data/abc/sign_nonsign_enh_enh_match.txt",quote=F,sep="\t",row.names=F)


