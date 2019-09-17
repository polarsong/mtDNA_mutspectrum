rm(list=ls(all=TRUE))

library(seqinr)

############ Syn mut
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")){file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")}

names(SynNuc)

SynNuc = SynNuc[!(SynNuc$Gene %in% c('ND1', 'ND2', 'ND6')),]

GT = read.table("../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt", header = TRUE, sep = '\t')
GT$Species = gsub(' ','_',GT$Scientific_name)
GT = GT[,c(11,13)]
shortestGT = GT[GT$GenerationLength_d == min(GT$GenerationLength_d), 'Species']
longestGT = GT[GT$GenerationLength_d == max(GT$GenerationLength_d), 'Species']


VecOfSynFourFoldDegenerateSites <- c('CTT', 'CTC', 'CTA', 'CTG', 
                                     'GTT', 'GTC', 'GTA', 'GTG', 
                                     'TCT', 'TCC', 'TCA', 'TCG', 
                                     'CCT', 'CCC', 'CCA', 'CCG', 
                                     'ACT', 'ACC', 'ACA', 'ACG', 
                                     'GCT', 'GCC', 'GCA', 'GCG', 
                                     'CGT', 'CGC', 'CGA', 'CGG', 
                                     'GGT', 'GGC', 'GGA', 'GGG')


getVecOfThirdNucleotides4f = function(x){
  codonsVec = splitseq(s2c(x))
  codons4f = codonsVec[codonsVec %in% VecOfSynFourFoldDegenerateSites]
  thirdPosList = gsub("CTA|GTA|TCA|CCA|ACA|GCA|CGA|GGA", 'A', codons4f)
  thirdPosList = gsub("CTT|GTT|TCT|CCT|ACT|GCT|CGT|GGT", 'T', thirdPosList)
  thirdPosList = gsub("CTG|GTG|TCG|CCG|ACG|GCG|CGG|GGG", 'G', thirdPosList)
  thirdPosList = gsub("CTC|GTC|TCC|CCC|ACC|GCC|CGC|GGC", 'C', thirdPosList)
  thirdPosVec = paste(as.vector(thirdPosList), collapse = '')
  return(thirdPosVec)
}

SynNuc$ThirdPos4f = lapply(as.character(SynNuc$CodonsNoOverlap), getVecOfThirdNucleotides4f)

M = merge(SynNuc, GT, by ='Species')

M$Gene =  ordered(M$Gene, levels = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','CytB'))
M = M[order(M$Gene),]

shortestGT = M[M$GenerationLength_d == min(M$GenerationLength_d),]
longestGT = M[M$GenerationLength_d == max(M$GenerationLength_d),]

shortGlVec4f = paste(shortestGT$ThirdPos4f, collapse = '')
longGlVec4f = paste(longestGT$ThirdPos4f, collapse = '')

shortSeq = splitseq(s2c(shortGlVec4f), word = 20)
longSeq = splitseq(s2c(longGlVec4f), word = 20)

nuclCount = function(seq, char){
  cnt = nchar(as.character(seq)) - nchar(gsub(char, "", seq))
  return(cnt)
}

aCount = sapply(shortSeq, nuclCount, char='A')
gCount = sapply(shortSeq, nuclCount, char='G')

pdf('../../Body/4Figures/WholeGenomeAnalyses.NoOverlapAgAlongGenome.pdf')
plot(1:length(aCount), aCount)
plot(1:length(gCount), gCount)
dev.off()
