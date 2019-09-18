rm(list=ls(all=TRUE))

library(seqinr)
library(ggplot2)

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

tCount = sapply(shortSeq, nuclCount, char='T')
cCount = sapply(shortSeq, nuclCount, char='C')

tCountlong = sapply(longSeq, nuclCount, char='T')
cCountlong = sapply(longSeq, nuclCount, char='C')

plot(1:length(tCountlong), tCountlong)
plot(1:length(cCountlong), cCountlong)

shortTC = as.data.frame(cbind(cCount, tCount))
shortTC$cCount = - (shortTC$cCount / 20)
shortTC$tCount = shortTC$tCount / 20
num = 1:nrow(shortTC)
shortTC = cbind(shortTC, num)

longTC = as.data.frame(cbind(cCountlong, tCountlong))
longTC$cCountlong = - (longTC$cCountlong / 20)
longTC$tCountlong = longTC$tCountlong / 20
num = 1:nrow(longTC)
longTC = cbind(longTC, num)


pdf('../../Body/4Figures/WholeGenomeAnalyses.NoOverlapAgAlongGenome.pdf')

a = ggplot(shortTC, aes(num, tCount)) +
  geom_bar(aes(fill = 'red'), stat = "identity") +
  geom_bar(aes(num, cCount, fill = 'green'), stat = "identity") +
  ggtitle('low Generation time') + xlab('Position') + ylab('') +
  scale_fill_discrete(name = "Nucleotide", labels = c("G", "A")); a

a = ggplot(longTC, aes(num, tCountlong)) +
  geom_bar(aes(fill = 'red'), stat = "identity") +
  geom_bar(aes(num, cCountlong, fill = 'green'), stat = "identity") +
  ggtitle('high Generation time') + xlab('Position') + ylab('') +
  scale_fill_discrete(name = "Nucleotide", labels = c("G", "A")); a

# dev.off()

shortTC = as.data.frame(cbind(cCount, tCount))
shortTC$cCount = - shortTC$cCount
num = 1:nrow(shortTC)
shortTC = cbind(shortTC, num)

longTC = as.data.frame(cbind(cCountlong, tCountlong))
longTC$cCountlong = - longTC$cCountlong
num = 1:nrow(longTC)
longTC = cbind(longTC, num)

p <- ggplot(shortTC, aes(x=num, y=tCount)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a blue color
  geom_bar(aes(fill = 'red'), stat="identity") +
  geom_bar(aes(num, cCount, fill='green'), stat = "identity") +
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-20, 20) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-2,4), "cm")     # This remove unnecessary margin around plot
  ) +
  
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = 0) + scale_fill_discrete(name = "Nucleotide", labels = c("G", "A"))
p

p <- ggplot(longTC, aes(x=num, y=tCountlong)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a blue color
  geom_bar(stat="identity", aes(fill = 'red')) +
  geom_bar(aes(num, cCountlong, fill='green'), stat = "identity") +
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-20, 20) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-2,4), "cm")     # This remove unnecessary margin around plot
  ) +
  
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = 0) +
  scale_fill_discrete(name = "Nucleotide", labels = c("G", "A"))
p

############################################################
### add zeros

shortTC = as.data.frame(cbind(cCount, tCount))
shortTC$cCount = - shortTC$cCount
a = as.data.frame(matrix(0, ncol = 2, nrow = 35))
names(a) = names(shortTC)
shortTC = rbind(shortTC, a)
num = 1:nrow(shortTC)
shortTC = cbind(shortTC, num)

longTC = as.data.frame(cbind(cCountlong, tCountlong))
longTC$cCountlong = - longTC$cCountlong
a = as.data.frame(matrix(0, ncol = 2, nrow = 37))
names(a) = names(longTC)
longTC = rbind(longTC, a)
num = 1:nrow(longTC)
longTC = cbind(longTC, num)

p <- ggplot(shortTC, aes(x=num, y=tCount)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a blue color
  geom_bar(stat="identity", aes(fill = 'red')) +
  geom_bar(aes(num, cCount, fill='green'), stat = "identity") +
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-20, 20) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-2,4), "cm")     # This remove unnecessary margin around plot
  ) +
  
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = 0) +
  scale_fill_discrete(name = "Nucleotide", labels = c("G", "A"))
p


p <- ggplot(longTC, aes(x=num, y=tCountlong)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a blue color
  geom_bar(stat="identity", aes(fill = 'red')) +
  geom_bar(aes(num, cCountlong, fill='green'), stat = "identity") +
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-20, 20) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-2,4), "cm")     # This remove unnecessary margin around plot
  ) +
  
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = 0) +
  scale_fill_discrete(name = "Nucleotide", labels = c("G", "A"))
p

dev.off()
