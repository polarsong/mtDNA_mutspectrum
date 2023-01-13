rm(list=ls(all=TRUE))

library(seqinr)
library(ggplot2)

############ Syn mut
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")){file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")}

names(SynNuc)

SynNuc = SynNuc[!(SynNuc$Gene %in% c('ND1', 'ND2', 'ND6')),]


VecOfSynFourFoldDegenerateSites <- c('CTT', 'CTC', 'CTA', 'CTG', 
                                     'GTT', 'GTC', 'GTA', 'GTG', 
                                     'TCT', 'TCC', 'TCA', 'TCG', 
                                     'CCT', 'CCC', 'CCA', 'CCG', 
                                     'ACT', 'ACC', 'ACA', 'ACG', 
                                     'GCT', 'GCC', 'GCA', 'GCG', 
                                     'CGT', 'CGC', 'CGA', 'CGG', 
                                     'GGT', 'GGC', 'GGA', 'GGG')

df_mtdna = read.csv('../../Head/2Scripts/Birds_dataset_paper.csv')
df_shift = read.csv('../../Head/2Scripts/Aminoacids_shift_birds.csv')

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

SynNuc$ThirdPos4f = lapply(as.character(SynNuc$CodonsNoOverlap), getVecOfThirdNucleotides4f) #CodonsNoOverlap -> sequence
df_mtdna$ThirdPos4f = lapply(as.character(df_mtdna$sequence), getVecOfThirdNucleotides4f)
df_need = df_mtdna[,c('species_name',"gene_name", "sequence", "ThirdPos4f", "realm")]

#M = merge(SynNuc, GT, by ='Species')

df_need$gene_name =  ordered(df_need$gene_name, levels = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','CYTB', "ND1", "ND2"))
df_need = df_need[order(df_need$gene_name),]

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

tCount = sapply(df_need, nuclCount, char='T')
cCount = sapply(shortSeq, nuclCount, char='C')

tCountlong = sapply(longSeq, nuclCount, char='T')
cCountlong = sapply(longSeq, nuclCount, char='C')

plot(1:length(tCountlong), tCountlong)
plot(1:length(cCountlong), cCountlong)

shortTC = as.data.frame(cbind(cCount, tCount))
shortTC$cCount = - shortTC$cCount
num = 1:nrow(shortTC)
shortTC = cbind(shortTC, num)

longTC = as.data.frame(cbind(cCountlong, tCountlong))
longTC$cCountlong = - longTC$cCountlong
num = 1:nrow(longTC)
longTC = cbind(longTC, num)


pdf('../../Body/4Figures/WholeGenomeAnalyses.NoOverlapAgAlongGenome.pdf')

ColG = rgb(0.1,0.1,0.1,0.5)
ColA = rgb(1,0.1,0.1,0.5)

a = ggplot(shortTC, aes(num, tCount)) +
  geom_bar(fill = ColA, color = 'white', stat = "identity") +
  geom_bar(aes(num, cCount), fill = ColG, color = 'white', stat = "identity") +
  ggtitle('low Generation time') + xlab('Position') + ylab('') +
  scale_fill_discrete(name = "Nucleotide", labels = c("G", "A")) +
  theme_classic(); a

a = ggplot(longTC, aes(num, tCountlong)) +
  geom_bar(fill = ColA, color = 'white', stat = "identity") +
  geom_bar(aes(num, cCountlong), fill = ColG, color = 'white', stat = "identity") +
  ggtitle('high Generation time') + xlab('Position') + ylab('') +
  scale_fill_discrete(name = "Nucleotide", labels = c("G", "A")) +
  theme_classic(); a

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
  geom_bar(fill = ColA, color = 'white', stat="identity") +
  geom_bar(aes(num, cCount), fill = ColG, color = 'white', stat = "identity") +
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
# p

p <- ggplot(longTC, aes(x=num, y=tCountlong)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a blue color
  geom_bar(stat="identity", fill = ColA, color = 'white') +
  geom_bar(aes(num, cCountlong), fill = ColG, color = 'white', stat = "identity") +
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
# p

############################################################
### add zeros

shortTC = as.data.frame(cbind(cCount, tCount))
shortTC$cCount = - shortTC$cCount
a = as.data.frame(matrix(0, ncol = 2, nrow = 40))
names(a) = names(shortTC)
shortTC = rbind(shortTC, a)
num = 1:nrow(shortTC)
shortTC = cbind(shortTC, num)

longTC = as.data.frame(cbind(cCountlong, tCountlong))
longTC$cCountlong = - longTC$cCountlong
a = as.data.frame(matrix(0, ncol = 2, nrow = 35))
names(a) = names(longTC)
longTC = rbind(longTC, a)
num = 1:nrow(longTC)
longTC = cbind(longTC, num)

p <- ggplot(shortTC, aes(x=num, y=tCount)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a blue color
  geom_bar(stat="identity", fill = ColA, color = 'white') +
  geom_bar(aes(num, cCount), fill = ColG, color = 'white', stat = "identity") +
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-20, 20) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    #plot.title = element_text('Tarsipes rostratus, 341.275')
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
    # plot.margin = unit(rep(-1,4), "cm")     # This remove unnecessary margin around plot
  ) +
  
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = 3*pi/4, direction = -1) +
  scale_fill_discrete(name = "Nucleotide", labels = c("G", "A"))
print(p + labs(title = 'Tarsipes rostratus, 341 days'))

p <- ggplot(longTC, aes(x=num, y=tCountlong)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a blue color
  geom_bar(stat="identity", fill = ColA, color = 'white') +
  geom_bar(aes(num, cCountlong), fill = ColG, color = 'white', stat = "identity") +
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-20, 20) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
    # plot.margin = unit(rep(-2,4), "cm")     # This remove unnecessary margin around plot
  ) +
  
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = 2*pi/3, direction = -1) +
  scale_fill_discrete(name = "Nucleotide", labels = c("G", "A"))
p + ggtitle('Balaena mysticetus, 18980 days')

dev.off()


#work from internet
install.packages('reshape2')
library(reshape2)
df_shift = df_shift[,c(2:67)]
df_shift1 <- melt(df_shift ,  id.vars = 'gene_name', variable.name = 'gene')
df_shift2 = df_shift[,c(2:66)]
df_shift3 = melt(df_shift2, id.vars = 'gene_name')


df <- data.frame(index=c(1, 2, 3, 4, 5, 6),
                 var1=c(4, 4, 5, 4, 3, 2),
                 var2=c(1, 2, 4, 4, 6, 9),
                 var3=c(9, 9, 9, 5, 5, 3))

#melt data frame into long format
df1 <- melt(df ,  id.vars = 'index', variable.name = 'series')

ggplot(df_shift3, aes(x = variable, y = value, fill = gene_name))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90)) 

names(df_shift3)

#plotly work
library(plotly)
plot_ly(
  data = df_shift3,
  y = ~value,
  x = ~variable,
  type = "box",
  color = ~gene_name,
  showlegend = FALSE
)
