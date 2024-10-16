rm(list = ls(all=TRUE))
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(ggbiplot)
library(phytools)
library(nlme)
library(geiger)
library(ggtree)
library(stringr)
library(dplyr)
#preparing datasets
df_mtdna = read.csv('../Work_with_Andrey/Birds_dataset_paper.csv', header = TRUE, sep = ';')
df_nd6 = read.csv('../Birds_mtDNA_data.csv')
df_nd6$GhAhSkew = (df_nd6$neutral_g - df_nd6$neutral_A)/(df_nd6$neutral_g + df_nd6$neutral_A)
df_nd6$ThChSkew = (df_nd6$neutral_T - df_nd6$neutral_c)/(df_nd6$neutral_T + df_nd6$neutral_c)
df_nd6$fTn = df_nd6$neutral_T/df_nd6$neutral_amount
df_nd6$fAn = df_nd6$neutral_A/df_nd6$neutral_amount
df_nd6$fCn = df_nd6$neutral_c/df_nd6$neutral_amount
df_nd6$fGn = df_nd6$neutral_g/df_nd6$neutral_amount

df_nd6$GhAhSkew = (df_nd6$neutral_c- df_nd6$neutral_T)/(df_nd6$neutral_c + df_nd6$neutral_T)
df_nd6$ThChSkew = (df_nd6$neutral_A - df_nd6$neutral_g)/(df_nd6$neutral_A + df_nd6$neutral_g)
df_nd6$fTn = df_nd6$neutral_A/df_nd6$neutral_amount
df_nd6$fAn = df_nd6$neutral_T/df_nd6$neutral_amount
df_nd6$fCn = df_nd6$neutral_g/df_nd6$neutral_amount
df_nd6$fGn = df_nd6$neutral_c/df_nd6$neutral_amount

#picture 1 variant 1

graph1 = ggplot(data = df_nd6, aes(x = gene_name, y = fTn))+
  geom_boxplot(notch = TRUE)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))+
  ylim(0, 0.8)+
  xlab('Mitochondrial genes')+
  ylab('Thymine frequency')+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 


graph2 = ggplot(data = df_nd6, aes(x = gene_name, y = fCn))+
  geom_boxplot(notch = TRUE)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))+
  ylim(0, 0.8)+
  xlab('Mitochondrial genes')+
  ylab('Cytosine frequency')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

graph3 = ggplot(data = df_nd6, aes(x = gene_name, y = fAn))+
  geom_boxplot(notch = TRUE)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))+
  ylim(0, 0.8)+
  xlab('Mitochondrial genes')+
  ylab('Adenine frequency')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

graph4 = ggplot(data = df_nd6, aes(x = gene_name, y = fGn))+
  geom_boxplot(notch = TRUE)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))+
  ylim(0, 0.8)+
  xlab('Mitochondrial genes')+
  ylab('Guanine frequency')+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

graph5 = ggplot(data = df_nd6, aes(x = gene_name, y = GhAhSkew))+
  geom_boxplot(notch = TRUE)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))+
  ylim(-1,1)+
  xlab('Mitochondrial genes')+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

graph6 = ggplot(data = df_nd6, aes(x = gene_name, y = ThChSkew))+
  geom_boxplot(notch = TRUE)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))+
  ylim(-1,1)+
  xlab('Mitochondrial genes')+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

unzip("../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip", exdir = "../../Body/3Results/")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")) file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")
SynNuc$ghahSkew = ((SynNuc$NeutralC - SynNuc$NeutralT))/((SynNuc$NeutralC + SynNuc$NeutralT))
SynNuc$chthSkew = ((SynNuc$NeutralA - SynNuc$NeutralG))/((SynNuc$NeutralA + SynNuc$NeutralG))
new_mam = SynNuc[, c(1, 2, 79, 80)]
new_mam$Сlass = 'Mammalia'
new_bird = df_nd6[, c('species_name', 'gene_name', 'GhAhSkew','ThChSkew')]
new_bird$Сlass = 'Aves'
new_bird$species_name = gsub(' ', '_', new_bird$species_name)
new_mam$Gene[new_mam$Gene == 'CytB'] = 'CYTB'
names(new_mam) = c('species_name', 'gene_name', 'GhAhSkew', 'ThChSkew', 'Class')
names(new_bird) = c('species_name', 'gene_name', 'GhAhSkew', 'ThChSkew', 'Class')

new_big = rbind(new_mam, new_bird)
graph7 = ggplot(new_big, aes(x = gene_name, y = GhAhSkew, fill = Class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Mitochondrial genes')+
  ylab('GhAhSkew')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5",'CYTB',"ND6","ND1","ND2"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")

graph8 = ggplot(new_big, aes(x = gene_name, y = ThChSkew, fill = Class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Mitochondrial genes')+
  ylab('ThChSkew')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5",'CYTB',"ND6","ND1","ND2"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")

ggarrange(graph1, graph4,graph5, graph6, graph2, graph3, graph7, graph8,
          ncol = 4, nrow = 2)

#picture one variant 2

graph1 = ggplot(data = df_nd6, aes(x = gene_name, y = fTn))+
  geom_boxplot(notch = TRUE)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))+
  ylim(0, 0.8)+
  xlab('Mitochondrial genes')+
  ylab('Thymine frequency')+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) 


graph2 = ggplot(data = df_nd6, aes(x = gene_name, y = fCn))+
  geom_boxplot(notch = TRUE)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))+
  ylim(0, 0.8)+
  xlab('Mitochondrial genes')+
  ylab('Cytosine frequency')+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

graph3 = ggplot(data = df_nd6, aes(x = gene_name, y = fAn))+
  geom_boxplot(notch = TRUE)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))+
  ylim(0, 0.8)+
  xlab('Mitochondrial genes')+
  ylab('Adenine frequency')+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

graph4 = ggplot(data = df_nd6, aes(x = gene_name, y = fGn))+
  geom_boxplot(notch = TRUE)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))+
  ylim(0, 0.8)+
  xlab('Mitochondrial genes')+
  ylab('Guanine frequency')+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

graph5 = ggplot(data = df_nd6, aes(x = gene_name, y = GhAhSkew))+
  geom_boxplot(notch = TRUE)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))+
  ylim(-1,1)+
  xlab('Mitochondrial genes')+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

graph6 = ggplot(data = df_nd6, aes(x = gene_name, y = ThChSkew))+
  geom_boxplot(notch = TRUE)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))+
  ylim(-1,1)+
  xlab('Mitochondrial genes')+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

graph7 = ggplot(new_big, aes(x = gene_name, y = GhAhSkew, fill = Class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Mitochondrial genes')+
  ylab('GhAhSkew')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5",'CYTB',"ND6","ND1","ND2"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")

graph8 = ggplot(new_big, aes(x = gene_name, y = ThChSkew, fill = Class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Mitochondrial genes')+
  ylab('ThChSkew')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5",'CYTB',"ND6","ND1","ND2"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")

ggarrange(graph1, graph4, graph2, graph3, graph5, graph6, graph7, graph8,
          ncol = 2, nrow = 4)

#TBSS statistics
col1 = df_nd6[df_nd6$gene_name == 'COX1',]
col2 = df_nd6[df_nd6$gene_name == 'CYTB',]
wilcox.test(col1$GhAhSkew,col2$GhAhSkew)
wilcox.test(col1$ThChSkew, col2$ThChSkew)

col1 = df_nd6[df_nd6$gene_name == 'ND1',]
col2 = df_nd6[df_nd6$gene_name == 'ND2',]
wilcox.test(col1$GhAhSkew,col2$GhAhSkew)
wilcox.test(col1$ThChSkew, col2$ThChSkew)

df_mtdna$ghahSkew = gsub(',', '.', df_mtdna$ghahSkew)
df_mtdna$chthSkew = gsub(',', '.', df_mtdna$chthSkew)
df_mtdna$ghahSkew = as.numeric(as.character(df_mtdna$ghahSkew))
df_mtdna$chthSkew = as.numeric(as.character(df_mtdna$chthSkew))
df_mtdna_cut = df_mtdna[df_mtdna$gene_name != 'ND1',]
df_mtdna_cut = df_mtdna_cut[df_mtdna_cut$gene_name != 'ND2',]
df_mtdna_cut[df_mtdna_cut$Species == "Drepanis coccinea",]$Species = "Vestiaria coccinea"  
b_names = unique(df_mtdna_cut$Species)
spearman_rhos_ghahskew = data.frame()
spearman_rhos_thchskew = data.frame()
tbss = c(1,2,3,4,5,6,7,8,9,10)
for (i in b_names)
{
  df_bird = df_mtdna_cut[df_mtdna_cut$Species == i,]
  speart = cor.test(df_bird$ghahSkew, tbss)
  spearman_rhos_ghahskew = rbind(spearman_rhos, c(i, speart$p.value))
}

for (i in b_names)
{
  df_bird = df_mtdna_cut[df_mtdna_cut$Species == i,]
  speart = cor.test(df_bird$chthSkew, tbss)
  spearman_rhos_thchskew = rbind(spearman_rhos, c(i, speart$p.value))
}