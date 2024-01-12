#Picture one
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

df_mtdna = read.csv('../../Paper_materials_2024/Birds_dataset_paper.csv')

df_nd6 = read.csv('../../Paper_materials_2024/Birds_mtDNA_data.csv')
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
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

graph6 = ggplot(data = df_nd6, aes(x = gene_name, y = ThChSkew))+
  geom_boxplot(notch = TRUE)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))+
  ylim(-1,1)+
  xlab('Mitochondrial genes')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
skew_graph= ggarrange(graph5, graph6,
                      ncol = 2, nrow = 1)

unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip", exdir = "../../Body/3Results/")
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
new_b_and_m = ggplot(new_big, aes(x = gene_name, y = GhAhSkew, fill = Class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Mitochondrial genes')+
  ylab('GhAhSkew')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5",'CYTB',"ND6","ND1","ND2"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")

new_b_and_m_one = ggplot(new_big, aes(x = gene_name, y = ThChSkew, fill = Class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Mitochondrial genes')+
  ylab('ThChSkew')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5",'CYTB',"ND6","ND1","ND2"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")

for_article = SynNuc[, c(1, 2, 73, 74, 75, 76)]
for_article$class = 'Mammalia'
new_bird1 = df_nd6[, c('species_name', 'gene_name', 'neutral_A','neutral_g', 'neutral_c', 'neutral_T')]
new_bird1$class = 'Aves'
new_bird1$species_name = gsub(' ', '_', new_bird1$species_name)
for_article$Gene[for_article$Gene == 'CytB'] = 'CYTB'
names(for_article) = c('species_name', 'gene_name', 'neutral_A', 'neutral_T', 'neutral_g', 'neutral_c', 'class')

new_big1 = rbind(for_article, new_bird1)
new_b_and_m1 = ggplot(new_big1, aes(x = gene_name, y = neutral_A, fill = class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Gene name')+
  ylab('Neutral T')+
  ylim(0, 200)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5",'CYTB',"ND6","ND1","ND2"))+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        legend.position = "none")

new_b_and_m2 = ggplot(new_big1, aes(x = gene_name, y = neutral_g, fill = class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  ylab('Neutral C')+
  xlab('Mitochondrial genes')+
  ylim(0, 200)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5",'CYTB',"ND6","ND1","ND2"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")

new_b_and_m3 = ggplot(new_big1, aes(x = gene_name, y = neutral_c, fill = class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Gene name')+
  ylab('Neutral G')+
  ylim(0, 200)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5",'CYTB',"ND6","ND1","ND2"))+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        legend.position = "none")


new_b_and_m4 = ggplot(new_big1, aes(x = gene_name, y = neutral_T, fill = class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Mitochondrial genes')+
  ylab('Neutral A')+
  ylim(0, 200)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5",'CYTB',"ND6","ND1","ND2"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")


freq_graph= ggarrange(graph1, graph4, graph2, graph3,
                      ncol = 2, nrow = 2)
freq_graph

skew_graph

mam_vs_aves_fr = ggarrange(new_b_and_m3, new_b_and_m1, new_b_and_m4, new_b_and_m2,
                           ncol = 2, nrow = 2)

mam_vs_aves_fr

mam_vs_aves_skew = ggarrange(new_b_and_m, new_b_and_m_one,
                             ncol = 2, nrow = 1)
mam_vs_aves_skew
