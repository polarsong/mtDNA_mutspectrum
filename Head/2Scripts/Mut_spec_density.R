rm(list = ls(all=TRUE))
library(ggbiplot)
library(ggplot2)
library(ggpubr)
library(plotly)
df_mtdna = read.csv('../../Head/2Scripts/Birds_dataset_paper.csv')


df_pca = df_mtdna[c('species_name','gene_name','fAn', 'fGn', 'fCn', 'fTn', 'Stg_Sac','ghahSkew', 'chthSkew', 'med_T', 'med_G')]
gene_vector = c('fAn', 'fGn', 'fCn', 'fTn', 'Stg_Sac','ghahSkew', 'chthSkew', 'med_T', 'med_G')
gene_stats = data.frame(unique(df_pca$species_name))
for ( char in gene_vector){
  
  stats1 = aggregate(df_pca[,char], by = list(df_pca$species_name), FUN = 'sum')[2]
  stats1 = stats1/12
  gene_stats = cbind(gene_stats, stats1)
  
}
names(gene_stats) = c('species_name', gene_vector)
df_realm = df_mtdna[c('species_name', 'realm', 'Trophic_niche')]
gene_stats = merge(gene_stats, df_realm, by = 'species_name')
gene_stats = unique(gene_stats)
row.names(gene_stats) = gene_stats$species_name
#gene_stats$species_name = NA
gene_stats = gene_stats[, colSums(is.na(gene_stats)) < nrow(gene_stats)]
stats_pca = prcomp(gene_stats[c(2,3,4,5,6,7,8,9,10)], center = TRUE, scale. = TRUE)
summary(stats_pca)

birds_pca = data.frame(stats_pca$x)
birds_pca = birds_pca[,c(1,2)]
birds_pca$species_name = row.names(birds_pca)
gene_stats$species_name = row.names(gene_stats)
gene_stats = merge(gene_stats, birds_pca, by = 'species_name')
row.names(gene_stats) = gene_stats$species_name
gene_stats = gene_stats[,c(2:14)]
b_bipl = ggbiplot(stats_pca, groups = gene_stats$Trophic_niche, labels = gene_stats$species_name,labels.size = 2)
b_bipl
b_bipl1 = ggbiplot(stats_pca, groups = gene_stats$realm, labels = gene_stats$species_name,labels.size = 2)
b_bipl1


g1 = ggplot(gene_stats, aes(x=PC1, color=realm)) +
  geom_density(size = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (48.0%)')

g1
g2 = ggplot(gene_stats, aes(x=PC2, color=realm)) +
  geom_density(size = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC2 (32.2%)')


g2

g3 = ggplot(gene_stats, aes(x=PC1, color=Trophic_niche)) +
  geom_density(size = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (48.0%)')

g3
g4 = ggplot(gene_stats, aes(x=PC2, color=Trophic_niche)) +
  geom_density(size = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC2 (32.2%)')


g4

g5 = ggplot(gene_stats, aes(x=PC1, color=realm)) +
  geom_density(size = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (48.0%)')+
  scale_colour_manual(name="Origin", values= c("black", "red", "black", "black", "black", "black", "black", "black", "black"))

g5

g6 = ggplot(gene_stats, aes(x=PC2, color=realm)) +
  geom_density(size = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC2 (32.2%)')+
  scale_colour_manual(name="Origin", values= c("black", "red", "black", "black", "black", "black", "black", "black", "black"))


g6

g7 = ggplot(gene_stats, aes(x=PC1, color=Trophic_niche)) +
  geom_density(size = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (48.0%)')+
  scale_colour_manual(name="Origin", values= c("black", "black", "black", "red", "black", "black", "black", "black", "black", 'black'))

g7
g8 = ggplot(gene_stats, aes(x=PC2, color=Trophic_niche)) +
  geom_density(size = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC2 (32.2%)')+
  scale_colour_manual(name="Origin", values= c("black", "black", "black", "red", "black", "black", "black", "black", "black", 'black'))


g8

df_mtdna = read.csv('../../Head/2Scripts/Birds_dataset_paper.csv')
mut_data = read.table("C:/Users/User/Desktop/Birds mutspec results from Bogdan/mutspec12.tsv", header = TRUE, fill = TRUE)
unique(mut_data$AltNode)
mut_data_ff = mut_data[mut_data$Label == 'ff',]
mut_data_ff = mut_data_ff[,c(1,2,3,4,5,7,8)]
ecozone_data = df_mtdna[,c('species_name', 'realm', 'Trophic_niche')]
ecozone_data = unique(ecozone_data)
ecozone_data$species_name = gsub(' ', '_', ecozone_data$species_name)
mut_data_ff = mut_data_ff[!grepl('Node', mut_data_ff$AltNode),]
names(mut_data_ff) = c('Mut', 'ObsNum', 'ExpNum', 'RawMutSpec', 'MutSpec', 'species_name', 'Label') 
mut_data_ff = merge(mut_data_ff, ecozone_data, by = 'species_name')
mut_data_ff$MutSpec =  as.numeric(mut_data_ff$MutSpec)

AC = mut_data_ff[mut_data_ff$Mut == 'A>C',] 
AG = mut_data_ff[mut_data_ff$Mut == 'A>G',]
AT = mut_data_ff[mut_data_ff$Mut == 'A>T',]
GC = mut_data_ff[mut_data_ff$Mut == 'G>C',]
GT = mut_data_ff[mut_data_ff$Mut == 'G>T',]
GA = mut_data_ff[mut_data_ff$Mut == 'G>A',]
CG = mut_data_ff[mut_data_ff$Mut == 'C>G',]
CT = mut_data_ff[mut_data_ff$Mut == 'C>T',]
CA = mut_data_ff[mut_data_ff$Mut == 'C>A',]
TG = mut_data_ff[mut_data_ff$Mut == 'T>G',]
TC = mut_data_ff[mut_data_ff$Mut == 'T>C',]
TA = mut_data_ff[mut_data_ff$Mut == 'T>A',]

AC = replace(AC, 'A>C', 'T>G')
AC = AC[,c(1,3,4,5,6,7,8,9,10)]
names(AC) = c('species_name', 'ObsNum', 'ExpNum', 'RawMutSpec', 'MutSpec', 'Label', 'realm', 'Trophic_niche', 'Mut')

AG = replace(AG, 'A>G', 'T>C')
AG = AG[,c(1,3,4,5,6,7,8,9,10)]
names(AG) = c('species_name', 'ObsNum', 'ExpNum', 'RawMutSpec', 'MutSpec', 'Label', 'realm', 'Trophic_niche', 'Mut')

AT = replace(AT, 'A>T', 'T>A')
AT = AT[,c(1,3,4,5,6,7,8,9,10)]
names(AT) = c('species_name', 'ObsNum', 'ExpNum', 'RawMutSpec', 'MutSpec', 'Label', 'realm', 'Trophic_niche', 'Mut')

GC = replace(GC, 'G>C', 'C>G')
GC = GC[,c(1,3,4,5,6,7,8,9,10)]
names(GC) = c('species_name', 'ObsNum', 'ExpNum', 'RawMutSpec', 'MutSpec', 'Label', 'realm', 'Trophic_niche', 'Mut')

GT = replace(GT, 'G>T', 'C>A')
GT = GT[,c(1,3,4,5,6,7,8,9,10)]
names(GT) = c('species_name', 'ObsNum', 'ExpNum', 'RawMutSpec', 'MutSpec', 'Label', 'realm', 'Trophic_niche', 'Mut')

GA = replace(GA, 'G>A', 'C>T')
GA = GA[,c(1,3,4,5,6,7,8,9,10)]
names(GA) = c('species_name', 'ObsNum', 'ExpNum', 'RawMutSpec', 'MutSpec', 'Label', 'realm', 'Trophic_niche', 'Mut')

CG = replace(CG, 'C>G', 'G>C')
CG = CG[,c(1,3,4,5,6,7,8,9,10)]
names(CG) = c('species_name', 'ObsNum', 'ExpNum', 'RawMutSpec', 'MutSpec', 'Label', 'realm', 'Trophic_niche', 'Mut')

CT = replace(CT, 'C>T', 'G>A')
CT = CT[,c(1,3,4,5,6,7,8,9,10)]
names(CT) = c('species_name', 'ObsNum', 'ExpNum', 'RawMutSpec', 'MutSpec', 'Label', 'realm', 'Trophic_niche', 'Mut')

CA = replace(CA, 'C>A', 'G>T')
CA = CA[,c(1,3,4,5,6,7,8,9,10)]
names(CA) = c('species_name', 'ObsNum', 'ExpNum', 'RawMutSpec', 'MutSpec', 'Label', 'realm', 'Trophic_niche', 'Mut')

TG = replace(TG, 'T>G', 'A>C')
TG = TG[,c(1,3,4,5,6,7,8,9,10)]
names(TG) = c('species_name', 'ObsNum', 'ExpNum', 'RawMutSpec', 'MutSpec', 'Label', 'realm', 'Trophic_niche', 'Mut')

TC = replace(TC, 'T>C', 'A>G')
TC = TC[,c(1,3,4,5,6,7,8,9,10)]
names(TC) = c('species_name', 'ObsNum', 'ExpNum', 'RawMutSpec', 'MutSpec', 'Label', 'realm', 'Trophic_niche', 'Mut')

TA = replace(TA, 'T>A', 'A>T')
TA = TA[,c(1,3,4,5,6,7,8,9,10)]
names(TA) = c('species_name', 'ObsNum', 'ExpNum', 'RawMutSpec', 'MutSpec', 'Label', 'realm', 'Trophic_niche', 'Mut')

mut_data_ff1 = rbind(AC, AG, AT, GC, GT, GA, CT, CA, CG, TA, TG, TC)

library(tidyr)
library(data.table)

pca_data = mut_data_ff1[,c(1,5,9)]
ex = reshape(data = pca_data, idvar = 'species_name',
             timevar = 'Mut',
             direction = 'wide')
names(ex) = c('species_name', 'T>G', 'T>C', 'T>A', 'C>G', 'C>A', 'C>T', 'G>A','G>T','G>C', 'A>T', 'A>C', 'A>G')
ex = merge(ex, ecozone_data, by = 'species_name')
ex = ex[,c(1,14,15,12,13,11,6,5,7,8,10,9,4,3,2)]
row.names(ex) = ex$species_name


stats_pca1 = prcomp(ex[,c(4,5,6,7,8,9,10,11,12,13,14,15)], center = TRUE, scale. = TRUE)
summary(stats_pca1)



bipl = ggbiplot(stats_pca1, groups = ex$realm, labels.size = 2)
bipl

bipl1 = ggbiplot(stats_pca1, groups = pca_data_shaped$Trophic_niche, labels.size = 2)
bipl1


birds_ms_pca = data.frame(stats_pca1$x)
birds_ms_pca = birds_ms_pca[,c(1,2)]
birds_ms_pca$species_name = row.names(birds_ms_pca)
ex2 = merge(ex, birds_ms_pca, by = 'species_name')

g11 = ggplot(ex2, aes(x=PC1, color=realm)) +
  geom_density(size = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (37.9%)')

g11
g12 = ggplot(ex2, aes(x=PC2, color=realm)) +
  geom_density(size = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (10.0%)')

g12

g13 = ggplot(ex2, aes(x=PC1, color=Trophic_niche)) +
  geom_density(size = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (37.9%)')

g13
g14 = ggplot(ex2, aes(x=PC2, color=Trophic_niche)) +
  geom_density(size = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (10.0%)')

g14

g15 = ggplot(ex2, aes(x=PC1, color=realm)) +
  geom_density(size = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (37.9%)')+
  scale_colour_manual(name="Origin", values= c("black", "red", "black", "black", "black", "black", "black", "black", "black"))


g15
g16 = ggplot(ex2, aes(x=PC2, color=realm)) +
  geom_density(size = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (10.0%)')+
  scale_colour_manual(name="Origin", values= c("black", "red", "black", "black", "black", "black", "black", "black", "black"))


g16

g17 = ggplot(ex2, aes(x=PC1, color=Trophic_niche)) +
  geom_density(size = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (37.9%)')+
  scale_colour_manual(name="Origin", values= c("black", "black", "black", "red", "black", "black", "black", "black", "black", 'black'))

g17
g18 = ggplot(ex2, aes(x=PC2, color=Trophic_niche)) +
  geom_density(size = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (10.0%)')+
  scale_colour_manual(name="Origin", values= c("black", "black", "black", "red", "black", "black", "black", "black", "black", 'black'))

g18

