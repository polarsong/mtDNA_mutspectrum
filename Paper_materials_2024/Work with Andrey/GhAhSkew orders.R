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
unique(df_mtdna$Habitat)
df_eco = read.csv('final_birds_list_with_no_mistakes.csv')
df_eco = df_eco[,c(2,6)]
names(df_eco) = c('species_name', 'foraging_niche')
df_eco = na.omit(df_eco)
df_mtdna1 = merge(df_mtdna, df_eco, by = 'species_name')
df_mtdna1 = unique(df_mtdna1)
ggplot(df_mtdna1, aes(x = foraging_niche, y = ghahSkew))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
df_mtdna2 = df_mtdna1[df_mtdna1$foraging_niche == 'Vertivore ground',]
names_v = unique(df_mtdna1$species_name)
df_short = data.frame()
for (i in names_v)
{
  df1 = df_mtdna1[df_mtdna1$species_name == i,]
  a = sum(df1$ghahSkew)/12
  b = sum(df1$chthSkew)/12
  v = unique(df1$foraging_niche)
  ab = c(i, a, b, v)
  df_short = rbind(df_short, ab)
}
names(df_short) = c('species_name', 'GhAhSkew', 'ThChSkew', 'Foraging_niche')
df_short$GhAhSkew = as.numeric(df_short$GhAhSkew)
ggplot(df_short, aes(x = Foraging_niche, y = GhAhSkew))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))