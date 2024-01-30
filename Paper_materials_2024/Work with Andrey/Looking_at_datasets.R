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
df11 = read.csv ('Species_life-histories.csv')
df12 = read.csv('GlobalBMRbase.csv')
df_mtdna = read.csv('../../Paper_materials_2024/Birds_dataset_paper.csv')
names_v = unique(df_mtdna$species_name)
df_short = data.frame()
for (i in names_v)
{
  df1 = df_mtdna[df_mtdna$species_name == i,]
  a = sum(df1$ghahSkew)/12
  b = sum(df1$chthSkew)/12
  ab = c(i, a, b)
  df_short = rbind(df_short, ab)
}
names(df_short) = c('Species', 'GhAhSkew', 'ThChSkew')
df_short$GhAhSkew = as.numeric(df_short$GhAhSkew)
df_mtdna11 = merge(df_short, df11, by = 'Species')
ggplot(df_mtdna11, aes(x = Migration, y = GhAhSkew))+
  geom_boxplot()
ggplot(df_mtdna11, aes(x = Habitat, y = GhAhSkew))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(df_mtdna11, aes(x = Foraging_environment, y = GhAhSkew))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(df_mtdna11, aes(x = Daily_activity, y = GhAhSkew))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(df_mtdna11, aes(x = Diet, y = GhAhSkew))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(df_mtdna11, aes(x = Clutch, y = GhAhSkew))+
  geom_point()
ggplot(df_mtdna11, aes(x = NestingPeriod, y = GhAhSkew))+
  geom_point()
