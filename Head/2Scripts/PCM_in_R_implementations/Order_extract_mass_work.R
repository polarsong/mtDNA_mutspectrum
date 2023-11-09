rm(list=ls(all=T))

library(ggtree)
library(ggplot2)
library(treeio)
library(dplyr)
library(ape)
library(geiger)
library(ggpubr)
df_mtdna = read.csv('Birds_dataset_paper.csv')
df_mtdna1 = df_mtdna
df_mtdna1$taxonomy = sub('formes.*', '', df_mtdna1$taxonomy)
df_mtdna1$taxonomy = sub(".* ", "", df_mtdna1$taxonomy)
df_mtdna1$taxonomy = paste(df_mtdna1$taxonomy, 'formes', sep = '')
df_tax = unique(df_mtdna1[,c(2,5,117)])
df_tax$logMass = log10(df_tax$Mass)

df_need = data.frame()
for (i in unique(df_mtdna$species_name))
{
  a = df_mtdna[df_mtdna$species_name == i,]
  b = sum(a$ghahSkew)/12
  ab = c(i,b)
  df_need = rbind(df_need, ab)
}
names(df_need) = c('species_name', 'GhAhSkew')
df_tax = merge(df_tax, df_need)
df_tax$GhAhSkew = as.numeric(df_tax$GhAhSkew)
counter = 0
for (i in unique(df_tax$taxonomy))
{
  a = df_tax[df_tax$taxonomy == i,]
  b = ggplot(a, aes(x = logMass, y = GhAhSkew))+
    geom_point()+
    ggtitle(i)+
    geom_smooth(method = 'lm', formula = y~x)
  assign(paste0("Order_", i), b)
  counter = counter + 1
}

ggarrange(Order_Accipitriformes, Order_Anseriformes, Order_Apodiformes, Order_Apterygiformes, Order_Bucerotiformes, Order_Caprimulgiformes, Order_Casuariiformes, Order_Charadriiformes, Order_Ciconiiformes, Order_Coliiformes, Order_Columbiformes, Order_Coraciiformes, Order_Cuculiformes, Order_Falconiformes, Order_Galbuliformes, Order_Galliformes, Order_Gaviiformes, Order_Gruiformes, Order_Musophagiformes, Order_Passeriformes, Order_Pelecaniformes, Order_Phoenicopteriformes, Order_Piciformes, Order_Podicipediformes, Order_Procellariiformes, Order_Psittaciformes, Order_Rheiformes, Order_Sphenisciformes, Order_Strigiformes, Order_Struthioniformes, Order_Tinamiformes, Order_Trogoniformes, Order_Upupiformes, 
          ncol = 3, nrow = 11)
