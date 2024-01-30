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
#foraging niche
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
#Migration
df_int = read.csv('../../Body/1Raw/Avonet_data.csv')
df_migr = df_int[,c('Species3', 'Migration')]
names(df_migr) = c('species_name', 'migration')
df_migr_mtdna = merge(df_mtdna, df_migr, by = 'species_name')
df_migr_mtdna$migration = as.character(df_migr_mtdna$migration)
names_v = unique(df_mtdna$species_name)
df_short1 = data.frame()
for (i in names_v)
{
  df1 = df_mtdna[df_mtdna$species_name == i,]
  a = sum(df1$ghahSkew)/12
  b = sum(df1$chthSkew)/12
  c = sum(df1$Mass)/12
  ab = c(i, a, b,c)
  df_short1 = rbind(df_short1, ab)
}
names(df_short1) = c('species_name', 'GhAhSkew', 'ThChSkew', 'Mass')
df_short1$GhAhSkew = as.numeric(df_short1$GhAhSkew)

df_short1 = merge(df_short1, df_migr)
df_short1$migration = as.character(df_short1$migration)
ggplot(df_migr_mtdna, aes(x = migration, y = ghahSkew))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(df_short1, aes(x = migration, y = GhAhSkew))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
wilcox.test(df_short1[df_short1$migration == '2',]$GhAhSkew, df_short1[df_short1$migration == '1',]$GhAhSkew)
wilcox.test(df_short1[df_short1$migration == '2',]$GhAhSkew, df_short1[df_short1$migration == '3',]$GhAhSkew)
wilcox.test(df_short1[df_short1$migration == '3',]$GhAhSkew, df_short1[df_short1$migration == '1',]$GhAhSkew)
#Orders
df_mtdna_orders = unique(df_mtdna[,c(2,5)])
df_mtdna_orders$taxonomy = sub('formes.*', '', df_mtdna_orders$taxonomy)
df_mtdna_orders$taxonomy = sub('.*\\ ', '', df_mtdna_orders$taxonomy)
df_mtdna_orders$taxonomy = paste(df_mtdna_orders$taxonomy, 'formes', sep = '')
#HWI and Kipp's distance
df_metrics = unique(df_mtdna[,c('species_name', 'realm', 'Trophic_niche', 'Hand_wing_index', "Kipps_distance")])
df_short2 = merge(df_short1, df_metrics)
df_short2 = merge(df_short2, df_mtdna_orders)
ggplot(df_short2, aes(x = Hand_wing_index, y = GhAhSkew))+
  geom_point()+
  geom_smooth(method=lm)
ggplot(df_short2, aes(x = log10(Kipps_distance), y = GhAhSkew))+
  geom_point()+
  geom_smooth(method=lm)
ggplot(df_short2, aes(x = Hand_wing_index, y = GhAhSkew))+
  geom_point(aes(colour = factor(taxonomy)))+
  geom_smooth(method=lm)
ggplot(df_short2, aes(x = taxonomy, y = GhAhSkew))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

counter = 0
for (i in unique(df_short2$taxonomy))
{
  a = df_short2[df_short2$taxonomy == i,]
  linearMod <- lm(Hand_wing_index~GhAhSkew, data=a)
  modelSummary = summary(linearMod)
  c = modelSummary$coefficients[,4][2]
  b = ggplot(a, aes(x = Hand_wing_index, y = GhAhSkew))+
    geom_point()+
    ggtitle(i)+
    geom_smooth(method = 'lm', formula = y~x)+                            
    annotate('text', x = (min(a$Hand_wing_index)+1), y = min(a$GhAhSkew), label = round(c, digits = 4))
  assign(paste0("Order_", i), b)
  counter = counter + 1
}
pdf('Birds_HWI.pdf',         
    width = 30,height = 100)
ggarrange(Order_Accipitriformes, Order_Anseriformes, Order_Apodiformes, Order_Apterygiformes, Order_Bucerotiformes, Order_Caprimulgiformes, Order_Casuariiformes, Order_Charadriiformes, Order_Ciconiiformes, Order_Coliiformes, Order_Columbiformes, Order_Coraciiformes, Order_Cuculiformes, Order_Falconiformes, Order_Galbuliformes, Order_Galliformes, Order_Gaviiformes, Order_Gruiformes, Order_Musophagiformes, Order_Passeriformes, Order_Pelecaniformes, Order_Phoenicopteriformes, Order_Piciformes, Order_Podicipediformes, Order_Procellariiformes, Order_Psittaciformes, Order_Rheiformes, Order_Sphenisciformes, Order_Strigiformes, Order_Struthioniformes, Order_Tinamiformes, Order_Trogoniformes, Order_Upupiformes, 
          ncol = 3, nrow = 11)
dev.off()
 
df_short2$Mass = as.numeric(df_short2$Mass)
counter = 0
for (i in unique(df_short2$taxonomy))
{
  a = df_short2[df_short2$taxonomy == i,]
  linearMod <- lm(Mass~GhAhSkew, data=a)
  modelSummary = summary(linearMod)
  c = modelSummary$coefficients[,4][2]
  b = ggplot(a, aes(x = Mass, y = GhAhSkew))+
    geom_point()+
    ggtitle(i)+
    geom_smooth(method = 'lm', formula = y~x)+                            
    annotate('text', x = (min(a$Mass)+1), y = min(a$GhAhSkew), label = round(c, digits = 4))
  assign(paste0("Order_", i), b)
  counter = counter + 1
}
pdf('Birds_mass.pdf',         
    width = 30,height = 100)
ggarrange(Order_Accipitriformes, Order_Anseriformes, Order_Apodiformes, Order_Apterygiformes, Order_Bucerotiformes, Order_Caprimulgiformes, Order_Casuariiformes, Order_Charadriiformes, Order_Ciconiiformes, Order_Coliiformes, Order_Columbiformes, Order_Coraciiformes, Order_Cuculiformes, Order_Falconiformes, Order_Galbuliformes, Order_Galliformes, Order_Gaviiformes, Order_Gruiformes, Order_Musophagiformes, Order_Passeriformes, Order_Pelecaniformes, Order_Phoenicopteriformes, Order_Piciformes, Order_Podicipediformes, Order_Procellariiformes, Order_Psittaciformes, Order_Rheiformes, Order_Sphenisciformes, Order_Strigiformes, Order_Struthioniformes, Order_Tinamiformes, Order_Trogoniformes, Order_Upupiformes, 
          ncol = 3, nrow = 11)
dev.off()

counter = 0
for (i in unique(df_short2$taxonomy))
{
  a = df_short2[df_short2$taxonomy == i,]
  linearMod <- lm(Kipps_distance~GhAhSkew, data=a)
  modelSummary = summary(linearMod)
  c = modelSummary$coefficients[,4][2]
  b = ggplot(a, aes(x = Kipps_distance, y = GhAhSkew))+
    geom_point()+
    ggtitle(i)+
    geom_smooth(method = 'lm', formula = y~x)+                            
    annotate('text', x = (min(a$Kipps_distance)+1), y = min(a$GhAhSkew), label = round(c, digits = 4))
  assign(paste0("Order_", i), b)
  counter = counter + 1
}
pdf('Birds_Kipps_distance.pdf',         
    width = 30,height = 100)
ggarrange(Order_Accipitriformes, Order_Anseriformes, Order_Apodiformes, Order_Apterygiformes, Order_Bucerotiformes, Order_Caprimulgiformes, Order_Casuariiformes, Order_Charadriiformes, Order_Ciconiiformes, Order_Coliiformes, Order_Columbiformes, Order_Coraciiformes, Order_Cuculiformes, Order_Falconiformes, Order_Galbuliformes, Order_Galliformes, Order_Gaviiformes, Order_Gruiformes, Order_Musophagiformes, Order_Passeriformes, Order_Pelecaniformes, Order_Phoenicopteriformes, Order_Piciformes, Order_Podicipediformes, Order_Procellariiformes, Order_Psittaciformes, Order_Rheiformes, Order_Sphenisciformes, Order_Strigiformes, Order_Struthioniformes, Order_Tinamiformes, Order_Trogoniformes, Order_Upupiformes, 
          ncol = 3, nrow = 11)
dev.off()
