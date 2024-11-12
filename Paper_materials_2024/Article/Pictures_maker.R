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
  annotate('text', x = 11, y = 0.75, label = 'N = 766')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


graph2 = ggplot(data = df_nd6, aes(x = gene_name, y = fCn))+
  geom_boxplot(notch = TRUE)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))+
  ylim(0, 0.8)+
  xlab('Mitochondrial genes')+
  ylab('Cytosine frequency')+
  annotate('text', x = 11, y = 0.75, label = 'N = 766')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

graph3 = ggplot(data = df_nd6, aes(x = gene_name, y = fAn))+
  geom_boxplot(notch = TRUE)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))+
  ylim(0, 0.8)+
  xlab('Mitochondrial genes')+
  ylab('Adenine frequency')+
  annotate('text', x = 11, y = 0.75, label = 'N = 766')+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

graph4 = ggplot(data = df_nd6, aes(x = gene_name, y = fGn))+
  geom_boxplot(notch = TRUE)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))+
  ylim(0, 0.8)+
  xlab('Mitochondrial genes')+
  ylab('Guanine frequency')+
  annotate('text', x = 11, y = 0.75, label = 'N = 766')+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

graph1_1 = ggarrange(graph3, graph4, graph2, graph1,
                     ncol = 2, nrow = 2)

graph1_1
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

#unzip("../Body/3Results/AllGenesCodonUsageNoOverlap.zip", exdir = "../../Body/3Results/")
#SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
#if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")) file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")
SynNuc = read.table('AllGenesCodonUsageNoOverlap.txt', header = TRUE, sep = '\t')
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
  ylim(-1,1)+
  annotate('text', x = 4.5, y = -0.75, label = 'N birds = 766')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")

graph8 = ggplot(new_big, aes(x = gene_name, y = ThChSkew, fill = Class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Mitochondrial genes')+
  ylab('ThChSkew')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5",'CYTB',"ND6","ND1","ND2"))+
  ylim(-1,1)+
  annotate('text', x = 4.5, y = -0.75, label = 'N mammals = 4356')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")

graph1_2 = ggarrange(graph5, graph6, graph7, graph8,
          ncol = 2, nrow = 2)
graph1_2
#picture one variant 2


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
df_mtdna[df_mtdna$Species == "Drepanis coccinea",]$Species = "Vestiaria coccinea"
df_mtdna[df_mtdna$Species == "Vestiaria coccinea",]$Species = "Drepanis coccinea"
df_mtdna_cut = df_mtdna[df_mtdna$gene_name != 'ND1',]
df_mtdna_cut = df_mtdna_cut[df_mtdna_cut$gene_name != 'ND2',]
b_names = unique(df_mtdna_cut$Species)
spearman_rhos_ghahskew = data.frame()
spearman_rhos_thchskew = data.frame()
tbss = c(1,2,3,4,5,6,7,8,9,10)
for (i in b_names)
{
  df_bird = df_mtdna_cut[df_mtdna_cut$Species == i,]
  speart = cor.test(df_bird$ghahSkew, tbss)
  spearman_rhos_ghahskew = rbind(spearman_rhos_ghahskew, c(i, speart$p.value))
}
tbss_sampl = sample(tbss, 10, replace = TRUE)
spearman_rhos_ghahskew_sample = data.frame()
for (i in b_names)
{
  df_bird = df_mtdna_cut[df_mtdna_cut$Species == i,]
  speart = cor.test(df_bird$ghahSkew, tbss_sampl)
  spearman_rhos_ghahskew_sample = rbind(spearman_rhos_ghahskew_sample, c(i, speart$p.value))
}

names(spearman_rhos_ghahskew) = c('species_name', 'rho_value')
names(spearman_rhos_ghahskew_sample) = c('species_name', 'rho_value')
spearman_rhos_ghahskew$rho_value = as.numeric(as.character(spearman_rhos_ghahskew$rho_value))
spearman_rhos_ghahskew_sample$rho_value = as.numeric(as.character(spearman_rhos_ghahskew_sample$rho_value))
spearman_rhos_ghahskew$rho_log = log10(spearman_rhos_ghahskew$rho_value)


rhogh1 = ggplot(spearman_rhos_ghahskew, aes(x = rho_value))+
  geom_histogram()+
  xlab("Rho value TBSS for GhAhSkew")
rhogh2 =ggplot(spearman_rhos_ghahskew_sample, aes(x = rho_value))+
  geom_histogram()+
  theme(axis.title.y=element_blank())+
  xlab("Rho value TBSS sample for GhAhSkew")
sup1 = ggarrange(rhogh1, rhogh2,
                 nrow = 1, ncol = 2)
sup1  

for (i in b_names)
{
  df_bird = df_mtdna_cut[df_mtdna_cut$Species == i,]
  speart = cor.test(df_bird$chthSkew, tbss)
  spearman_rhos_thchskew = rbind(spearman_rhos_thchskew, c(i, speart$p.value))
}
tbss_sampl = sample(tbss, 10, replace = TRUE)
spearman_rhos_thchskew_sample = data.frame()
for (i in b_names)
{
  df_bird = df_mtdna_cut[df_mtdna_cut$Species == i,]
  speart = cor.test(df_bird$chthSkew, tbss_sampl)
  spearman_rhos_thchskew_sample = rbind(spearman_rhos_thchskew_sample, c(i, speart$p.value))
}

names(spearman_rhos_thchskew) = c('species_name', 'rho_value')
names(spearman_rhos_thchskew_sample) = c('species_name', 'rho_value')
spearman_rhos_thchskew$rho_value = as.numeric(as.character(spearman_rhos_thchskew$rho_value))
spearman_rhos_thchskew_sample$rho_value = as.numeric(as.character(spearman_rhos_thchskew_sample$rho_value))
spearman_rhos_thchskew$rho_log = log10(spearman_rhos_thchskew$rho_value)


rhogh3 = ggplot(spearman_rhos_thchskew, aes(x = rho_value))+
  geom_histogram()+
  xlab("Rho value TBSS for ThChSkew")
rhogh4 = ggplot(spearman_rhos_thchskew_sample, aes(x = rho_value))+
  geom_histogram()+
  theme(axis.title.y=element_blank())+
  xlab("Rho value TBSS sample for ThChSkew")
sup2 = ggarrange(rhogh3, rhogh4,
                 nrow = 1, ncol = 2)
sup2

#Supp materials

#Mass
df_mtdna$Mass = gsub(',', '.', df_mtdna$Mass)
df_mtdna$Mass = as.numeric(as.character(df_mtdna$Mass))
names_v = unique(df_mtdna$Species)
df_short = data.frame()
for (i in names_v)
{
  df1 = df_mtdna[df_mtdna$Species == i,]
  a = sum(df1$ghahSkew)/12
  b = sum(df1$chthSkew)/12
  v = sum(df1$Mass)/12
  ab = c(i, a, b, v)
  df_short = rbind(df_short, ab)
}
names(df_short) = c('Species', 'GhAhSkew', 'ThChSkew', 'Mass')
df_short$Mass = as.numeric(df_short$Mass)
df_short$GhAhSkew = as.numeric(df_short$GhAhSkew)
df_short$ThChSkew = as.numeric(df_short$ThChSkew)
df_short$log_mass = log10(df_short$Mass)
mass_ghskew = ggplot(df_short, aes(x = log_mass, y = GhAhSkew))+
  geom_point()+
  annotate('text', x = 4, y = 0.1, label = 'N = 766')+
  xlab('Decimal logarithm of mass')+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

mass_thskew = ggplot(df_short, aes(x = log_mass, y = ThChSkew))+
  geom_point()+
  xlab('Decimal logarithm of mass')

#Clutch
df_par = read.csv('../Work_with_Andrey/Species_life-histories.csv')
df_mtdna_par = merge(df_short, df_par, by = 'Species')
df_clutch = df_mtdna_par[,c(1,2,3,15)]
df_clutch = na.omit(df_clutch)

clutch_ghskew = ggplot(df_mtdna_par, aes(x = Clutch, y = GhAhSkew))+
  geom_point()+
  annotate('text', x = 4, y = 0.1, label = 'N = 203')+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank())

clutch_thskew = ggplot(df_mtdna_par, aes(x = Clutch, y = ThChSkew))+
  geom_point()+
  theme(axis.title.y=element_blank())


#BMR
df_bmr = read.csv('../Work_with_Andrey/GlobalBMRbase.csv', sep = ';')
df_bmr_e = df_bmr[df_bmr$Trait == 'BMR',]
names_v = unique(df_bmr_e$Species)
df_short_1 = data.frame()
df_bmr_e$TraitValue = gsub(',', '.', df_bmr_e$TraitValue)
df_bmr_e$TraitValue = suppressWarnings(as.numeric(df_bmr_e$TraitValue))
for (i in names_v)
{
  df1 = df_bmr_e[df_bmr_e$Species == i,]
  a1 = sum(df1$TraitValue)
  a2 = nrow(df1)
  a = a1/a2
  b = 'BMR'
  ab = c(i, b, a)
  df_short_1 = rbind(df_short_1, ab)
}
names(df_short_1) = c('Species', 'Trait', 'BMR_value')
df_mtdna_bmr = merge(df_short, df_short_1)
df_mtdna_bmr$BMR_value = as.numeric(df_mtdna_bmr$BMR_value)
bmr_ghskew = ggplot(df_mtdna_bmr, aes(x = BMR_value, y = GhAhSkew))+
  geom_point()+
  annotate('text', x = 2500, y = 0.1, label = 'N = 186')+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank())

bmr_thskew = ggplot(df_mtdna_bmr, aes(x = BMR_value, y = ThChSkew))+
  geom_point()+
  xlab('BMR value')+
  theme(axis.title.y=element_blank())

#Longevity

df_long = read.csv('../Work_with_Andrey/AVES_longevity.csv')
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
df_long$scinam = firstup(df_long$scinam)
names(df_long) = c('Species', 'Longevity', 'Origin', 'Data')
df_long_correct = data.frame()
long_birds = unique(df_long$Species)
for (i in long_birds)
{
  bird = df_long[df_long$Species == i,]
  a = sum(bird$Longevity)/nrow(bird)
  df_long_correct = rbind(df_long_correct, c(i,a))
}
names(df_long_correct) = c('Species', 'Longevity')
df_long_mtdna = merge(df_long_correct, df_short)
df_long_mtdna$Longevity = as.numeric(as.character(df_long_mtdna$Longevity))
long_ghskew = ggplot(df_long_mtdna, aes(x = Longevity, y = GhAhSkew))+
  geom_point()+
  annotate('text', x = 10, y = 0.1, label = 'N = 264')+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank())
long_thskew = ggplot(df_long_mtdna, aes(x = Longevity, y = ThChSkew))+
  geom_point()+
  theme(axis.title.y=element_blank())

#var 1
sup3 = ggarrange(mass_ghskew, clutch_ghskew, bmr_ghskew, long_ghskew, mass_thskew, clutch_thskew, bmr_thskew, long_thskew,
          ncol = 4, nrow = 2)


#picture 2
#ecozone
skew_eco_ghahskew = ggplot(data = df_mtdna, aes(x = realm, y = ghahSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Birds realms')+
  ylab('GhAhSkew')+
  xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  annotate('text', x = 2, y = -0.25, label = 'N = 766')

skew_eco_thchskew = ggplot(data = df_mtdna, aes(x = realm, y = chthSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Birds realms')+
  ylab('ThChSkew')+
  xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#trophic niche

skew_niche_ghahskew = ggplot(data = df_mtdna, aes(x = Trophic_niche, y = ghahSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Trophic niche')+
  ylab('GhAhSkew')+
  xlim(c('Herbivore aquatic', 'Scavenger', 'Vertivore', 'Granivore', 'Herbivore terrestrial', 'Invertivore', 'Aquatic predator', 'Nectarivore', 'Omnivore', 'Frugivore'))+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank())+
  annotate('text', x = 2, y = -0.25, label = 'N = 766')

skew_niche_thchskew = ggplot(data = df_mtdna, aes(x = Trophic_niche, y = chthSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Trophic niche')+
  ylab('ThChSkew')+
  xlim(c('Herbivore aquatic', 'Scavenger', 'Vertivore', 'Granivore', 'Herbivore terrestrial', 'Invertivore', 'Aquatic predator', 'Nectarivore', 'Omnivore', 'Frugivore'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.title.y=element_blank())

ggarrange(skew_eco_ghahskew, skew_niche_ghahskew, skew_eco_thchskew, skew_niche_thchskew,
          ncol = 2, nrow = 2)


#migration
df_int = read.csv('../../Body/1Raw/Avonet_data.csv')
df_migr = df_int[,c('Species3', 'Migration')]
names(df_migr) = c('Species', 'migration')
df_migr_mtdna = merge(df_short, df_migr, by = 'Species')
df_migr_mtdna$migration = as.character(df_migr_mtdna$migration)
df_migr_mtdna[df_migr_mtdna$migration == "1",]$migration = "Resident"
df_migr_mtdna[df_migr_mtdna$migration == "2",]$migration = "Short-distance migration"
df_migr_mtdna[df_migr_mtdna$migration == "3",]$migration = "Long-distance migration"

skew_migr_ghahskew = ggplot(df_migr_mtdna, aes(x = migration, y = GhAhSkew))+
  geom_boxplot()+
  xlim('Resident','Short-distance migration', 'Long-distance migration')+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  annotate('text', x = 3, y = 0, label = 'N = 763')
skew_migr_thchskew = ggplot(df_migr_mtdna, aes(x = migration, y = ThChSkew))+
  geom_boxplot()+
  xlim('Resident','Short-distance migration', 'Long-distance migration')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('Migration')

#daily_activity
skew_act_ghahskew = ggplot(df_mtdna_par, aes(x = Daily_activity, y = GhAhSkew))+
  geom_boxplot()+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank())+
  annotate('text', x = 1.5, y = 0, label = 'N = 210')
skew_act_thchskew = ggplot(df_mtdna_par, aes(x = Daily_activity, y = ThChSkew))+
  geom_boxplot()+
  theme(axis.title.y=element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('Daily activity')

#TNZ
df_temp = read.csv('../Work_with_Andrey/temp_new_data.csv')
df_mtdna_temp = merge(df_short, df_temp, by = 'Species')

skew_tnz_ghahskew = ggplot(df_mtdna_temp, aes(x = TNZ, y = GhAhSkew))+
  geom_point()+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank())+
  annotate('text', x = 25, y = 0.1, label = 'N = 32')
skew_tnz_thchskew = ggplot(df_mtdna_temp, aes(x = TNZ, y = ThChSkew))+
  geom_point()+
  theme(axis.title.y=element_blank())

ggarrange(skew_migr_ghahskew, skew_act_ghahskew, skew_tnz_ghahskew, 
          skew_migr_thchskew, skew_act_thchskew, skew_tnz_thchskew,
          ncol = 3, nrow = 2)

#fly and dive
df_fly = read.csv('../flying_birds.csv')
df_fly = df_fly[,c(2,3,4)]
names(df_fly) = c('species_name', 'flightless', 'diving')
df_fly_clean1 = df_fly[df_fly$flightless =='Flightless',]
df_fly_clean= df_fly[df_fly$flightless == 'Almost_flightless',]
df_fly_clean = na.omit(df_fly_clean)
df_fly_clean1 = na.omit(df_fly_clean1)
df_dive = df_fly
df_fly = df_fly[df_fly$flightless != 'Flightless',]
df_fly = df_fly[df_fly$flightless != 'Almost_flightless',]
df_fly_clean$flightless = 'Tinamiformes'
df_fly_clean1$flightless = 'Tinamiformes'
df_fly_big = rbind(df_fly, df_fly_clean, df_fly_clean1)
names(df_fly_big) = c("Species", 'flightless', 'diving')
df_fly_final = merge(df_fly_big, df_short)
df_fly_final = df_fly_final[df_fly_final$flightless != 'Galliformes',]
df_fly_final[df_fly_final$flightless == '0',]$flightless = 'Flying birds'
df_fly_final$flightless1 = factor(df_fly_final$flightless, levels = c('Flying birds', 'Tinamiformes', 'Apterygiformes', 'Casuariiformes', 'Struthioniformes', 'Rheiformes', "Psittaciformes", "Columbiformes", "Eurypygiformes", "Gruiformes", "Sphenisciformes"))
fly_skew = ggplot(df_fly_final, aes(x = flightless, y = GhAhSkew, color = flightless1))+
  geom_point(position = position_jitter(width = 0.2))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('Birds groups')+
  xlim('Flying birds', 'Tinamiformes', 'Apterygiformes', 'Casuariiformes', 'Struthioniformes', 'Rheiformes', "Psittaciformes", "Columbiformes", "Eurypygiformes", "Gruiformes", "Sphenisciformes")




#trying regression
all_data = merge(df_long_mtdna, df_short_1,  by = 'Species')
all_data = merge(all_data, df_short, by = 'Species')
all_data = merge(all_data, df_mtdna_par, by = 'Species')
all_data$BMR_value = as.numeric(as.character(all_data$BMR_value))
all_data_with_temp = merge(all_data,  df_mtdna_temp, by = 'Species')
all_data = all_data[,c('Longevity', 'BMR_value', 'Mass', 'Migration', 'Clutch', 'GhAhSkew')] #add 1 if needed
all_data$Migration_value = factor(all_data$Migration,
                       levels = c('resident', 'short-distance migrant', 'long-distance migrant'),
                       labels = c(1, 2, 3))
birds_model = lm(GhAhSkew ~ Longevity + BMR_value + Mass + Clutch + Migration_value, data = all_data)
summary(birds_model)
birds_model = lm(GhAhSkew ~ BMR_value + Mass + Clutch + Migration_value, data = all_data)
summary(birds_model)
birds_model = lm(GhAhSkew ~ BMR_value + Mass + Migration_value, data = all_data)
summary(birds_model)
birds_model = lm(GhAhSkew ~ BMR_value + Mass, data = all_data)
summary(birds_model)
birds_model = lm(GhAhSkew ~ Mass, data = all_data)
summary(birds_model)
plot(all_data, col="navy", main="Matrix Scatterplot")


#PGLS for sup and picture 2

#Sup

library(ape); library(phytools);  library(geiger)
pgls_res_table = data.frame()
feathertree <- read.nexus("../Work_with_Andrey/Ultrametric_feathertree.nex")
feathertree$node.label <- NULL # Remove internal node labels (if any)
is.ultrametric(feathertree)
is.binary(feathertree)
is.rooted(feathertree)
df_short$Species = gsub(' ', '_', df_short$Species)
listSkew = df_short$Species
listTree <- feathertree$tip.label
SpeciesToDrop <- setdiff(listTree, listSkew)
#drop.tip(feathertree, SpeciesToDrop) -> Fly_skew_tree

#Mass
rownames(df_short) <- df_short[,1] 
name.check(feathertree, df_short)
df_short[df_short$Species == "Agapornis_pullarius",] = NA
df_short = na.omit(df_short)
df_short[df_short$Species == "Mergus_squamatus",] = NA
df_short = na.omit(df_short)
name.check(feathertree, df_short)

spp = rownames(df_short)
corLambda<-corPagel(value=1,phy=feathertree,form=~spp)
pgls_mass = gls(GhAhSkew~Mass,
                  data=df_short,correlation=corLambda)
a = as.data.frame(summary(pgls_mass)$tTable)
a$lambda_value = summary(pgls_mass)$modelStruct
pgls_res_table = rbind(pgls_res_table, a)

#Clutch
df_clutch_pgls = df_mtdna_par[,c(1,2,3,14)]
df_clutch_pgls$Species = gsub(' ', '_', df_clutch_pgls$Species)
df_clutch_pgls = na.omit(df_clutch_pgls)
listSkew = df_clutch_pgls$Species
listTree <- feathertree$tip.label
SpeciesToDrop <- setdiff(listTree, listSkew)
drop.tip(feathertree, SpeciesToDrop) -> Clutch_skew_tree
rownames(df_clutch_pgls) = df_clutch_pgls[,1]
name.check(Clutch_skew_tree, df_clutch_pgls)
spp = rownames(df_clutch_pgls)
corLambda<-corPagel(value=1,phy=Clutch_skew_tree,form=~spp)
pgls_clutch = gls(GhAhSkew~Clutch,
                data=df_clutch_pgls,correlation=corLambda)
a = as.data.frame(summary(pgls_clutch)$tTable)
a$lambda_value = summary(pgls_clutch)$modelStruct
pgls_res_table = rbind(pgls_res_table, a)

#BMR
df_mtdna_bmr$Species = gsub(' ', '_', df_mtdna_bmr$Species)
listSkew = df_mtdna_bmr$Species
listTree <- feathertree$tip.label
SpeciesToDrop <- setdiff(listTree, listSkew)
drop.tip(feathertree, SpeciesToDrop) -> Bmr_skew_tree
rownames(df_mtdna_bmr) = df_mtdna_bmr[,1]
name.check(Bmr_skew_tree, df_mtdna_bmr)

spp = rownames(df_mtdna_bmr)
corLambda<-corPagel(value=1,phy=Bmr_skew_tree,form=~spp)
pgls_bmr = gls(GhAhSkew~BMR_value,
               data=df_mtdna_bmr,correlation=corLambda)
a = as.data.frame(summary(pgls_bmr)$tTable)
a$lambda_value = summary(pgls_bmr)$modelStruct
pgls_res_table = rbind(pgls_res_table, a)

#longevity
df_long_pgls = df_long_mtdna[,c(1,2,3,4)]
df_long_pgls$Species = gsub(' ', '_', df_long_pgls$Species)
listSkew = df_long_pgls$Species
listTree <- feathertree$tip.label
SpeciesToDrop <- setdiff(listTree, listSkew)
drop.tip(feathertree, SpeciesToDrop) -> Long_skew_tree
rownames(df_long_pgls) = df_long_pgls[,1]
name.check(Long_skew_tree, df_long_pgls)
df_long_pgls[df_long_pgls$Species == "Agapornis_pullarius",] = NA
df_long_pgls = na.omit(df_long_pgls)
name.check(Long_skew_tree, df_long_pgls)

spp = rownames(df_long_pgls)
corLambda<-corPagel(value=1,phy=Long_skew_tree,form=~spp)
pgls_long = gls(GhAhSkew~Longevity,
                data=df_long_pgls,correlation=corLambda)
a = as.data.frame(summary(pgls_long)$tTable)
a$lambda_value = summary(pgls_long)$modelStruct
pgls_res_table = rbind(pgls_res_table, a)

pgls_sup_res_table = pgls_res_table[-c(1,3,5,7),]
pgls_sup_res_table$lambda_value = as.numeric(as.character(pgls_sup_res_table$lambda_value))
write.csv(pgls_sup_res_table, 'Sup_pgls_results.csv')

#Picture 2

#Ecozone+niche
pgls_res_table = data.frame()
df_econiche = df_mtdna[,c('Species', 'realm', 'Trophic_niche')]
df_econiche = unique(df_econiche)
df_econiche$ant_1_other_0 = 0
df_econiche[df_econiche$realm == 'Antarctic',]$ant_1_other_0 = 1
df_econiche$ha_1_other_0 = 0
df_econiche[df_econiche$Trophic_niche == 'Herbivore aquatic',]$ha_1_other_0 = 1
df_econiche$Species = gsub(' ', '_', df_econiche$Species)
df_econiche_analyse = merge(df_econiche, df_short, by = 'Species')
rownames(df_econiche_analyse) <- df_econiche_analyse[,1] 
name.check(feathertree, df_econiche_analyse)
spp = rownames(df_econiche_analyse)
corLambda<-corPagel(value=1,phy=feathertree,form=~spp)
pgls_ecozone = gls(GhAhSkew~ant_1_other_0,
                  data=df_econiche_analyse,correlation=corLambda)
a = as.data.frame(summary(pgls_ecozone)$tTable)
a$lambda_value = summary(pgls_ecozone)$modelStruct
pgls_res_table = rbind(pgls_res_table, a)

pgls_niche = gls(GhAhSkew~ha_1_other_0,
                   data=df_econiche_analyse,correlation=corLambda)
a = as.data.frame(summary(pgls_niche)$tTable)
a$lambda_value = summary(pgls_niche)$modelStruct
pgls_res_table = rbind(pgls_res_table, a)

#Migration
df_migr_pgls = df_migr_mtdna[,c(1,2,6)]
df_migr_pgls$Species = gsub(' ', '_', df_migr_pgls$Species)
df_migr_pgls$res_1_oth_0 = 0
df_migr_pgls[df_migr_pgls$migration == 'Resident',]$res_1_oth_0 = 1 
df_migr_pgls$sd_1_oth_0 = 0
df_migr_pgls[df_migr_pgls$migration == 'Short-distance migration',]$sd_1_oth_0 = 1 
df_migr_pgls$ld_1_oth_0 = 0
df_migr_pgls[df_migr_pgls$migration == 'Long-distance migration',]$ld_1_oth_0 = 1 

rownames(df_migr_pgls) <- df_migr_pgls[,1] 
listSkew = df_migr_pgls$Species
listTree <- feathertree$tip.label
name.check(feathertree, df_migr_pgls)
SpeciesToDrop <- setdiff(listTree, listSkew)
drop.tip(feathertree, SpeciesToDrop) -> Migr_skew_tree
name.check(Migr_skew_tree, df_migr_pgls)
spp = rownames(df_migr_pgls)
corLambda<-corPagel(value=1,phy=Migr_skew_tree, form=~spp)
pgls_migration_res = gls(GhAhSkew~res_1_oth_0,
                   data=df_migr_pgls, correlation=corLambda)
a = as.data.frame(summary(pgls_migration_res)$tTable)
a$lambda_value = summary(pgls_migration_res)$modelStruct
pgls_res_table = rbind(pgls_res_table, a)

spp = rownames(df_migr_pgls)
corLambda<-corPagel(value=1,phy=Migr_skew_tree, form=~spp)
pgls_migration_sd = gls(GhAhSkew~sd_1_oth_0,
                         data=df_migr_pgls, correlation=corLambda)
a = as.data.frame(summary(pgls_migration_sd)$tTable)
a$lambda_value = summary(pgls_migration_sd)$modelStruct
pgls_res_table = rbind(pgls_res_table, a)

spp = rownames(df_migr_pgls)
corLambda<-corPagel(value=1,phy=Migr_skew_tree, form=~spp)
pgls_migration_ld = gls(GhAhSkew~ld_1_oth_0,
                        data=df_migr_pgls, correlation=corLambda)
a = as.data.frame(summary(pgls_migration_ld)$tTable)
a$lambda_value = summary(pgls_migration_ld)$modelStruct
pgls_res_table = rbind(pgls_res_table, a)

df_temp_pgls = df_mtdna_temp[,c(1,2,11)]
df_temp_pgls$Species = gsub(' ', '_', df_temp_pgls$Species)
rownames(df_temp_pgls) <- df_temp_pgls[,1] 
listSkew = df_temp_pgls$Species
listTree <- feathertree$tip.label
name.check(feathertree, df_temp_pgls)
SpeciesToDrop <- setdiff(listTree, listSkew)
drop.tip(feathertree, SpeciesToDrop) -> Temp_skew_tree
name.check(Temp_skew_tree, df_temp_pgls)
spp = rownames(df_temp_pgls)
corLambda<-corPagel(value=1,phy=Temp_skew_tree, form=~spp)
pgls_temp = gls(GhAhSkew~TNZ,
                        data=df_temp_pgls, correlation=corLambda)
a = as.data.frame(summary(pgls_temp)$tTable)
a$lambda_value = summary(pgls_temp)$modelStruct
pgls_res_table = rbind(pgls_res_table, a)

#daily act
df_da_pgls = df_mtdna_par[,c(1,2,3,13)]
df_da_pgls$Species = gsub(' ', '_', df_da_pgls$Species)
df_da_pgls = na.omit(df_da_pgls)
listSkew = df_da_pgls$Species
listTree <- feathertree$tip.label
SpeciesToDrop <- setdiff(listTree, listSkew)
drop.tip(feathertree, SpeciesToDrop) -> Da_skew_tree
rownames(df_da_pgls) = df_da_pgls[,1]
name.check(Da_skew_tree, df_da_pgls)
df_da_pgls$cath_1_other_0 = 0
df_da_pgls[df_da_pgls$Daily_activity == 'cathemeral',]$cath_1_other_0 = 1
df_da_pgls$crep_1_other_0 = 0
df_da_pgls[df_da_pgls$Daily_activity == 'crepuscular',]$crep_1_other_0 = 1
df_da_pgls$diu_1_other_0 = 0
df_da_pgls[df_da_pgls$Daily_activity == 'diurnal',]$diu_1_other_0 = 1
df_da_pgls$noc_1_other_0 = 0
df_da_pgls[df_da_pgls$Daily_activity == 'nocturnal',]$noc_1_other_0 = 1

spp = rownames(df_da_pgls)
corLambda<-corPagel(value=1,phy=Da_skew_tree,form=~spp)
pgls_da_cath = gls(GhAhSkew~cath_1_other_0,
                  data=df_da_pgls,correlation=corLambda)
a = as.data.frame(summary(pgls_da_cath)$tTable)
a$lambda_value = summary(pgls_da_cath)$modelStruct
pgls_res_table = rbind(pgls_res_table, a)
spp = rownames(df_da_pgls)
corLambda<-corPagel(value=1,phy=Da_skew_tree,form=~spp)
pgls_da_crep = gls(GhAhSkew~crep_1_other_0,
                   data=df_da_pgls,correlation=corLambda)
a = as.data.frame(summary(pgls_da_crep)$tTable)
a$lambda_value = summary(pgls_da_crep)$modelStruct
pgls_res_table = rbind(pgls_res_table, a)
spp = rownames(df_da_pgls)
corLambda<-corPagel(value=1,phy=Da_skew_tree,form=~spp)
pgls_da_diu = gls(GhAhSkew~diu_1_other_0,
                   data=df_da_pgls,correlation=corLambda)
a = as.data.frame(summary(pgls_da_diu)$tTable)
a$lambda_value = summary(pgls_da_diu)$modelStruct
pgls_res_table = rbind(pgls_res_table, a)
spp = rownames(df_da_pgls)
corLambda<-corPagel(value=1,phy=Da_skew_tree,form=~spp)
pgls_da_noc = gls(GhAhSkew~noc_1_other_0,
                   data=df_da_pgls,correlation=corLambda)
a = as.data.frame(summary(pgls_da_noc)$tTable)
a$lambda_value = summary(pgls_da_noc)$modelStruct
pgls_res_table = rbind(pgls_res_table, a)

pgls_pict2_res_table = pgls_res_table[-c(1,3,5,7,9,11,13,15,17,19),]
pgls_pict2_res_table$lambda_value = as.numeric(as.character(pgls_pict2_res_table$lambda_value))
write.csv(pgls_pict2_res_table, 'Pict2_pgls_results.csv')


#MutSpec easy
df_mut = read.csv('MutSpecVertebrates12.csv')
df_ac = df_mut[df_mut$Mut == 'A>C',]
df_ac[df_ac$Mut == 'A>C',]$Mut = 'T>G'
df_ag = df_mut[df_mut$Mut == 'A>G',]
df_ag[df_ag$Mut == 'A>G',]$Mut = 'T>C'
df_at = df_mut[df_mut$Mut == 'A>T',]
df_at[df_at$Mut == 'A>T',]$Mut = 'T>A'
df_ca = df_mut[df_mut$Mut == 'C>A',]
df_ca[df_ca$Mut == 'C>A',]$Mut = 'G>T'
df_cg = df_mut[df_mut$Mut == 'C>G',]
df_cg[df_cg$Mut == 'C>G',]$Mut = 'G>C'
df_ct = df_mut[df_mut$Mut == 'C>T',]
df_ct[df_ct$Mut == 'C>T',]$Mut = 'G>A'
df_ga = df_mut[df_mut$Mut == 'G>A',]
df_ga[df_ga$Mut == 'G>A',]$Mut = 'C>T'
df_gc = df_mut[df_mut$Mut == 'G>C',]
df_gc[df_gc$Mut == 'G>C',]$Mut = 'C>G'
df_gt = df_mut[df_mut$Mut == 'G>T',]
df_gt[df_gt$Mut == 'G>T',]$Mut = 'C>A'
df_ta = df_mut[df_mut$Mut == 'T>A',]
df_ta[df_ta$Mut == 'T>A',]$Mut = 'A>T'
df_tc = df_mut[df_mut$Mut == 'T>C',]
df_tc[df_tc$Mut == 'T>C',]$Mut = 'A>G'
df_tg = df_mut[df_mut$Mut == 'T>G',]
df_tg[df_tg$Mut == 'T>G',]$Mut = 'A>C'

df_mut_cor = rbind(df_ac, df_ag, df_at, df_ca, df_cg, df_ct, df_ga, df_gc, df_gt, df_ta, df_tc, df_ag)

df_mut_aves = df_mut_cor[df_mut_cor$Class == 'Aves',]
df_cytb = df_mut_aves[df_mut_aves$Gene == 'Cytb',]

ggplot(df_cytb, aes(x = Mut, y = MutSpec))+
  geom_boxplot()+
  ylab('Mutspec for CytB')

df_cox1 = df_mut_aves[df_mut_aves$Gene == 'CO1',] 
ggplot(df_cox1, aes(x = Mut, y = MutSpec))+
  geom_boxplot()+
  ylab('Mutspec for COX1')
