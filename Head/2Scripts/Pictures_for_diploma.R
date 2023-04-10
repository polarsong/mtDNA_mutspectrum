#Картинки для дипома/практики
rm(list = ls(all=TRUE))
library(ggbiplot)
library(ggplot2)
library(ggpubr)
#Частоты
df_mtdna = read.csv('../../Head/2Scripts/Birds_dataset_paper.csv')
f1 = ggplot(data = df_mtdna, aes(x = gene_name, y = fTn))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Митохондриальные гены')+
  ylab('Частоты тимина')+
  ylim(0, 0.7)
f1 = f1 +xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND1","ND2"))
f1 = f1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
f2 = ggplot(data = df_mtdna, aes(x = gene_name, y = fCn))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Митохондриальные гены')+
  ylab('Частоты цитозина')+
  ylim(0, 0.7)
f2 = f2 +xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND1","ND2"))
f2 = f2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
f3 = ggplot(data = df_mtdna, aes(x = gene_name, y = fGn))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Митохондриальные гены')+
  ylab('Частоты гуанина')+
  ylim(0, 0.7)
f3 = f3 +xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND1","ND2"))
f3 = f3 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
f4 = ggplot(data = df_mtdna, aes(x = gene_name, y = fAn))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Митохондриальные гены')+
  ylab('Частоты аденина')+
  ylim(0, 0.7)
f4 = f4 +xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND1","ND2"))
f4 = f4 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
f5 = ggarrange(f1, f3, f2, f4, 
               ncol = 2, nrow = 2)
f5


#Скосы
skew_all = ggplot(data = df_mtdna, aes(x = gene_name, y = ghahSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Митохондриальные гены')+
  ylab('GhAhSkew')
skew_all = skew_all + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND1","ND2"))
skew_all

skew_all1 = ggplot(data = df_mtdna, aes(x = gene_name, y = chthSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Митохондриальные гены')+
  ylab('ThChSkew')
skew_all1 = skew_all1 + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND1","ND2"))
skew_all1

stg_all = ggplot(data = df_mtdna, aes(x = gene_name, y = Stg_Sac))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Митохондриальные гены')+
  ylab('Stg-Sac')
stg_all = stg_all + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND1","ND2"))
stg_all

#Экология
skew_eco = ggplot(data = df_mtdna, aes(x = realm, y = ghahSkew))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Экозоны')+
  ylab('GhAhSkew')
skew_eco = skew_eco + xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))
skew_eco = skew_eco + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
skew_eco


skew_eco12 = ggplot(data = df_mtdna, aes(x = Trophic_niche, y = ghahSkew))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Трофические ниши')+
  ylab('GhAhSkew')
skew_eco12 = skew_eco12 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
skew_eco12 = skew_eco12 + xlim(c('Herbivore aquatic', 'Scavenger', 'Vertivore', 'Granivore', 'Herbivore terrestrial', 'Invertivore', 'Aquatic predator', 'Nectarivore', 'Omnivore', 'Frugivore'))
skew_eco12


skew_eco2 = ggplot(data = df_mtdna, aes(x = realm, y = chthSkew))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Экозоны')+
  ylab('ThChSkew')
skew_eco2 = skew_eco2 + xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))
skew_eco2 = skew_eco2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
skew_eco2

skew_eco31 = ggplot(data = df_mtdna, aes(x = Trophic_niche, y = chthSkew))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Трофические ниши')+
  ylab('ThChSkew')
skew_eco31 = skew_eco31 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
skew_eco31 = skew_eco31 + xlim(c('Herbivore aquatic', 'Scavenger', 'Vertivore', 'Granivore', 'Herbivore terrestrial', 'Invertivore', 'Aquatic predator', 'Nectarivore', 'Omnivore', 'Frugivore'))
skew_eco31

stg_eco = ggplot(data = df_mtdna, aes(x = realm, y = Stg_Sac))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Экозоны')+
  ylab('Stg-Sac')
stg_eco = stg_eco + xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))
stg_eco = stg_eco + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
stg_eco


stg_eco2 = ggplot(data = df_mtdna, aes(x = Trophic_niche, y = Stg_Sac))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Трофические ниши')+
  ylab('Stg-Sac')
stg_eco2 = stg_eco2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
stg_eco2 = stg_eco2 + xlim(c('Herbivore aquatic', 'Scavenger', 'Vertivore', 'Granivore', 'Herbivore terrestrial', 'Invertivore', 'Aquatic predator', 'Nectarivore', 'Omnivore', 'Frugivore'))
stg_eco2


medT_realm = ggplot(data = df_mtdna, aes(x = realm, y = med_T))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Экозоны')+
  ylab('Асимметрия тимина')
medT_realm = medT_realm + xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))
medT_realm = medT_realm + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
medT_realm

medG_realm = ggplot(data = df_mtdna, aes(x = realm, y = med_G))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Экозоны')+
  ylab('Асимметрия гуанина')
medG_realm = medG_realm + xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))
medG_realm = medG_realm + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
medG_realm

medG_tn = ggplot(data = df_mtdna, aes(x = Trophic_niche, y = med_G))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Трофические ниши')+
  ylab('Асимметрия гуанина')
medG_tn = medG_tn + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
medG_tn = medG_tn + xlim(c('Herbivore aquatic', 'Scavenger', 'Vertivore', 'Granivore', 'Herbivore terrestrial', 'Invertivore', 'Aquatic predator', 'Nectarivore', 'Omnivore', 'Frugivore'))
medG_tn

medT_tn = ggplot(data = df_mtdna, aes(x = Trophic_niche, y = med_T))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Трофические ниши')+
  ylab('Асимметрия тимина')
medT_tn = medT_tn + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
medT_tn = medT_tn + xlim(c('Herbivore aquatic', 'Scavenger', 'Vertivore', 'Granivore', 'Herbivore terrestrial', 'Invertivore', 'Aquatic predator', 'Nectarivore', 'Omnivore', 'Frugivore'))
medT_tn

df_pca = df_mtdna[c('species_name','gene_name','fAn', 'fGn', 'fCn', 'fTn', 'Stg_Sac','ghahSkew', 'chthSkew', 'med_T', 'med_G')]
gene_vector = c('fAn', 'fGn', 'fCn', 'fTn', 'Stg_Sac','ghahSkew', 'chthSkew', 'med_T', 'med_G')
gene_stats = data.frame(unique(df_pca$species_name))
for ( char in gene_vector){
  
  stats1 = aggregate(df_pca[,char], by = list(df_pca$species_name), FUN = 'sum')[2]
  stats1 = stats1/12
  gene_stats = cbind(gene_stats, stats1)
  
}
names(gene_stats) = c('species_name', gene_vector)
df_realm = df_mtdna[c('species_name', 'realm')]
gene_stats = merge(gene_stats, df_realm, by = 'species_name')
gene_stats = unique(gene_stats)
#row.names(gene_stats) = gene_stats$species_name
#gene_stats$species_name = NA
gene_stats = gene_stats[, colSums(is.na(gene_stats)) < nrow(gene_stats)]
stats_pca = prcomp(gene_stats[c(2,3,4,5,6,7,8,9,10)], center = TRUE, scale. = TRUE)
summary(stats_pca)

bipl = ggbiplot(stats_pca, x)
bipl

#Млеки и птицы
unzip("../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip", exdir = "../../Body/3Results/")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")) file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")
names(SynNuc)
SynNuc$ghahSkew = ((SynNuc$NeutralC - SynNuc$NeutralT))/((SynNuc$NeutralC + SynNuc$NeutralT))
SynNuc$chthSkew = ((SynNuc$NeutralA - SynNuc$NeutralG))/((SynNuc$NeutralA + SynNuc$NeutralG))
new_mam = SynNuc[, c(1, 2, 79, 80)]
new_mam = new_mam[new_mam$Gene != 'ND6',]
new_mam$class = 'Mammalia'
new_bird = df_mtdna[, c('species_name', 'gene_name', 'ghahSkew','chthSkew')]
new_bird$class = 'Aves'
new_bird$species_name = gsub(' ', '_', new_bird$species_name)
new_bird$gene_name[new_bird$gene_name == 'CYTB'] = 'CytB'
names(new_mam) = c('species_name', 'gene_name', 'ghahSkew', 'chthSkew', 'class')

new_big = rbind(new_mam, new_bird)
new_b_and_m = ggplot(new_big, aes(x = gene_name, y = ghahSkew, fill = class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Mitochondrial genes')+
  ylab('GhAhSkew')
new_b_and_m = new_b_and_m + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CytB","ND1","ND2"))
new_b_and_m = new_b_and_m + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
new_b_and_m

new_b_and_m_one = ggplot(new_big, aes(x = gene_name, y = chthSkew, fill = class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Mitochondrial genes')+
  ylab('ThChSkew')
new_b_and_m_one = new_b_and_m_one + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CytB","ND1","ND2"))
new_b_and_m_one = new_b_and_m_one + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
new_b_and_m_one
two_big = ggarrange(new_b_and_m, new_b_and_m_one,
                    ncol = 2,
                    nrow = 1)
two_big

for_article = SynNuc[, c(1, 2, 73, 74, 75, 76)]
for_article = for_article[for_article$Gene != 'ND6',]
for_article$class = 'Mammalia'
new_bird1 = df_mtdna[, c('species_name', 'gene_name', 'neutral_A','neutral_g', 'neutral_c', 'neutral_T')]
new_bird1$class = 'Aves'
new_bird1$species_name = gsub(' ', '_', new_bird1$species_name)
new_bird1$gene_name[new_bird1$gene_name == 'CYTB'] = 'CytB'
names(for_article) = c('species_name', 'gene_name', 'neutral_A', 'neutral_T', 'neutral_g', 'neutral_c', 'class')

new_big1 = rbind(for_article, new_bird1)
new_b_and_m1 = ggplot(new_big1, aes(x = gene_name, y = neutral_A, fill = class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Gene name')+
  ylab('Neutral T')+
  ylim(0, 200)
new_b_and_m1 = new_b_and_m1 + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CytB","ND1","ND2"))
new_b_and_m1 = new_b_and_m1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
new_b_and_m1

new_b_and_m2 = ggplot(new_big1, aes(x = gene_name, y = neutral_g, fill = class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Gene name')+
  ylab('Neutral C')+
  ylim(0, 200)
new_b_and_m2 = new_b_and_m2 + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CytB","ND1","ND2"))
new_b_and_m2 = new_b_and_m2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
new_b_and_m2

new_b_and_m3 = ggplot(new_big1, aes(x = gene_name, y = neutral_c, fill = class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Gene name')+
  ylab('Neutral G')+
  ylim(0, 200)
new_b_and_m3 = new_b_and_m3 + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CytB","ND1","ND2"))
new_b_and_m3 = new_b_and_m3 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
new_b_and_m3

new_b_and_m4 = ggplot(new_big1, aes(x = gene_name, y = neutral_T, fill = class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Gene name')+
  ylab('Neutral A')+
  ylim(0, 200)
new_b_and_m4 = new_b_and_m4 + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CytB","ND1","ND2"))
new_b_and_m4 = new_b_and_m4 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
new_b_and_m4

f10 = ggarrange(new_b_and_m1, new_b_and_m3, new_b_and_m2, new_b_and_m4, 
                ncol = 2, nrow = 2)
f10