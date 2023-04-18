#Картинки для дипома/практики
rm(list = ls(all=TRUE))
library(ggbiplot)
library(ggplot2)
library(ggpubr)
library(plotly)
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
df_for_diploma = df_mtdna
df_for_diploma$russian_eco = ''
df_for_diploma$russian_tn = ''
dfr1 = df_for_diploma[df_for_diploma$realm == 'Antarctic',]
dfr2 = df_for_diploma[df_for_diploma$realm == 'Nearctic',]
dfr3 = df_for_diploma[df_for_diploma$realm == 'Palearctic',]
dfr4 = df_for_diploma[df_for_diploma$realm == 'Indo_Malay',]
dfr5 = df_for_diploma[df_for_diploma$realm == 'Afrotropic',]
dfr6 = df_for_diploma[df_for_diploma$realm == 'Neotropic',]
dfr7 = df_for_diploma[df_for_diploma$realm == 'Madagascar',]
dfr8 = df_for_diploma[df_for_diploma$realm == 'Australian',]
dfr9 = df_for_diploma[df_for_diploma$realm == 'Oceania',]

dfr1$russian_eco = 'Антарктика'
dfr2$russian_eco = 'Неарктика'
dfr3$russian_eco = 'Палеарктика'
dfr4$russian_eco = 'Индо-Малайзия'
dfr5$russian_eco = 'Афротропики'
dfr6$russian_eco = 'Неотропики'
dfr7$russian_eco = 'Мадагаскар'
dfr8$russian_eco = 'Австралия'
dfr9$russian_eco = 'Океания'

df_for_diploma = rbind(dfr1, dfr2, dfr3, dfr4, dfr5, dfr6, dfr7, dfr8, dfr9)

dfr1 = df_for_diploma[df_for_diploma$Trophic_niche == 'Herbivore aquatic',]
dfr2 = df_for_diploma[df_for_diploma$Trophic_niche == 'Scavenger',]
dfr3 = df_for_diploma[df_for_diploma$Trophic_niche == 'Vertivore',]
dfr4 = df_for_diploma[df_for_diploma$Trophic_niche == 'Granivore',]
dfr5 = df_for_diploma[df_for_diploma$Trophic_niche == 'Herbivore terrestrial',]
dfr6 = df_for_diploma[df_for_diploma$Trophic_niche == 'Invertivore',]
dfr7 = df_for_diploma[df_for_diploma$Trophic_niche == 'Aquatic predator',]
dfr8 = df_for_diploma[df_for_diploma$Trophic_niche == 'Nectarivore',]
dfr9 = df_for_diploma[df_for_diploma$Trophic_niche == 'Omnivore',]
dfr10 = df_for_diploma[df_for_diploma$Trophic_niche == 'Frugivore',]


dfr1$russian_tn = 'Водные травоядные'
dfr2$russian_tn = 'Падальщики'
dfr3$russian_tn = 'Позвоночные'
dfr4$russian_tn = 'Зерноядные'
dfr5$russian_tn = 'Наземные травоядные'
dfr6$russian_tn = 'Беспозвоночные'
dfr7$russian_tn = 'Водные хищники'
dfr8$russian_tn = 'Нектароядные'
dfr9$russian_tn = 'Всеядные'
dfr10$russian_tn = 'Плодоядные'

df_for_diploma = rbind(dfr1, dfr2, dfr3, dfr4, dfr5, dfr6, dfr7, dfr8, dfr9, dfr10)


skew_eco = ggplot(data = df_for_diploma, aes(x = russian_eco, y = ghahSkew))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Экозоны')+
  ylab('GhAhSkew')
skew_eco = skew_eco + xlim(c('Антарктика', 'Неарктика', 'Палеарктика', 'Индо-Малайзия', 'Афротропики', 'Мадагаскар', 'Неотропики', 'Австралия', 'Океания'))
skew_eco = skew_eco + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
skew_eco


skew_eco12 = ggplot(data = df_for_diploma, aes(x = russian_tn, y = ghahSkew))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Трофические ниши')+
  ylab('GhAhSkew')
skew_eco12 = skew_eco12 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
skew_eco12 = skew_eco12 + xlim(c('Водные травоядные', 'Падальщики', 'Позвоночные', 'Зерноядные', 'Наземные травоядные', 'Беспозвоночные', 'Водные хищники', 'Нектароядные', 'Всеядные', 'Плодоядные'))
skew_eco12


skew_eco2 = ggplot(data = df_for_diploma, aes(x = russian_eco, y = chthSkew))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Экозоны')+
  ylab('ThChSkew')
skew_eco2 = skew_eco2 + xlim(c('Антарктика', 'Неарктика', 'Палеарктика', 'Индо-Малайзия', 'Афротропики', 'Мадагаскар', 'Неотропики', 'Австралия', 'Океания'))
skew_eco2 = skew_eco2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
skew_eco2

skew_eco31 = ggplot(data = df_for_diploma, aes(x = russian_tn, y = chthSkew))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Трофические ниши')+
  ylab('ThChSkew')
skew_eco31 = skew_eco31 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
skew_eco31 = skew_eco31 + xlim(c('Водные травоядные', 'Падальщики', 'Позвоночные', 'Зерноядные', 'Наземные травоядные', 'Беспозвоночные', 'Водные хищники', 'Нектароядные', 'Всеядные', 'Плодоядные'))
skew_eco31

stg_eco = ggplot(data = df_for_diploma, aes(x = russian_eco, y = Stg_Sac))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Экозоны')+
  ylab('Stg-Sac')
stg_eco = stg_eco + xlim(c('Антарктика', 'Неарктика', 'Палеарктика', 'Индо-Малайзия', 'Афротропики', 'Мадагаскар', 'Неотропики', 'Австралия', 'Океания'))
stg_eco = stg_eco + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
stg_eco


stg_eco2 = ggplot(data = df_for_diploma, aes(x = russian_tn, y = Stg_Sac))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Трофические ниши')+
  ylab('Stg-Sac')
stg_eco2 = stg_eco2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
stg_eco2 = stg_eco2 + xlim(c('Водные травоядные', 'Падальщики', 'Позвоночные', 'Зерноядные', 'Наземные травоядные', 'Беспозвоночные', 'Водные хищники', 'Нектароядные', 'Всеядные', 'Плодоядные'))
stg_eco2


medT_realm = ggplot(data = df_for_diploma, aes(x = russian_eco, y = med_T))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Экозоны')+
  ylab('Асимметрия тимина')
medT_realm = medT_realm + xlim(c('Антарктика', 'Неарктика', 'Палеарктика', 'Индо-Малайзия', 'Афротропики', 'Мадагаскар', 'Неотропики', 'Австралия', 'Океания'))
medT_realm = medT_realm + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
medT_realm

medG_realm = ggplot(data = df_for_diploma, aes(x = russian_eco, y = med_G))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Экозоны')+
  ylab('Асимметрия гуанина')
medG_realm = medG_realm + xlim(c('Антарктика', 'Неарктика', 'Палеарктика', 'Индо-Малайзия', 'Афротропики', 'Мадагаскар', 'Неотропики', 'Австралия', 'Океания'))
medG_realm = medG_realm + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
medG_realm

medG_tn = ggplot(data = df_for_diploma, aes(x = russian_tn, y = med_G))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Трофические ниши')+
  ylab('Асимметрия гуанина')
medG_tn = medG_tn + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
medG_tn = medG_tn + xlim(c('Водные травоядные', 'Падальщики', 'Позвоночные', 'Зерноядные', 'Наземные травоядные', 'Беспозвоночные', 'Водные хищники', 'Нектароядные', 'Всеядные', 'Плодоядные'))
medG_tn

medT_tn = ggplot(data = df_for_diploma, aes(x = russian_tn, y = med_T))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Трофические ниши')+
  ylab('Асимметрия тимина')
medT_tn = medT_tn + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
medT_tn = medT_tn + xlim(c('Водные травоядные', 'Падальщики', 'Позвоночные', 'Зерноядные', 'Наземные травоядные', 'Беспозвоночные', 'Водные хищники', 'Нектароядные', 'Всеядные', 'Плодоядные'))
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



#Млеки и птицы
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip", exdir = "../../Body/3Results/")
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
  xlab('Митохондриальные гены')+
  ylab('GhAhSkew')+
  scale_fill_hue(labels = c("Птицы", "Млекопитающие"))
new_b_and_m = new_b_and_m + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CytB","ND1","ND2"))
new_b_and_m = new_b_and_m + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
new_b_and_m

new_b_and_m_one = ggplot(new_big, aes(x = gene_name, y = chthSkew, fill = class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Митохондриальные гены')+
  ylab('ThChSkew')+
  scale_fill_hue(labels = c("Птицы", "Млекопитающие"))
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
  xlab('Митохондриальные гены')+
  ylab('Тимин')+
  ylim(0, 200)+
  scale_fill_hue(labels = c("Птицы", "Млекопитающие"))
new_b_and_m1 = new_b_and_m1 + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CytB","ND1","ND2"))
new_b_and_m1 = new_b_and_m1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
new_b_and_m1

new_b_and_m2 = ggplot(new_big1, aes(x = gene_name, y = neutral_g, fill = class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Митохондриальные гены')+
  ylab('Цитозин')+
  ylim(0, 200)+
  scale_fill_hue(labels = c("Птицы", "Млекопитающие"))
new_b_and_m2 = new_b_and_m2 + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CytB","ND1","ND2"))
new_b_and_m2 = new_b_and_m2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
new_b_and_m2

new_b_and_m3 = ggplot(new_big1, aes(x = gene_name, y = neutral_c, fill = class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Митохондриальные гены')+
  ylab('Гуанин')+
  ylim(0, 200)+
  scale_fill_hue(labels = c("Птицы", "Млекопитающие"))
new_b_and_m3 = new_b_and_m3 + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CytB","ND1","ND2"))
new_b_and_m3 = new_b_and_m3 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
new_b_and_m3

new_b_and_m4 = ggplot(new_big1, aes(x = gene_name, y = neutral_T, fill = class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Митохондриальные гены')+
  ylab('Аденин')+
  ylim(0, 200)+
  scale_fill_hue(labels = c("Птицы", "Млекопитающие"))
new_b_and_m4 = new_b_and_m4 + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CytB","ND1","ND2"))
new_b_and_m4 = new_b_and_m4 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
new_b_and_m4

f10 = ggarrange(new_b_and_m1, new_b_and_m3, new_b_and_m2, new_b_and_m4, 
                ncol = 2, nrow = 2)
f10


#Density для первого PCA
df_pca = df_mtdna[c('species_name','gene_name','fAn', 'fGn', 'fCn', 'fTn', 'Stg_Sac','ghahSkew', 'chthSkew', 'med_T', 'med_G')]
gene_vector = c('fAn', 'fGn', 'fCn', 'fTn', 'Stg_Sac','ghahSkew', 'chthSkew', 'med_T', 'med_G')
gene_stats = data.frame(unique(df_pca$species_name))
for ( char in gene_vector){
  
  stats1 = aggregate(df_pca[,char], by = list(df_pca$species_name), FUN = 'sum')[2]
  stats1 = stats1/12
  gene_stats = cbind(gene_stats, stats1)
  
}
names(gene_stats) = c('species_name', gene_vector)
df_realm = df_for_diploma[c('species_name', 'russian_eco', 'russian_tn')]
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

g5 = ggplot(gene_stats, aes(x=PC1, color=russian_eco)) +
  geom_density()+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (48.0%)')+
  scale_colour_manual(name="Origin", values= c("black", "red", "black", "black", "black", "black", "black", "black", "black"))+
  ylab('Плотность')


g5

g6 = ggplot(gene_stats, aes(x=PC2, color=russian_eco)) +
  geom_density()+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC2 (32.2%)')+
  scale_colour_manual(name="Origin", values= c("black", "red", "black", "black", "black", "black", "black", "black", "black"))+
  ylab('Плотность')

g6

g7 = ggplot(gene_stats, aes(x=PC1, color=russian_tn)) +
  geom_density()+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (48.0%)')+
  scale_colour_manual(name="Origin", values= c("black", "red", "black", "black", "black", "black", "black", "black", "black", 'black'))+
  ylab('Плотность')

g7
g8 = ggplot(gene_stats, aes(x=PC2, color=russian_tn)) +
  geom_density()+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC2 (32.2%)')+
  scale_colour_manual(name="Origin", values= c("black", "red", "black", "black", "black", "black", "black", "black", "black", 'black'))


g8



#MutSpec
mut_data = read.table("C:/Users/User/Desktop/Birds mutspec results from Bogdan/mutspec12.tsv", header = TRUE, fill = TRUE)
mut_data_syn = mut_data[mut_data$Label == 'syn',]
mut_data_syn = mut_data_syn[,c(1,2,3,4,5,7,8)]
ecozone_data = df_for_diploma[,c('species_name', 'russian_eco', 'russian_tn')]
ecozone_data = unique(ecozone_data)
ecozone_data$species_name = gsub(' ', '_', ecozone_data$species_name)
mut_data_syn = mut_data_syn[!grepl('Node', mut_data_syn$AltNode),]
names(mut_data_syn) = c('Mut', 'ObsNum', 'ExpNum', 'RawMutSpec', 'MutSpec', 'species_name', 'Label') 
mut_data_syn = merge(mut_data_syn, ecozone_data, by = 'species_name')
mut_data_syn$MutSpec =  as.numeric(mut_data_syn$MutSpec)

AC = mut_data_syn[mut_data_syn$Mut == 'A>C',] 
AG = mut_data_syn[mut_data_syn$Mut == 'A>G',]
AT = mut_data_syn[mut_data_syn$Mut == 'A>T',]
GC = mut_data_syn[mut_data_syn$Mut == 'G>C',]
GT = mut_data_syn[mut_data_syn$Mut == 'G>T',]
GA = mut_data_syn[mut_data_syn$Mut == 'G>A',]
CG = mut_data_syn[mut_data_syn$Mut == 'C>G',]
CT = mut_data_syn[mut_data_syn$Mut == 'C>T',]
CA = mut_data_syn[mut_data_syn$Mut == 'C>A',]
TG = mut_data_syn[mut_data_syn$Mut == 'T>G',]
TC = mut_data_syn[mut_data_syn$Mut == 'T>C',]
TA = mut_data_syn[mut_data_syn$Mut == 'T>A',]

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

mut_data_syn1 = rbind(AC, AG, AT, GC, GT, GA, CT, CA, CG, TA, TG, TC)

ggplot(data = mut_data_syn1, aes(x = Mut, y = MutSpec)) +
  geom_boxplot()+
  ylab('Мутационный спектр')

ggplot(data = mut_data_syn1, aes(x = Mut, y = MutSpec, fill = realm)) +
  geom_boxplot()+
  ylab('Мутационный спектр')

ggplot(data = mut_data_syn1, aes(x = Mut, y = MutSpec, fill = Trophic_niche)) +
  geom_boxplot()+
  ylab('Мутационный спектр')

