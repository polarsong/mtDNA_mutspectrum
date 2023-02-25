rm(list = ls(all=TRUE))
install.packages("ggbiplot")
install.packages("ggplot2")
install.packages("ggpubr")
library(ggbiplot)
library(ggplot2)
library(ggpubr)
df_mtdna = read.csv('../../Head/2Scripts/Birds_dataset_paper.csv')
f1 = ggplot(data = df_mtdna, aes(x = gene_name, y = fTn))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Mitochondrial genes')+
  ylab('Thymine frequence')+
  ylim(0, 0.7)
f1 = f1 +xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND1","ND2"))
f1 = f1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
f2 = ggplot(data = df_mtdna, aes(x = gene_name, y = fCn))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Mitochondrial genes')+
  ylab('Cytosine frequence')+
  ylim(0, 0.7)
f2 = f2 +xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND1","ND2"))
f2 = f2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
f3 = ggplot(data = df_mtdna, aes(x = gene_name, y = fGn))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Mitochondrial genes')+
  ylab('Guanine frequence')+
  ylim(0, 0.7)
f3 = f3 +xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND1","ND2"))
f3 = f3 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
f4 = ggplot(data = df_mtdna, aes(x = gene_name, y = fAn))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Mitochondrial genes')+
  ylab('Adenine frequence')+
  ylim(0, 0.7)
f4 = f4 +xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND1","ND2"))
f4 = f4 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
f5 = ggarrange(f1, f3, f2, f4, 
               ncol = 2, nrow = 2)
f5
f51 = ggarrange(f1,f2,
                 ncol = 1, nrow = 2)
f51
f52 = ggarrange(f3,f4,
                ncol = 1, nrow = 2)
f52


skew_all = ggplot(data = df_mtdna, aes(x = gene_name, y = ghahSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Gene names')+
  ylab('GhAhSkew')
skew_all = skew_all + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND1","ND2"))
skew_all

skew_all1 = ggplot(data = df_mtdna, aes(x = gene_name, y = chthSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Gene names')+
  ylab('ThChSkew')
skew_all1 = skew_all1 + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND1","ND2"))
skew_all1

stg_all = ggplot(data = df_mtdna, aes(x = gene_name, y = Stg_Sac))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Gene names')+
  ylab('Stg-Sac')
stg_all = stg_all + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND1","ND2"))
stg_all

skew_eco = ggplot(data = df_mtdna, aes(x = realm, y = ghahSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Birds realms')+
  ylab('GhAhSkew')
skew_eco = skew_eco + xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))
skew_eco = skew_eco + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
skew_eco

skew_eco1 = ggplot(data = df_mtdna, aes(x = Trophic_level, y = ghahSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Trophic level')+
  ylab('GhAhSkew')
skew_eco1 = skew_eco1 + xlim(c('Carnivore', 'Omnivore', 'Herbivore', 'Scavenger'))
skew_eco1

skew_eco12 = ggplot(data = df_mtdna, aes(x = Trophic_niche, y = ghahSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Trophic niche')+
  ylab('GhAhSkew')
skew_eco12 = skew_eco12 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
skew_eco12 = skew_eco12 + xlim(c('Herbivore aquatic', 'Scavenger', 'Vertivore', 'Granivore', 'Herbivore terrestrial', 'Invertivore', 'Aquatic predator', 'Nectarivore', 'Omnivore', 'Frugivore'))
skew_eco12

skew_eco2 = ggplot(data = df_mtdna, aes(x = realm, y = chthSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Birds realms')+
  ylab('ThChSkew')
skew_eco2 = skew_eco2 + xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))
skew_eco2 = skew_eco2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
skew_eco2

skew_eco3 = ggplot(data = df_mtdna, aes(x = Trophic_level, y = chthSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Trophic level')+
  ylab('ThChSkew')
skew_eco3 = skew_eco3 + xlim(c('Carnivore', 'Omnivore', 'Herbivore', 'Scavenger'))
skew_eco3

skew_eco31 = ggplot(data = df_mtdna, aes(x = Trophic_niche, y = chthSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Trophic niche')+
  ylab('ThChSkew')
skew_eco31 = skew_eco31 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
skew_eco31 = skew_eco31 + xlim(c('Herbivore aquatic', 'Scavenger', 'Vertivore', 'Granivore', 'Herbivore terrestrial', 'Invertivore', 'Aquatic predator', 'Nectarivore', 'Omnivore', 'Frugivore'))
skew_eco31

stg_eco = ggplot(data = df_mtdna, aes(x = realm, y = Stg_Sac))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Birds realm')+
  ylab('Stg-Sac')
stg_eco = stg_eco + xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))
stg_eco = stg_eco + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
stg_eco

stg_eco1 = ggplot(data = df_mtdna, aes(x = Trophic_level, y = Stg_Sac))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Trophic level')+
  ylab('Stg-Sac')
stg_eco1 = stg_eco1 + xlim(c('Carnivore', 'Omnivore', 'Herbivore', 'Scavenger'))
stg_eco1

stg_eco2 = ggplot(data = df_mtdna, aes(x = Trophic_niche, y = Stg_Sac))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Trophic niche')+
  ylab('Stg-Sac')
stg_eco2 = stg_eco2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
stg_eco2 = stg_eco2 + xlim(c('Herbivore aquatic', 'Scavenger', 'Vertivore', 'Granivore', 'Herbivore terrestrial', 'Invertivore', 'Aquatic predator', 'Nectarivore', 'Omnivore', 'Frugivore'))
stg_eco2

medT_realm = ggplot(data = df_mtdna, aes(x = realm, y = med_T))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Birds realm')+
  ylab('Thymine asymmetry')
medT_realm = medT_realm + xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))
medT_realm = medT_realm + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
medT_realm

medG_realm = ggplot(data = df_mtdna, aes(x = realm, y = med_G))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Birds realm')+
  ylab('Guanine asymmetry')
medG_realm = medG_realm + xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))
medG_realm = medG_realm + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
medG_realm

medG_tl = ggplot(data = df_mtdna, aes(x = Trophic_level, y = med_G))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Trophic level')+
  ylab('Guanine asymmetry')
medG_tl = medG_tl + xlim(c('Carnivore', 'Omnivore', 'Herbivore', 'Scavenger'))
medG_tl

medT_tl = ggplot(data = df_mtdna, aes(x = Trophic_level, y = med_T))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Trophic level')+
  ylab('Thymine asymmetry')
medT_tl = medT_tl + xlim(c('Carnivore', 'Omnivore', 'Herbivore', 'Scavenger'))
medT_tl

medG_tn = ggplot(data = df_mtdna, aes(x = Trophic_niche, y = med_G))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Trophic niche')+
  ylab('Guanine asymmetry')
medG_tn = medG_tn + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
medG_tn = medG_tn + xlim(c('Herbivore aquatic', 'Scavenger', 'Vertivore', 'Granivore', 'Herbivore terrestrial', 'Invertivore', 'Aquatic predator', 'Nectarivore', 'Omnivore', 'Frugivore'))
medG_tn

medT_tn = ggplot(data = df_mtdna, aes(x = Trophic_niche, y = med_T))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Trophic niche')+
  ylab('Thymine asymmetry')
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

bipl = ggbiplot(stats_pca)
bipl

#Doing statistics
df_antarctic = df_mtdna[df_mtdna$realm == 'Antarctic',]
df_else = df_mtdna[df_mtdna$realm != 'Antarctic',]
df_oceania = df_mtdna[df_mtdna$realm == 'Oceania',]
wilcox.test(df_antarctic$ghahSkew, df_else$ghahSkew)
wilcox.test(df_antarctic$ghahSkew, df_oceania$ghahSkew)
wilcox.test(df_antarctic$chthSkew, df_else$chthSkew)
wilcox.test(df_antarctic$chthSkew, df_oceania$chthSkew)
wilcox.test(df_antarctic$Stg_Sac, df_else$Stg_Sac)
wilcox.test(df_antarctic$Stg_Sac, df_oceania$Stg_Sac)
df_ha = df_mtdna[df_mtdna$Trophic_niche == 'Herbivore aquatic',]
df_noha = df_mtdna[df_mtdna$Trophic_niche != 'Herbivore aquatic',]
wilcox.test(df_ha$ghahSkew, df_noha$ghahSkew)
wilcox.test(df_ha$chthSkew, df_noha$chthSkew)

#PCA coloring 
library('ggbiplot')
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

bipl = ggbiplot(stats_pca, groups = gene_stats$realm, labels = gene_stats$species_name, labels.size = 2)+
  scale_colour_manual(name="Origin", values= c("black", "green", "black", "black", "black", "black", "black", "red", "black"))
bipl
ggplotly(bipl)

bipl1 = ggbiplot(stats_pca, choices = 2:3, groups = gene_stats$realm, labels = gene_stats$species_name, labels.size = 2)+
  scale_colour_manual(name="Origin", values= c("black", "green", "black", "black", "black", "black", "black", "red", "black"))
bipl1
ggplotly(bipl1)


bipl3 = ggbiplot(stats_pca, choices = 3:4, groups = gene_stats$realm, labels = gene_stats$species_name, labels.size = 2)+
  scale_colour_manual(name="Origin", values= c("black", "green", "black", "black", "black", "black", "black", "red", "black"))
bipl3
ggplotly(bipl3)

#install.packages('plotly')
#install.packages('dplyr')
#install.packages("car")
#install.packages("babynames")
#install.packages("gapminder")

library(plotly)
library(dplyr)
#library(carData)
library(gapminder)
library(babynames)

bipl_niche = ggbiplot(stats_pca, groups = gene_stats$Trophic_niche, labels = gene_stats$species_name, labels.size = 2)+
  scale_colour_manual(name='Origin', values = c('black', 'black', 'black', 'red', 'black','black', 'black', 'black', 'black', 'black'))
bipl_niche
ggplotly(bipl_niche)

ggplotly(bipl)
ggplotly(bipl2)
ggplotly(bipl3)
spec_pca = ggplotly(bipl1)

direct.label(bipl)
birds_pca = data.frame(stats_pca$x)
birds_pca = birds_pca[,c(1,2)]
birds_pca$species_name = row.names(birds_pca)
gene_stats$species_name = row.names(gene_stats)
gene_stats = merge(gene_stats, birds_pca, by = 'species_name')
row.names(gene_stats) = gene_stats$species_name
gene_stats = gene_stats[,c(2:13)]

#Lera's plot
library(ggfortify)
library(dabestr)
library(ggrepel)


g1 = ggplot(gene_stats, aes(x=PC1, color=realm)) +
  geom_density(size = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (44.02%)')

g1
g2 = ggplot(gene_stats, aes(x=PC2, color=realm)) +
  geom_density(size = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC2 (35.28%)')
  

g2


#Phylo PCA
install.packages("phytools")
library(ape)
library(geiger)
library(caper)
library(phytools)
tree = read.tree("../../Body/3Results/phylo.treefile")
df_tree = read.csv('../../Head/2Scripts/Table_for_PGLS.csv')
row.names(df_tree) = df_tree$species_name
phy=multi2di(tree)
name.check(phy, df_tree)
df_tree = df_tree[df_tree$species_name != "Agapornis_pullarius",]
df_tree = df_tree[df_tree$species_name != "Mergus_squamatus",]
name.check(phy, df_tree)
df_gene_realms = df_tree[,c(1:3,53:58)]
df_realms = unique(as.data.frame(df_mtdna[,c('species_name', 'realm')]))
df_realms = df_realms[df_realms$species_name != "Agapornis pullarius",]
df_realms = df_realms[df_realms$species_name != "Mergus squamatus",]
df_realms$species_name = gsub(" ", "_", df_realms$species_name)
df_gene_realms = merge(df_realms, df_gene_realms, by = 'species_name')
row.names(df_gene_realms) = df_gene_realms$species_name
name.check(df_gene_realms, phy)


df_try_pca = df_gene_realms[,c(3:10)]
phylo_pca = phyl.pca(phy, df_try_pca)
phylo_pca
par(mar = c(4.1,4.1,2.1,1.1),las=1)
plot(phylo_pca, main="")

par(cex.axis=0.8,mar=c(5.1,5.1,1.1,1.1))
phylomorphospace(phy,
                 scores(phylo_pca)[,1:2],
                 ftype="off",node.size=c(0,1),bty="n",las=1,
                 xlab="PC1 ()",
                 ylab=expression(paste("PC2 ("%up%"lamellae number, "
                                       %down%"tail length)")))
row.names(df_realms) = df_realms$species_name
df_realms1 = df_realms[,c(2)]

eco<-setNames(df_realms[,2],rownames(df_realms))
levels(eco)
ECO=to.matrix(eco)
tiplabels(pie=ECO[phy$tip.label,],cex=0.5)
legend(x="bottomright",legend=levels(eco),cex=0.8,pch=21,
       pt.bg=rainbow(n=length(levels(eco))),pt.cex=1.5)



#What birds are in realms?
unique(unique(df_mtdna[df_mtdna$realm == 'Antarctic',])$species_name)
unique(unique(df_mtdna[df_mtdna$Trophic_niche == 'Herbivore aquatic',])$species_name)

#Dima's dat
df_all = read.table('../../Head/2Scripts/VertebratePolymorphisms.MutSpecData.txt')
unique(df_all$Class)
df_mammals = filter(df_all, df_all$Class == 'Mammalia')
unique(df_mammals$Species)
df_what = filter(df_all, df_all$Species == 'Abrothrix_longipilis')
unique(df_mammals$Class)
install.packages('plyr')
library(plyr)
a = (count(df_mammals, 'Species') >= 12)
df_mammals$ghahSkew = ((df_mammals$FrC - df_mammals$FrT))/((df_mammals$FrC + df_mammals$FrT))
mam_box = ggplot(df_mammals, aes(x = Gene, y = ghahSkew))+
  geom_boxplot()
mam_box = mam_box + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","CytB","ND1","ND2"))
mam_box
bird_only = df_mtdna[c('species_name', 'gene_name', 'ghahSkew')]
bird_only$class = 'Aves'
mam_only = df_mammals[c('Species', 'Gene', 'ghahSkew')]
names(mam_only) = c('species_name', 'gene_name', 'ghahSkew')
mam_only$class = 'Mammalia'
bird_only = bird_only[bird_only$gene_name != 'ND5',]
bird_only$gene_name[bird_only$gene_name == 'CYTB'] = 'CytB'
unique(bird_only$gene_name)
mam_and_birds = rbind(bird_only, mam_only)
b_and_m = ggplot(mam_and_birds, aes(x = gene_name, y = ghahSkew, fill = class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)
b_and_m = b_and_m + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","CytB","ND1","ND2"))
b_and_m
df_stats = data.frame()
for (i in c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","CytB","ND1","ND2"))
{
  a = wilcox.test(bird_only[bird_only$gene_name == i,]$ghahSkew, mam_only[mam_only$gene_name == i,]$ghahSkew)
  b = c(i, a$p.value)
  df_stats = rbind(df_stats, b)
}
names(df_stats) = c('gene_name', 'p_value')
#a = wilcox.test(bird_only[bird_only$gene_name == 'CytB',]$ghahSkew, mam_only[mam_only$gene_name == 'CytB',]$ghahSkew)
write.csv(df_stats, file = 'ghahSkew_Aves_against_Mammalia_p_values.csv')




#Birds sizes + 3D
install.packages("plot3D")
library("plot3D")
bm = unique(df_mtdna[c(2,108,109,110,111,112,113,114,115,116,117)])
b_3d = merge(bm, gene_stats, by = 'species_name')
x = b_3d$Mass
y = b_3d$fAn
z = b_3d$fGn
scatter3D(x, y, z, colvar = y, theta = 15, phi = 20)
stats_pca_big = prcomp(b_3d[c(2:11)], center = TRUE, scale. = TRUE)
summary(stats_pca_big)

b_3d_bipl = ggbiplot(stats_pca_big, groups = b_3d$Trophic_niche, labels = b_3d$species_name, labels.size = 2)
b_3d_bipl
b_3d_bipl2 = ggbiplot(stats_pca_big, groups = b_3d$realm, labels = b_3d$species_name, labels.size = 2)
ggplotly(b_3d_bipl)
ggplotly(b_3d_bipl2)

#Alya's dataset
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
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)
new_b_and_m = new_b_and_m + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CytB","ND1","ND2"))
new_b_and_m

new_b_and_m_one = ggplot(new_big, aes(x = gene_name, y = chthSkew, fill = class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)
new_b_and_m_one = new_b_and_m_one + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CytB","ND1","ND2"))
new_b_and_m_one

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

#Work with Valya's data
valya_data = read.csv('../../Head/2Scripts/valyadata_final.csv')
valya_data = na.omit(valya_data)
gene_data = df_mtdna[,c('species_name', 'ghahSkew', 'chthSkew', 'fAn', 'fGn', 'fCn', 'fTn', 'Stg_Sac', 'med_G', 'med_T')]
gene_data$species_name = gsub(' ', '_', gene_data$species_name)
valya_data["species_name"][valya_data["species_name"] == "Strigops_habroptilus"] = "Strigops_habroptila"
valya_gene = merge(gene_data, valya_data, by = 'species_name')

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
gene_stats$species_name = gsub(' ', '_', gene_stats$species_name)
gene_stats = merge(gene_stats, valya_data, by = 'species_name')
#gene_stats$species_name = NA
gene_stats = gene_stats[, colSums(is.na(gene_stats)) < nrow(gene_stats)]
stats_pca = prcomp(gene_stats[c(2,3,4,5,6,7,8,9,10)], center = TRUE, scale. = TRUE)
summary(stats_pca)
bipl = ggbiplot(stats_pca, groups = gene_stats$reference2, labels = gene_stats$species_name, labels.size = 2)
bipl
ggplotly(bipl)


#Valya's data boxplots
fly_box1 = ggplot(valya_gene, aes(x = flying, y = ghahSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Ability to fly')+
  ylab('GhAhSkew')
fly_box1

fly_box2 = ggplot(valya_gene, aes(x = wintering, y = ghahSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Wintering')+
  ylab('GhAhSkew')
fly_box2

fly_box3 = ggplot(valya_gene, aes(x = diving, y = ghahSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Diving')+
  ylab('GhAhSkew')
fly_box3

fly_box4 = ggplot(valya_gene, aes(x = far_migration, y = ghahSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Far_migration')+
  ylab('GhAhSkew')
fly_box4

fly_box5 = ggplot(valya_gene, aes(x = flying, y = chthSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Ability to fly')+
  ylab('ChThSkew')
fly_box5

fly_box6 = ggplot(valya_gene, aes(x = wintering, y = chthSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Wintering')+
  ylab('ChThSkew')
fly_box6

fly_box7 = ggplot(valya_gene, aes(x = diving, y = chthSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Diving')+
  ylab('ChThSkew')
fly_box7

fly_box8 = ggplot(valya_gene, aes(x = far_migration, y = chthSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Far_migration')+
  ylab('ChThSkew')
fly_box8

med1 = ggplot(valya_gene, aes(x = flying, y = med_G))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Ability to fly')+
  ylab('med_G')
med1

med2 = ggplot(valya_gene, aes(x = wintering, y = med_G))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Wintering')+
  ylab('med_G')
med2

med3 = ggplot(valya_gene, aes(x = diving, y = med_G))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Divers')+
  ylab('med_G')
med3

med4 = ggplot(valya_gene, aes(x = far_migration, y = med_G))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Far migration')+
  ylab('med_G')
med4

med5 = ggplot(valya_gene, aes(x = flying, y = med_T))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Ability to fly')+
  ylab('med_T')
med5

med6 = ggplot(valya_gene, aes(x = wintering, y = med_T))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Wintering')+
  ylab('med_T')
med6

med7 = ggplot(valya_gene, aes(x = diving, y = med_T))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Divers')+
  ylab('med_T')
med7

med8 = ggplot(valya_gene, aes(x = far_migration, y = med_T))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Far migration')+
  ylab('med_T')
med8

#Aminoacids shift

df_sgc = df_mtdna[,c(2, 6, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53,
               54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77,
               78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95)] #getting codon usage


anti_vec_all = c('TTC','TTT','TCC','TCT','TAC','TAT','TGC','TGT',
            'TTA','TTG','TCA','TCG','TAA','TAG','TGA','TGG',
            'CTC','CTT','CCC','CCT','CAC','CAT','CGC','CGT',
            'CTA','CTG','CCA','CCG','CAA','CAG','CGA','CGG',
            'ATC','ATT','ACC','ACT','AAC','AAT','AGC','AGT',
            'ATA','ATG','ACA','ACG','AAA','AAG','AGA','AGG',
            'GTC','GTT','GCC','GCT','GAC','GAT','GGC','GGT',
            'GTA','GTG','GCA','GCG','GAA','GAG','GGA','GGG')

vec_all = c('TTT','TTC','TCT','TCC','TAT','TAC','TGT','TGC',
            'TTG','TTA','TCG','TCA','TAG','TAA','TGG','TGA',
            'CTT','CTC','CCT','CCC','CAT','CAC','CGT','CGC',
            'CTG','CTA','CCG','CCA','CAG','CAA','CGG','CGA',
            'ATT','ATC','ACT','ACC','AAT','AAC','AGT','AGC',
            'ATG','ATA','ACG','ACA','AAG','AAA','AGG','AGA',
            'GTT','GTC','GCT','GCC','GAT','GAC','GGT','GGC',
            'GTG','GTA','GCG','GCA','GAG','GAA','GGG','GGA')

df_aa = data.frame()

for (i in 1:nrow(df_sgc)){
  org_gen = df_sgc[i,]
  org_gen = as.vector(org_gen)
  df_out= data.frame(df_sgc[i,]$species_name, df_sgc[i,]$gene_name) 
  for (codon in seq(from = 1, to = 64)){
    if (as.numeric(as.character(unlist(org_gen[vec_all[codon]]))) + as.numeric(as.character(unlist(org_gen[anti_vec_all[codon]]))) == 0)
    {
      norm_cod = 0
      df_out = cbind(df_out, norm_cod)
    }
    else
    {
      norm_cod = (as.numeric(as.character(unlist(org_gen[vec_all[codon]]))))/(as.numeric(as.character(unlist(org_gen[vec_all[codon]]))) + as.numeric(as.character(unlist(org_gen[anti_vec_all[codon]]))))
      df_out = cbind(df_out, norm_cod)
    }
  }
  df_aa = rbind(df_aa, df_out)
}

names(df_aa) = c('species_name', 'gene_name', vec_all)

write.csv(df_aa, file = 'Aminoacids_shift_birds.csv')
