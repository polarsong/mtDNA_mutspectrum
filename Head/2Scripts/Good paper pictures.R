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
  ylab('Thymine frequence')
f1 = f1 +xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND1","ND2"))
f1 = f1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
f2 = ggplot(data = df_mtdna, aes(x = gene_name, y = fCn))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Mitochondrial genes')+
  ylab('Cytosine frequence')
f2 = f2 +xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND1","ND2"))
f2 = f2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
f3 = ggplot(data = df_mtdna, aes(x = gene_name, y = fGn))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Mitochondrial genes')+
  ylab('Guanine frequence')
f3 = f3 +xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND1","ND2"))
f3 = f3 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
f4 = ggplot(data = df_mtdna, aes(x = gene_name, y = fAn))+
  geom_boxplot(outlier.shape = NA)+
  xlab('Mitochondrial genes')+
  ylab('Adenine frequence')
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
stats_pca = prcomp(gene_stats[c(1,2,3,4,5,6,7,8)], center = TRUE, scale. = TRUE)
summary(stats_pca)

bipl = ggbiplot(stats_pca, groups = gene_stats[gene_stats$realm == 'Antarctic',])
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
df_realm = df_mtdna[c('species_name', 'realm')]
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

install.packages('plotly')
install.packages('dplyr')
#install.packages("car")
install.packages("babynames")
install.packages("gapminder")

library(plotly)
library(dplyr)
#library(carData)
library(gapminder)
library(babynames)

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
pca_plot = autoplot(stats_pca, data = gene_stats, colour = 'gray', loadings = TRUE, loadings.label = TRUE, loadings.label.size = 5, scale = 0, 
                    loadings.colour = 'black', loadings.label.colour = 'black')+
  geom_point(data = gene_stats, aes(PC1, PC2, color = realm, alpha = 0.5),size = 3)
pca_plot

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

#Dima's data
df_all = read.table('../../Head/2Scripts/VertebratePolymorphisms.MutSpecData.txt')
df_mammals = filter(df_all, df_all$Class == 'Mammalia')
df_what = filter(df_all, df_all$Species == 'Abrothrix_longipilis')
unique(df_mammals$Class)
