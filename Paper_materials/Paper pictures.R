#Paper pictures
rm(list = ls(all=TRUE))
library(ggbiplot)
library(ggplot2)
library(ggpubr)
library(ggbiplot)
df_mtdna = read.csv('../Paper_materials/Birds_dataset_paper_extra.csv')

#Skews and frequencies with ND6

df_nd6 = read.csv('../Body/3Results/Birds_mtDNA_data.csv')
nd6_look = df_nd6[df_nd6$gene_name == 'ND6',]
nd6_nolook = df_nd6[df_nd6$gene_name != 'ND6',]
nd6_look$GhAhSkew = (nd6_look$neutral_g - nd6_look$neutral_A)/(nd6_look$neutral_g + nd6_look$neutral_A)
nd6_look$ThChSkew = (nd6_look$neutral_T - nd6_look$neutral_c)/(nd6_look$neutral_T + nd6_look$neutral_c)
nd6_look$fTn = nd6_look$neutral_T/nd6_look$neutral_amount
nd6_look$fAn = nd6_look$neutral_A/nd6_look$neutral_amount
nd6_look$fCn = nd6_look$neutral_c/nd6_look$neutral_amount
nd6_look$fGn = nd6_look$neutral_g/nd6_look$neutral_amount

nd6_nolook$GhAhSkew = (nd6_nolook$neutral_c- nd6_nolook$neutral_T)/(nd6_nolook$neutral_c + nd6_nolook$neutral_T)
nd6_nolook$ThChSkew = (nd6_nolook$neutral_A - nd6_nolook$neutral_g)/(nd6_nolook$neutral_A + nd6_nolook$neutral_g)
nd6_nolook$fTn = nd6_nolook$neutral_A/nd6_nolook$neutral_amount
nd6_nolook$fAn = nd6_nolook$neutral_T/nd6_nolook$neutral_amount
nd6_nolook$fCn = nd6_nolook$neutral_g/nd6_nolook$neutral_amount
nd6_nolook$fGn = nd6_nolook$neutral_c/nd6_nolook$neutral_amount

nd6_correct = rbind(nd6_look, nd6_nolook)

graph1 = ggplot(data = nd6_correct, aes(x = gene_name, y = fTn))+
  geom_boxplot(notch = TRUE, fill = 'blue')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB", "ND6", "ND1","ND2"))+
  ylim(0, 0.8)+
  annotate("text", x=11, y=0.78, label= "Th")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank())
graph1

graph2 = ggplot(data = nd6_correct, aes(x = gene_name, y = fCn))+
  geom_boxplot(notch = TRUE, fill = 'green')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB", "ND6", "ND1","ND2"))+
  ylim(0, 0.8)+
  annotate("text", x=11, y=0.78, label= "Ch")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank())
graph2

graph3 = ggplot(data = nd6_correct, aes(x = gene_name, y = fAn))+
  geom_boxplot(notch = TRUE, fill = 'yellow')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB", "ND6", "ND1","ND2"))+
  ylim(0, 0.8)+
  annotate("text", x=11, y=0.78, label= "Ah")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank())
graph3

graph4 = ggplot(data = nd6_correct, aes(x = gene_name, y = fGn))+
  geom_boxplot(notch = TRUE, fill = 'red')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB", "ND6", "ND1","ND2"))+
  ylim(0, 0.8)+
  annotate("text", x=11, y=0.78, label= "Gh")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank())
graph4

graph5 = ggplot(data = nd6_correct, aes(x = gene_name, y = GhAhSkew))+
  geom_boxplot(notch = TRUE, fill = 'orange')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB", "ND6", "ND1","ND2"))+
  ylim(-0.5,1)+
  annotate("text", x=11, y=-0.3, label= "GhAhSkew")+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
graph5

graph6 = ggplot(data = nd6_correct, aes(x = gene_name, y = ThChSkew))+
  geom_boxplot(notch = TRUE, fill = 'cyan')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB", "ND6", "ND1","ND2"))+
  ylim(-0.5,1)+
  annotate("text", x=11, y=-0.3, label= "ThChSkew")+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
graph6

nd6fr= ggarrange(graph1, graph4, graph2, graph3, graph6, graph5,
                 ncol = 2, nrow = 3)

nd6fr

#Mammalia against birds
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
  ylab('GhAhSkew')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CytB","ND1","ND2"))+
  annotate("text", x=8, y=-0.5, label= "GhAhSkew")+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
new_b_and_m

new_b_and_m_one = ggplot(new_big, aes(x = gene_name, y = chthSkew, fill = class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Mitochondrial genes')+
  ylab('ThChSkew')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CytB","ND1","ND2"))+
  annotate("text", x=8, y=0, label= "ThChSkew")+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
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
  ylim(0, 200)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CytB","ND1","ND2"))+
  annotate("text", x=8, y=175, label= "Th")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank())
new_b_and_m1

new_b_and_m2 = ggplot(new_big1, aes(x = gene_name, y = neutral_g, fill = class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Gene name')+
  ylim(0, 200)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CytB","ND1","ND2"))+
  annotate("text", x=8, y=150, label= "Ch")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank())
new_b_and_m2

new_b_and_m3 = ggplot(new_big1, aes(x = gene_name, y = neutral_c, fill = class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Gene name')+
  ylim(0, 200)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CytB","ND1","ND2"))+
  annotate("text", x=8, y=175, label= "Gh")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank())
new_b_and_m3

new_b_and_m4 = ggplot(new_big1, aes(x = gene_name, y = neutral_T, fill = class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Gene name')+
  ylim(0, 200)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CytB","ND1","ND2"))+
  annotate("text", x=8, y=150, label= "Ah")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank())
new_b_and_m4

f10 = ggarrange(new_b_and_m1, new_b_and_m3, new_b_and_m2, new_b_and_m4, new_b_and_m, new_b_and_m_one,
                ncol = 2, nrow = 3)
f10

#Ecology
#realm
skew_eco = ggplot(data = df_mtdna, aes(x = realm, y = ghahSkew))+
  geom_boxplot(outlier.shape = NA, notch = T, fill = 'blue')+
  xlab('Birds realms')+
  xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))+
  annotate("text", x=7, y=-0.25, label= "GhAhSkew")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank())
skew_eco

skew_eco2 = ggplot(data = df_mtdna, aes(x = realm, y = chthSkew))+
  geom_boxplot(outlier.shape = NA, notch = T, fill = 'red')+
  xlab('Birds realms')+
  xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))+
  annotate("text", x=7, y=0.35, label= "ThChSkew")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank())
skew_eco2

medT_realm = ggplot(data = df_mtdna, aes(x = realm, y = med_T))+
  geom_boxplot(outlier.shape = NA, notch = T, fill = 'yellow')+
  xlab('Birds realm')+
  xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))+
  annotate("text", x=7, y=0.825, label= "T asymmetry")+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
medT_realm

medG_realm = ggplot(data = df_mtdna, aes(x = realm, y = med_G))+
  geom_boxplot(outlier.shape = NA, notch = T, fill = 'green')+
  xlab('Birds realm')+
  xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))+
  annotate("text", x=7, y=0.6, label= "G asymmetry")+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
medG_realm

f_eco = ggarrange(skew_eco, skew_eco2, medG_realm, medT_realm,
                ncol = 2, nrow = 2)
f_eco


#trophic niche
skew_eco31 = ggplot(data = df_mtdna, aes(x = Trophic_niche, y = chthSkew))+
  geom_boxplot(outlier.shape = NA, notch = T, fill = 'red')+
  xlab('Trophic niche')+
  ylab('ThChSkew')+
  xlim(c('Herbivore aquatic', 'Scavenger', 'Vertivore', 'Granivore', 'Herbivore terrestrial', 'Invertivore', 'Aquatic predator', 'Nectarivore', 'Omnivore', 'Frugivore'))+
  annotate("text", x=7, y=0.32, label= "ThChSkew")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank())
skew_eco31

skew_eco12 = ggplot(data = df_mtdna, aes(x = Trophic_niche, y = ghahSkew))+
  geom_boxplot(outlier.shape = NA, notch = T, fill = 'blue')+
  xlab('Trophic niche')+
  ylab('GhAhSkew')+
  xlim(c('Herbivore aquatic', 'Scavenger', 'Vertivore', 'Granivore', 'Herbivore terrestrial', 'Invertivore', 'Aquatic predator', 'Nectarivore', 'Omnivore', 'Frugivore'))+
  annotate("text", x=7, y=-0.25, label= "GhAhSkew")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank())
skew_eco12

medG_tn = ggplot(data = df_mtdna, aes(x = Trophic_niche, y = med_G))+
  geom_boxplot(outlier.shape = NA, notch = T, fill = 'green')+
  xlab('Trophic niche')+
  ylab('Guanine asymmetry')+
  xlim(c('Herbivore aquatic', 'Scavenger', 'Vertivore', 'Granivore', 'Herbivore terrestrial', 'Invertivore', 'Aquatic predator', 'Nectarivore', 'Omnivore', 'Frugivore'))+
  annotate("text", x=7, y=0.6, label= "G asymmetry")+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
medG_tn

medT_tn = ggplot(data = df_mtdna, aes(x = Trophic_niche, y = med_T))+
  geom_boxplot(outlier.shape = NA, notch = T, fill = 'yellow')+
  xlab('Trophic niche')+
  ylab('Thymine asymmetry')+
  xlim(c('Herbivore aquatic', 'Scavenger', 'Vertivore', 'Granivore', 'Herbivore terrestrial', 'Invertivore', 'Aquatic predator', 'Nectarivore', 'Omnivore', 'Frugivore'))+
  annotate("text", x=7, y=0.825, label= "T asymmetry")+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
medT_tn

f_tn = ggarrange(skew_eco12, skew_eco31, medG_tn, medT_tn,
                  ncol = 2, nrow = 2)
f_tn


#PCA mtDNA metrics
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
bipl_com = ggbiplot(stats_pca)
bipl_com

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

#PCA density
g5 = ggplot(gene_stats, aes(x=PC1, color=realm)) +
  geom_density(size = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (48.0%)')+
  scale_colour_manual(name="Origin", values= c("black", "red", "black", "black", "black", "black", "black", "black", "black"))+
  theme(legend.position="none")

g5

g6 = ggplot(gene_stats, aes(x=PC2, color=realm)) +
  geom_density(size = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC2 (32.2%)')+
  scale_colour_manual(name="Origin", values= c("black", "red", "black", "black", "black", "black", "black", "black", "black"))+
  theme(axis.title.y=element_blank(), axis.ticks.y=element_blank())

g6

g7 = ggplot(gene_stats, aes(x=PC1, color=Trophic_niche)) +
  geom_density(size = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (48.0%)')+
  scale_colour_manual(name="Origin", values= c("black", "black", "black", "red", "black", "black", "black", "black", "black", 'black'))+
  theme(legend.position="none")

g7
g8 = ggplot(gene_stats, aes(x=PC2, color=Trophic_niche)) +
  geom_density(size = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC2 (32.2%)')+
  scale_colour_manual(name="Origin", values= c("black", "black", "black", "red", "black", "black", "black", "black", "black", 'black'))+
  theme(axis.title.y=element_blank(), axis.ticks.y=element_blank())
g8

f_den_stats1 = ggarrange(g5, g6,
                        ncol = 2, nrow = 1)
f_den_stats1

f_den_stats2 = ggarrange(g7, g8,
                         ncol = 2, nrow = 1)
f_den_stats2
