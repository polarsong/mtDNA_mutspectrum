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

#Stats for TBST
nd6_correct$tbst = ''
gene1 = nd6_correct[nd6_correct$gene_name == 'COX1',]
gene1$tbst = 1
gene2 = nd6_correct[nd6_correct$gene_name == 'COX2',]
gene2$tbst = 2
gene3 = nd6_correct[nd6_correct$gene_name == 'ATP8',]
gene3$tbst = 3
gene4 = nd6_correct[nd6_correct$gene_name == 'ATP6',]
gene4$tbst = 4
gene5 = nd6_correct[nd6_correct$gene_name == 'COX3',]
gene5$tbst = 5
gene6 = nd6_correct[nd6_correct$gene_name == 'ND3',]
gene6$tbst = 6
gene7 = nd6_correct[nd6_correct$gene_name == 'ND4L',]
gene7$tbst = 7
gene8 = nd6_correct[nd6_correct$gene_name == 'ND4',]
gene8$tbst = 8
gene9 = nd6_correct[nd6_correct$gene_name == 'ND5',]
gene9$tbst = 9
gene10 = nd6_correct[nd6_correct$gene_name == 'CYTB',]
gene10$tbst = 10

gene_tbst = rbind(gene1, gene2, gene3, gene4, gene5, gene6, gene7, gene8, gene9, gene10)
wilcox.test(gene_tbst$GhAhSkew, gene_tbst$tbst)
wilcox.test(gene_tbst$ThChSkew, gene_tbst$tbst)


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
  geom_density(linewidth = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (48.0%)')+
  scale_colour_manual(name="Origin", values= c("black", "red", "black", "black", "black", "black", "black", "black", "black"))+
  theme(legend.position="none")

g5

g6 = ggplot(gene_stats, aes(x=PC2, color=realm)) +
  geom_density(linewidth = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC2 (32.2%)')+
  scale_colour_manual(name="Origin", values= c("black", "red", "black", "black", "black", "black", "black", "black", "black"))+
  theme(axis.title.y=element_blank(), axis.ticks.y=element_blank())

g6

g7 = ggplot(gene_stats, aes(x=PC1, color=Trophic_niche)) +
  geom_density(linewidth = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (48.0%)')+
  scale_colour_manual(name="Origin", values= c("black", "black", "black", "red", "black", "black", "black", "black", "black", 'black'))+
  theme(legend.position="none")

g7
g8 = ggplot(gene_stats, aes(x=PC2, color=Trophic_niche)) +
  geom_density(linewidth = 1)+
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

#MutSpec graphs
#Boxplots
df_mtdna = read.csv('../Head/2Scripts/Birds_dataset_paper.csv')
mut_data = read.table("C:/Users/User/Desktop/Birds mutspec results from Bogdan/mutspec12.tsv", header = TRUE, fill = TRUE)
mut_data_ff = mut_data[mut_data$Label == 'syn',]
mut_data_ff = mut_data_ff[,c(1,2,3,4,5,7,8)]
ecozone_data = df_mtdna[,c('species_name', 'realm', 'Trophic_niche')]
ecozone_data = unique(ecozone_data)
ecozone_data$species_name = gsub(' ', '_', ecozone_data$species_name)
mut_data_ff = mut_data_ff[!grepl('Node', mut_data_ff$AltNode),]
names(mut_data_ff) = c('Mut', 'ObsNum', 'ExpNum', 'RawMutSpec', 'MutSpec', 'species_name', 'Label') 
mut_data_ff = merge(mut_data_ff, ecozone_data, by = 'species_name')
mut_data_ff$MutSpec =  as.numeric(mut_data_ff$MutSpec)

AC = mut_data_ff[mut_data_ff$Mut == 'A>C',] 
AG = mut_data_ff[mut_data_ff$Mut == 'A>G',]
AT = mut_data_ff[mut_data_ff$Mut == 'A>T',]
GC = mut_data_ff[mut_data_ff$Mut == 'G>C',]
GT = mut_data_ff[mut_data_ff$Mut == 'G>T',]
GA = mut_data_ff[mut_data_ff$Mut == 'G>A',]
CG = mut_data_ff[mut_data_ff$Mut == 'C>G',]
CT = mut_data_ff[mut_data_ff$Mut == 'C>T',]
CA = mut_data_ff[mut_data_ff$Mut == 'C>A',]
TG = mut_data_ff[mut_data_ff$Mut == 'T>G',]
TC = mut_data_ff[mut_data_ff$Mut == 'T>C',]
TA = mut_data_ff[mut_data_ff$Mut == 'T>A',]

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

mut_data_ff1 = rbind(AC, AG, AT, GC, GT, GA, CT, CA, CG, TA, TG, TC)

ggplot(mut_data_ff1, aes(x = Mut, y = MutSpec))+
  geom_boxplot(notch = TRUE)+
  xlab('Mutations')+
  ylab('Mutational spectrum')

extra_mut = mut_data_ff1
extra_mut$realm <- factor(extra_mut$realm, levels = c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))

ggplot(data = extra_mut, aes(x = Mut, y = MutSpec, fill = realm)) +
  geom_boxplot(notch = TRUE)+
  guides(fill = guide_legend(title = "Ecozone"))+
  xlab('Mutations')+
  ylab('Mutational spectrum')


ggplot(data = extra_mut, aes(x = Mut, y = MutSpec, fill = Trophic_niche)) +
  geom_boxplot(notch = TRUE)+
  guides(fill = guide_legend(title = "Trophic niche"))+
  xlab('Mutations')+
  ylab('Mutational spectrum')

#PCA biplot and density
library(tidyr)
library(data.table)

pca_data = mut_data_ff1[,c(1,5,9)]
ex = reshape(data = pca_data, idvar = 'species_name',
             timevar = 'Mut',
             direction = 'wide')
names(ex) = c('species_name', 'T>G', 'T>C', 'T>A', 'C>G', 'C>A', 'C>T', 'G>A','G>T','G>C', 'A>T', 'A>C', 'A>G')
ex = merge(ex, ecozone_data, by = 'species_name')
ex = ex[,c(1,14,15,12,13,11,6,5,7,8,10,9,4,3,2)]
row.names(ex) = ex$species_name


stats_pca1 = prcomp(ex[,c(4,5,6,7,8,9,10,11,12,13,14,15)], center = TRUE, scale. = TRUE)
summary(stats_pca1)

bipl = ggbiplot(stats_pca1, varname.size = 4, varname.adjust = 6)
bipl

birds_ms_pca = data.frame(stats_pca1$x)
birds_ms_pca = birds_ms_pca[,c(1,2)]
birds_ms_pca$species_name = row.names(birds_ms_pca)
ex2 = merge(ex, birds_ms_pca, by = 'species_name')

d1 = ggplot(ex2, aes(x=PC1, color=realm)) +
  geom_density(linewidth = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (38.2%)')+
  scale_colour_manual(name="Origin", values= c("black", "red", "black", "black", "black", "black", "black", "black", "black"))+
  theme(legend.position="none")

d1

d2 = ggplot(ex2, aes(x=PC2, color=realm)) +
  geom_density(linewidth = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC2 (10.4%)')+
  scale_colour_manual(name="Origin", values= c("black", "red", "black", "black", "black", "black", "black", "black", "black"))+
  theme(axis.title.y=element_blank(), axis.ticks.y=element_blank())

d2

d3 = ggplot(ex2, aes(x=PC1, color=Trophic_niche)) +
  geom_density(linewidth = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (38.2%)')+
  scale_colour_manual(name="Origin", values= c("black", "black", "black", "red", "black", "black", "black", "black", "black", 'black'))+
  theme(legend.position="none")

d3

d4 = ggplot(ex2, aes(x=PC2, color=Trophic_niche)) +
  geom_density(linewidth = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC2 (10.4%)')+
  scale_colour_manual(name="Origin", values= c("black", "black", "black", "red", "black", "black", "black", "black", "black", 'black'))+
  theme(axis.title.y=element_blank(), axis.ticks.y=element_blank())
d4

d5 = ggarrange(d1, d2,
                 ncol = 2, nrow = 1)
d5

d6 = ggarrange(d3, d4,
               ncol = 2, nrow = 1)
d6

#Transitions/Transversions
pca_data = mut_data_ff1[,c(1,5,7,8,9)]

pca_data_shaped = dcast.data.table(setDT(pca_data), species_name+realm+Trophic_niche~Mut,
                                   value.var='MutSpec')
pca_data_shaped$TR_TS = (pca_data_shaped$`A>G`+pca_data_shaped$`C>T`+pca_data_shaped$`G>A`+pca_data_shaped$`T>C`)/(pca_data_shaped$`A>C`+pca_data_shaped$`A>T`+pca_data_shaped$`C>A`+pca_data_shaped$`C>G`+pca_data_shaped$`G>C`+pca_data_shaped$`G>T`+pca_data_shaped$`T>A`+pca_data_shaped$`T>G`)
pca_data_shaped$CT_GA = (pca_data_shaped$`C>T`/pca_data_shaped$`G>A`)
pca_data_shaped$AG_TC = (pca_data_shaped$`A>G`/pca_data_shaped$`T>C`)

realm_tr_ts = ggplot(data = pca_data_shaped, aes(x = realm, y = TR_TS))+
  geom_boxplot(fill = 'cyan', notch = TRUE)+
  ylim(0,75)+
  annotate("text", x=6, y=50, label= "Tr/Ts")+
  xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
realm_tr_ts

niche_tr_ts = ggplot(data = pca_data_shaped, aes(x = Trophic_niche, y = TR_TS))+
  geom_boxplot(fill = 'orange', notch = TRUE)+
  ylim(0,75)+
  annotate("text", x=6, y=50, label= "Tr/Ts")+
  xlim(c('Herbivore aquatic', 'Scavenger', 'Vertivore', 'Granivore', 'Herbivore terrestrial', 'Invertivore', 'Aquatic predator', 'Nectarivore', 'Omnivore', 'Frugivore'))+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
niche_tr_ts

realm_tr_ts1 = ggplot(data = pca_data_shaped, aes(x = realm, y = CT_GA))+
  geom_boxplot(fill = 'green', notch = TRUE)+
  ylim(0,5)+
  annotate("text", x=6, y=3, label= "CT/GA")+
  xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank())
realm_tr_ts1

niche_tr_ts1 = ggplot(data = pca_data_shaped, aes(x = Trophic_niche, y = CT_GA))+
  geom_boxplot(fill = 'yellow', notch = TRUE)+
  ylim(0,5)+
  annotate("text", x=2, y=1, label= "CT/GA")+
  xlim(c('Herbivore aquatic', 'Scavenger', 'Vertivore', 'Granivore', 'Herbivore terrestrial', 'Invertivore', 'Aquatic predator', 'Nectarivore', 'Omnivore', 'Frugivore'))+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank())
niche_tr_ts1

realm_tr_ts2 = ggplot(data = pca_data_shaped, aes(x = realm, y = AG_TC))+
  geom_boxplot(fill = 'blue', notch = TRUE)+
  ylim(0,5)+
  annotate("text", x=6, y=1, label= "AG/TC")+
  xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank())
realm_tr_ts2

niche_tr_ts2 = ggplot(data = pca_data_shaped, aes(x = Trophic_niche, y = AG_TC))+
  geom_boxplot(fill = 'red', notch = TRUE)+
  ylim(0,5)+
  annotate("text", x=2, y=1, label= "AG/TC")+
  xlim(c('Herbivore aquatic', 'Scavenger', 'Vertivore', 'Granivore', 'Herbivore terrestrial', 'Invertivore', 'Aquatic predator', 'Nectarivore', 'Omnivore', 'Frugivore'))+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank())
niche_tr_ts2

tr_ts_graph= ggarrange(realm_tr_ts2, niche_tr_ts2, realm_tr_ts1, niche_tr_ts1, realm_tr_ts, niche_tr_ts,
                 ncol = 2, nrow = 3)

tr_ts_graph
