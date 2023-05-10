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
  geom_boxplot(fill = 'blue')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB", "ND6", "ND1","ND2"))+
  ylim(0, 0.8)+
  annotate("text", x=11, y=0.78, label= "Th")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank())
graph1

graph2 = ggplot(data = nd6_correct, aes(x = gene_name, y = fCn))+
  geom_boxplot(fill = 'green')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB", "ND6", "ND1","ND2"))+
  ylim(0, 0.8)+
  annotate("text", x=11, y=0.78, label= "Ch")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank())
graph2

graph3 = ggplot(data = nd6_correct, aes(x = gene_name, y = fAn))+
  geom_boxplot(fill = 'yellow')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB", "ND6", "ND1","ND2"))+
  ylim(0, 0.8)+
  annotate("text", x=11, y=0.78, label= "Ah")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank())
graph3

graph4 = ggplot(data = nd6_correct, aes(x = gene_name, y = fGn))+
  geom_boxplot(fill = 'red')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB", "ND6", "ND1","ND2"))+
  ylim(0, 0.8)+
  annotate("text", x=11, y=0.78, label= "Gh")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank())
graph4

graph5 = ggplot(data = nd6_correct, aes(x = gene_name, y = GhAhSkew))+
  geom_boxplot(fill = 'orange')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB", "ND6", "ND1","ND2"))+
  ylim(-0.5,1)+
  annotate("text", x=11, y=-0.3, label= "GhAhSkew")+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
graph5

graph6 = ggplot(data = nd6_correct, aes(x = gene_name, y = ThChSkew))+
  geom_boxplot(fill = 'cyan')+
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
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Birds realms')+
  ylab('GhAhSkew')
skew_eco = skew_eco + xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))
skew_eco = skew_eco + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
skew_eco

skew_eco2 = ggplot(data = df_mtdna, aes(x = realm, y = chthSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Birds realms')+
  ylab('ThChSkew')
skew_eco2 = skew_eco2 + xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))
skew_eco2 = skew_eco2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
skew_eco2

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

#trophic niche
skew_eco31 = ggplot(data = df_mtdna, aes(x = Trophic_niche, y = chthSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Trophic niche')+
  ylab('ThChSkew')
skew_eco31 = skew_eco31 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
skew_eco31 = skew_eco31 + xlim(c('Herbivore aquatic', 'Scavenger', 'Vertivore', 'Granivore', 'Herbivore terrestrial', 'Invertivore', 'Aquatic predator', 'Nectarivore', 'Omnivore', 'Frugivore'))
skew_eco31

skew_eco12 = ggplot(data = df_mtdna, aes(x = Trophic_niche, y = ghahSkew))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  xlab('Trophic niche')+
  ylab('GhAhSkew')
skew_eco12 = skew_eco12 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
skew_eco12 = skew_eco12 + xlim(c('Herbivore aquatic', 'Scavenger', 'Vertivore', 'Granivore', 'Herbivore terrestrial', 'Invertivore', 'Aquatic predator', 'Nectarivore', 'Omnivore', 'Frugivore'))
skew_eco12

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
