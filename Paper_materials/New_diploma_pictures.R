rm(list = ls(all=TRUE))
library(ggbiplot)
library(ggplot2)
library(ggpubr)
df_mtdna = read.csv('../Paper_materials/Birds_dataset_paper_extra.csv')
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

graph1 = ggplot(data = df_mtdna, aes(x = gene_name, y = fTn))+
  geom_boxplot(notch = TRUE, fill = 'blue')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB", "ND1","ND2"))+
  ylim(0, 0.8)+
  annotate("text", x=11, y=0.78, label= "Th")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank())
graph1

graph2 = ggplot(data = df_mtdna, aes(x = gene_name, y = fCn))+
  geom_boxplot(notch = TRUE, fill = 'green')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB", "ND1","ND2"))+
  ylim(0, 0.8)+
  annotate("text", x=11, y=0.78, label= "Ch")+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
graph2

graph3 = ggplot(data = df_mtdna, aes(x = gene_name, y = fAn))+
  geom_boxplot(notch = TRUE, fill = 'yellow')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB", "ND1","ND2"))+
  ylim(0, 0.8)+
  annotate("text", x=11, y=0.78, label= "Ah")+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
graph3

graph4 = ggplot(data = df_mtdna, aes(x = gene_name, y = fGn))+
  geom_boxplot(notch = TRUE, fill = 'red')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB", "ND1","ND2"))+
  ylim(0, 0.8)+
  annotate("text", x=11, y=0.78, label= "Gh")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank())
graph4

graph5 = ggplot(data = df_mtdna, aes(x = gene_name, y = ghahSkew))+
  geom_boxplot(notch = TRUE, fill = 'orange')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB", "ND1","ND2"))+
  ylim(-0.5,1)+
  annotate("text", x=11, y=-0.3, label= "GhAhSkew")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank())
graph5

graph6 = ggplot(data = df_mtdna, aes(x = gene_name, y = chthSkew))+
  geom_boxplot(notch = TRUE, fill = 'cyan')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB", "ND1","ND2"))+
  ylim(-0.5,1)+
  annotate("text", x=11, y=-0.3, label= "ThChSkew")+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
graph6

graph7 = ggplot(data = df_mtdna, aes(x = gene_name, y = Stg_Sac))+
  geom_boxplot(notch = TRUE, fill = 'purple')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB", "ND1","ND2"))+
  ylim(-0.5,1)+
  annotate("text", x=11, y=-0.3, label= "Stg-Sac")+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
graph7

nd6fr= ggarrange(graph1, graph4, graph5, graph2, graph3, graph6,
                 ncol = 3, nrow = 2)

nd6fr

nd6fr1 = ggarrange(graph1, graph4, graph2, graph3,
                    ncol = 2, nrow = 2)
nd6fr1

unzip("../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip", exdir = "../../Body/3Results/")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")) file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")
names(SynNuc)
SynNuc$ghahSkew = ((SynNuc$NeutralC - SynNuc$NeutralT))/((SynNuc$NeutralC + SynNuc$NeutralT))
SynNuc$chthSkew = ((SynNuc$NeutralA - SynNuc$NeutralG))/((SynNuc$NeutralA + SynNuc$NeutralG))
new_mam = SynNuc[, c(1, 2, 79, 80)]
new_mam = new_mam[new_mam$Gene != 'ND6',]
new_mam$class = 'Млекопитающие'
new_bird = df_mtdna[, c('species_name', 'gene_name', 'ghahSkew','chthSkew')]
new_bird$class = 'Птицы'
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
        axis.title.y=element_blank(), axis.ticks.y=element_blank(),
        legend.position = "none")
new_b_and_m

new_b_and_m_one = ggplot(new_big, aes(x = gene_name, y = chthSkew, fill = class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Mitochondrial genes')+
  ylab('ThChSkew')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CytB","ND1","ND2"))+
  annotate("text", x=8, y=0, label= "ThChSkew")+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())+
  guides(fill = guide_legend(title = "Класс позвоночных"))
new_b_and_m_one

for (i in unique(new_big$gene_name))
{
  print(i)
  print(wilcox.test(new_big[new_big$gene_name == i & new_big$class == 'Птицы',]$ghahSkew, 
                      new_big[new_big$gene_name == i & new_big$class == 'Млекопитающие',]$ghahSkew))
}

for (i in unique(new_big$gene_name))
{
  print(i)
  print(wilcox.test(new_big[new_big$gene_name == i & new_big$class == 'Птицы',]$chthSkew, 
                    new_big[new_big$gene_name == i & new_big$class == 'Млекопитающие',]$chthSkew))
}




etr = ggarrange(new_b_and_m, new_b_and_m_one,
                nrow = 1, ncol = 2)

etr

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
        axis.title.y=element_blank(),axis.ticks.y=element_blank(),
        legend.position = "none")
new_b_and_m1

new_b_and_m2 = ggplot(new_big1, aes(x = gene_name, y = neutral_g, fill = class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Gene name')+
  ylim(0, 200)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CytB","ND1","ND2"))+
  annotate("text", x=8, y=150, label= "Ch")+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank(),
        legend.position = "none")
new_b_and_m2

new_b_and_m3 = ggplot(new_big1, aes(x = gene_name, y = neutral_c, fill = class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Gene name')+
  ylim(0, 200)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CytB","ND1","ND2"))+
  annotate("text", x=8, y=175, label= "Gh")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank(),
        legend.position = "none")
new_b_and_m3

new_b_and_m4 = ggplot(new_big1, aes(x = gene_name, y = neutral_T, fill = class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Gene name')+
  ylim(0, 200)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CytB","ND1","ND2"))+
  annotate("text", x=8, y=150, label= "Ah")+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank(),
        legend.position = "none")
new_b_and_m4

f10 = ggarrange(new_b_and_m1, new_b_and_m3, new_b_and_m, new_b_and_m2, new_b_and_m4, new_b_and_m_one,
                ncol = 3, nrow = 2)
f10


#ecology
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
  geom_boxplot(outlier.shape = NA, notch = T, fill = 'blue')+
  xlab('Birds realms')+
  xlim(c('Антарктика', 'Неарктика', 'Палеарктика', 'Индо-Малайзия', 'Афротропики', 'Мадагаскар', 'Неотропики', 'Австралия', 'Океания'))+
  annotate("text", x=7, y=-0.25, label= "GhAhSkew")+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
skew_eco

skew_eco2 = ggplot(data = df_for_diploma, aes(x = russian_eco, y = chthSkew))+
  geom_boxplot(outlier.shape = NA, notch = T, fill = 'red')+
  xlab('Birds realms')+
  xlim(c('Антарктика', 'Неарктика', 'Палеарктика', 'Индо-Малайзия', 'Афротропики', 'Мадагаскар', 'Неотропики', 'Австралия', 'Океания'))+
  annotate("text", x=7, y=0.35, label= "ThChSkew")+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
skew_eco2

for (i in unique(df_for_diploma$russian_eco))
{
  if (i != 'Антарктика')
  {
    print(i)
    print(wilcox.test(df_for_diploma[df_for_diploma$russian_eco == 'Антарктика',]$ghahSkew, 
                      df_for_diploma[df_for_diploma$russian_eco == i,]$ghahSkew))
  }
}

for (i in unique(df_for_diploma$russian_eco))
{
  if (i != 'Антарктика')
  {
    print(i)
    print(wilcox.test(df_for_diploma[df_for_diploma$russian_eco == 'Антарктика',]$chthSkew, 
                      df_for_diploma[df_for_diploma$russian_eco == i,]$chthSkew))
  }
}


medT_realm = ggplot(data = df_for_diploma, aes(x = russian_eco, y = med_T))+
  geom_boxplot(outlier.shape = NA, notch = T, fill = 'yellow')+
  xlab('Birds realm')+
  xlim(c('Антарктика', 'Неарктика', 'Палеарктика', 'Индо-Малайзия', 'Афротропики', 'Мадагаскар', 'Неотропики', 'Австралия', 'Океания'))+
  annotate("text", x=7, y=0.825, label= "Асимметрия Т")+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
medT_realm

medG_realm = ggplot(data = df_for_diploma, aes(x = russian_eco, y = med_G))+
  geom_boxplot(outlier.shape = NA, notch = T, fill = 'green')+
  xlab('Birds realm')+
  xlim(c('Антарктика', 'Неарктика', 'Палеарктика', 'Индо-Малайзия', 'Афротропики', 'Мадагаскар', 'Неотропики', 'Австралия', 'Океания'))+
  annotate("text", x=7, y=0.6, label= "Асимметрия Г")+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
medG_realm

f_eco = ggarrange(skew_eco, skew_eco2,
                  ncol = 2, nrow = 1)
f_eco

skew_eco31 = ggplot(data = df_for_diploma, aes(x = russian_tn, y = chthSkew))+
  geom_boxplot(outlier.shape = NA, notch = T, fill = 'red')+
  xlab('Trophic niche')+
  ylab('ThChSkew')+
  xlim(c('Водные травоядные', 'Падальщики', 'Позвоночные', 'Зерноядные', 'Наземные травоядные', 'Беспозвоночные', 'Водные хищники', 'Нектароядные', 'Всеядные', 'Плодоядные'))+
  annotate("text", x=7, y=0.32, label= "ThChSkew")+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
skew_eco31

skew_eco12 = ggplot(data = df_for_diploma, aes(x = russian_tn, y = ghahSkew))+
  geom_boxplot(outlier.shape = NA, notch = T, fill = 'blue')+
  xlab('Trophic niche')+
  ylab('GhAhSkew')+
  xlim(c('Водные травоядные', 'Падальщики', 'Позвоночные', 'Зерноядные', 'Наземные травоядные', 'Беспозвоночные', 'Водные хищники', 'Нектароядные', 'Всеядные', 'Плодоядные'))+
  annotate("text", x=7, y=-0.25, label= "GhAhSkew")+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
skew_eco12

for (i in unique(df_for_diploma$russian_tn))
{
  if (i != 'Водные травоядные')
  {
    print(i)
    print(wilcox.test(df_for_diploma[df_for_diploma$russian_tn == 'Водные травоядные',]$ghahSkew, 
                      df_for_diploma[df_for_diploma$russian_tn == i,]$ghahSkew))
  }
}

for (i in unique(df_for_diploma$russian_tn))
{
  if (i != 'Водные травоядные')
  {
    print(i)
    print(wilcox.test(df_for_diploma[df_for_diploma$russian_tn == 'Водные травоядные',]$chthSkew, 
                      df_for_diploma[df_for_diploma$russian_tn == i,]$chthSkew))
  }
}
check1 = df_for_diploma[df_for_diploma$russian_tn == 'Падальщики',]

medG_tn = ggplot(data = df_for_diploma, aes(x = russian_tn, y = med_G))+
  geom_boxplot(outlier.shape = NA, notch = T, fill = 'green')+
  xlab('Trophic niche')+
  ylab('Guanine asymmetry')+
  xlim(c('Водные травоядные', 'Падальщики', 'Позвоночные', 'Зерноядные', 'Наземные травоядные', 'Беспозвоночные', 'Водные хищники', 'Нектароядные', 'Всеядные', 'Плодоядные'))+
  annotate("text", x=7, y=0.6, label= "Асимметрия Г")+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
medG_tn

medT_tn = ggplot(data = df_for_diploma, aes(x = russian_tn, y = med_T))+
  geom_boxplot(outlier.shape = NA, notch = T, fill = 'yellow')+
  xlab('Trophic niche')+
  ylab('Thymine asymmetry')+
  xlim(c('Водные травоядные', 'Падальщики', 'Позвоночные', 'Зерноядные', 'Наземные травоядные', 'Беспозвоночные', 'Водные хищники', 'Нектароядные', 'Всеядные', 'Плодоядные'))+
  annotate("text", x=7, y=0.825, label= "Асимметрия Т")+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
medT_tn

f_tn = ggarrange(skew_eco12, skew_eco31,
                 ncol = 2, nrow = 1)
f_tn


#PCA 
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
  geom_density(linewidth = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (48.0%)')+
  scale_colour_manual(name="Origin", values= c("black", "red", "black", "black", "black", "black", "black", "black", "black"))+
  ylab('Плотность')+
  theme(legend.position="none")

g5

g6 = ggplot(gene_stats, aes(x=PC2, color=russian_eco)) +
  geom_density(linewidth = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC2 (32.2%)')+
  scale_colour_manual(name="Origin", values= c("black", "red", "black", "black", "black", "black", "black", "black", "black"))+
  ylab('Плотность')+
  theme(axis.title.y=element_blank(), axis.ticks.y=element_blank())+
  guides(color = guide_legend(title = "Экозона"))


g6

g7 = ggplot(gene_stats, aes(x=PC1, color=russian_tn)) +
  geom_density(linewidth = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (48.0%)')+
  scale_colour_manual(name="Origin", values= c("black", "red", "black", "black", "black", "black", "black", "black", "black", 'black'))+
  ylab('Плотность')+
  theme(legend.position="none")

g7
g8 = ggplot(gene_stats, aes(x=PC2, color=russian_tn)) +
  geom_density(linewidth = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC2 (32.2%)')+
  ylab('Плотность')+
  scale_colour_manual(name="Origin", values= c("black", "red", "black", "black", "black", "black", "black", "black", "black", 'black'))+
  theme(axis.title.y=element_blank(), axis.ticks.y=element_blank())+
  guides(color = guide_legend(title = "Трофическая ниша"))



g8

f_den_stats1 = ggarrange(g5, g6, g7, g8,
                         ncol = 2, nrow = 2)
f_den_stats1




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
  ylab('Мутационный спектр')+
  xlab('Мутации')

wilcox.test(mut_data_syn1[mut_data_syn1$Mut == 'C>T',]$MutSpec, mut_data_syn1[mut_data_syn1$Mut == 'G>A',]$MutSpec)

wilcox.test(mut_data_syn1[mut_data_syn1$Mut == 'A>G',]$MutSpec ,mut_data_syn1[mut_data_syn1$Mut == 'T>C',]$MutSpec)

MutSpecAG = mut_data_syn1[mut_data_syn1$Mut == 'A>G' | mut_data_syn1$Mut == 'C>T',]
ggplot(data = MutSpecAG, aes(x = Mut, y = MutSpec, fill = realm)) +
  geom_boxplot()+
  ylab('Мутационный спектр')+
  xlab('Мутации')+
  guides(fill = guide_legend(title = "Экозона"))


ggplot(data = MutSpecAG, aes(x = Mut, y = MutSpec, fill = Trophic_niche)) +
  geom_boxplot()+
  ylab('Мутационный спектр')+
  xlab('Мутации')+
  guides(fill = guide_legend(title = "Трофическая ниша"))

MutSpecAGR = MutSpecAG[MutSpecAG$realm == 'Антарктика' & MutSpecAG$Mut == 'A>G',]

MutSpecAG[MutSpecAG$realm == 'Антарктика' & MutSpecAG$Mut == 'A>G',]$RawMutSpec

for (i in unique(MutSpecAG$realm))
{
  if (i != 'Антарктика')
  {
    print(i)
    print(wilcox.test(MutSpecAG[MutSpecAG$realm == 'Антарктика' & MutSpecAG$Mut == 'A>G',]$MutSpec, 
               MutSpecAG[MutSpecAG$realm == i & MutSpecAG$Mut == 'A>G',]$MutSpec))
  }
}

for (i in unique(MutSpecAG$realm))
{
  if (i != 'Антарктика')
  {
    print(i)
    print(wilcox.test(MutSpecAG[MutSpecAG$realm == 'Антарктика' & MutSpecAG$Mut == 'C>T',]$MutSpec, 
                      MutSpecAG[MutSpecAG$realm == i & MutSpecAG$Mut == 'C>T',]$MutSpec))
  }
}

for (i in unique(MutSpecAG$Trophic_niche))
{
  if (i != 'Водные травоядные')
  {
    print(i)
    print(wilcox.test(MutSpecAG[MutSpecAG$Trophic_niche == 'Водные травоядные' & MutSpecAG$Mut == 'A>G',]$MutSpec, 
                      MutSpecAG[MutSpecAG$Trophic_niche == i & MutSpecAG$Mut == 'A>G',]$MutSpec))
  }
}

for (i in unique(MutSpecAG$Trophic_niche))
{
  if (i != 'Водные травоядные' & i != 'Падальщики')
  {
    print(i)
    print(wilcox.test(MutSpecAG[MutSpecAG$Trophic_niche == 'Водные травоядные' & MutSpecAG$Mut == 'C>T',]$MutSpec, 
                      MutSpecAG[MutSpecAG$Trophic_niche == i & MutSpecAG$Mut == 'C>T',]$MutSpec))
  }
}

wilcox.test(MutSpecAG[MutSpecAG$realm == 'Антарктика' & MutSpecAG$Mut == 'A>G',]$RawMutSpec, 
            MutSpecAG[MutSpecAG$realm != 'Антарктика' & MutSpecAG$Mut == 'A>G',]$RawMutSpec)

wilcox.test(MutSpecAG[MutSpecAG$realm == 'Антарктика' & MutSpecAG$Mut == 'C>T',]$RawMutSpec, 
            MutSpecAG[MutSpecAG$realm != 'Антарктика' & MutSpecAG$Mut == 'C>T',]$RawMutSpec)

wilcox.test(MutSpecAG[MutSpecAG$Trophic_niche == 'Водные травоядные' & MutSpecAG$Mut == 'A>G',]$RawMutSpec, 
            MutSpecAG[MutSpecAG$Trophic_niche != 'Водные травоядные' & MutSpecAG$Mut == 'A>G',]$RawMutSpec)

wilcox.test(MutSpecAG[MutSpecAG$Trophic_niche == 'Водные травоядные' & MutSpecAG$Mut == 'C>T',]$RawMutSpec, 
            MutSpecAG[MutSpecAG$Trophic_niche != 'Водные травоядные' & MutSpecAG$Mut == 'C>T',]$RawMutSpec)



pca_data = mut_data_syn1[,c(1,5,9)]
ex = reshape(data = pca_data, idvar = 'species_name',
             timevar = 'Mut',
             direction = 'wide')
names(ex) = c('species_name', 'T>G', 'T>C', 'T>A', 'C>G', 'C>A', 'C>T', 'G>A','G>T','G>C', 'A>T', 'A>C', 'A>G')
ex = merge(ex, ecozone_data, by = 'species_name')
ex = ex[,c(1,14,15,12,13,11,6,5,7,8,10,9,4,3,2)]
row.names(ex) = ex$species_name


stats_pca1 = prcomp(ex[,c(4,5,6,7,8,9,10,11,12,13,14,15)], center = TRUE, scale. = TRUE)
summary(stats_pca1)

bipl = ggbiplot(stats_pca1, labels.size = 2)
bipl


birds_ms_pca = data.frame(stats_pca1$x)
birds_ms_pca = birds_ms_pca[,c(1,2)]
birds_ms_pca$species_name = row.names(birds_ms_pca)
ex2 = merge(ex, birds_ms_pca, by = 'species_name')


g15 = ggplot(ex2, aes(x=PC1, color=russian_eco)) +
  geom_density(linewidth = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (38.2%)')+
  ylab('Плотность')+
  scale_colour_manual(name="Origin", values= c("black", "red", "black", "black", "black", "black", "black", "black", "black"))+
  theme(legend.position="none")


g15
g16 = ggplot(ex2, aes(x=PC2, color=russian_eco)) +
  geom_density(linewidth = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC2 (10.4%)')+
  ylab('Плотность')+
  scale_colour_manual(name="Origin", values= c("black", "red", "black", "black", "black", "black", "black", "black", "black"))+
  theme(axis.title.y=element_blank(), axis.ticks.y=element_blank())+
  guides(color = guide_legend(title = "Экозона"))


g16

g17 = ggplot(ex2, aes(x=PC1, color=russian_tn)) +
  geom_density(linewidth = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC1 (38.2%)')+
  ylab('Плотность')+
  scale_colour_manual(name="Origin", values= c("black", "red", "black", "black", "black", "black", "black", "black", "black", 'black'))+
  theme(legend.position="none")

g17
g18 = ggplot(ex2, aes(x=PC2, color=russian_tn)) +
  geom_density(linewidth = 1)+
  theme_bw()+
  theme(axis.title = element_text(size = 29), axis.text = element_text(size = 25))+
  xlab('PC2 (10.4%)')+
  ylab('Плотность')+
  scale_colour_manual(name="Origin", values= c("black", "red", "black", "black", "black", "black", "black", "black", "black", 'black'))+
  theme(axis.title.y=element_blank(), axis.ticks.y=element_blank())+
  guides(color = guide_legend(title = "Трофическая ниша"))

g18

f_den_stats2 = ggarrange(g15, g16, g17, g18,
                         ncol = 2, nrow = 2)
f_den_stats2


library(data.table)
pca_data = mut_data_syn1[,c(1,5,7,8,9)]

pca_data_shaped = dcast.data.table(setDT(pca_data), species_name+realm+Trophic_niche~Mut,
                                   value.var='MutSpec')
pca_data_shaped$TR_TS = (pca_data_shaped$`A>G`+pca_data_shaped$`C>T`+pca_data_shaped$`G>A`+pca_data_shaped$`T>C`)/(pca_data_shaped$`A>C`+pca_data_shaped$`A>T`+pca_data_shaped$`C>A`+pca_data_shaped$`C>G`+pca_data_shaped$`G>C`+pca_data_shaped$`G>T`+pca_data_shaped$`T>A`+pca_data_shaped$`T>G`)
pca_data_shaped$CT_GA = (pca_data_shaped$`C>T`/pca_data_shaped$`G>A`)
pca_data_shaped$AG_TC = (pca_data_shaped$`A>G`/pca_data_shaped$`T>C`)


realm_tr_ts = ggplot(data = pca_data_shaped, aes(x = realm, y = TR_TS))+
  geom_boxplot(fill = 'cyan', notch = TRUE)+
  ylim(0,75)+
  annotate("text", x=6, y=50, label= "Tr/Ts")+
  xlim(c('Антарктика', 'Неарктика', 'Палеарктика', 'Индо-Малайзия', 'Афротропики', 'Мадагаскар', 'Неотропики', 'Австралия', 'Океания'))+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
realm_tr_ts

niche_tr_ts = ggplot(data = pca_data_shaped, aes(x = Trophic_niche, y = TR_TS))+
  geom_boxplot(fill = 'orange', notch = TRUE)+
  ylim(0,75)+
  annotate("text", x=6, y=50, label= "Tr/Ts")+
  xlim(c('Водные травоядные', 'Падальщики', 'Позвоночные', 'Зерноядные', 'Наземные травоядные', 'Беспозвоночные', 'Водные хищники', 'Нектароядные', 'Всеядные', 'Плодоядные'))+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
niche_tr_ts

realm_tr_ts1 = ggplot(data = pca_data_shaped, aes(x = realm, y = CT_GA))+
  geom_boxplot(fill = 'green', notch = TRUE)+
  ylim(0,5)+
  annotate("text", x=6, y=3, label= "CT/GA")+
  xlim(c('Антарктика', 'Неарктика', 'Палеарктика', 'Индо-Малайзия', 'Афротропики', 'Мадагаскар', 'Неотропики', 'Австралия', 'Океания'))+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
niche_tr_ts1 = ggplot(data = pca_data_shaped, aes(x = Trophic_niche, y = CT_GA))+
  geom_boxplot(fill = 'yellow', notch = TRUE)+
  ylim(0,5)+
  annotate("text", x=2, y=1, label= "CT/GA")+
  xlim(c('Водные травоядные', 'Падальщики', 'Позвоночные', 'Зерноядные', 'Наземные травоядные', 'Беспозвоночные', 'Водные хищники', 'Нектароядные', 'Всеядные', 'Плодоядные'))+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(), axis.ticks.y=element_blank())
niche_tr_ts1

realm_tr_ts2 = ggplot(data = pca_data_shaped, aes(x = realm, y = AG_TC))+
  geom_boxplot(fill = 'blue', notch = TRUE)+
  ylim(0,5)+
  annotate("text", x=6, y=1, label= "AG/TC")+
  xlim(c('Антарктика', 'Неарктика', 'Палеарктика', 'Индо-Малайзия', 'Афротропики', 'Мадагаскар', 'Неотропики', 'Австралия', 'Океания'))+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank())
realm_tr_ts2

niche_tr_ts2 = ggplot(data = pca_data_shaped, aes(x = Trophic_niche, y = AG_TC))+
  geom_boxplot(fill = 'red', notch = TRUE)+
  ylim(0,5)+
  annotate("text", x=2, y=1, label= "AG/TC")+
  xlim(c('Водные травоядные', 'Падальщики', 'Позвоночные', 'Зерноядные', 'Наземные травоядные', 'Беспозвоночные', 'Водные хищники', 'Нектароядные', 'Всеядные', 'Плодоядные'))+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.ticks.y=element_blank())
niche_tr_ts2

tr_ts_graph= ggarrange(realm_tr_ts2, realm_tr_ts1,
                       ncol = 1, nrow = 2)

tr_ts_graph


tr_ts_graph1= ggarrange(niche_tr_ts2, niche_tr_ts1,
                       ncol = 1, nrow = 2)

tr_ts_graph1


#Valya's ideas
library(plotly)
ggplot(df_mtdna, aes(x = ghahSkew, y = Mass, color = gene_name))+
  geom_point()
ggplot(df_mtdna, aes(x = chthSkew, y = Mass, color = gene_name))+
  geom_point()
df_lowmass = df_mtdna[df_mtdna$Mass < 15000,]
ggplot(df_lowmass, aes(x = ghahSkew, y = Mass, color = gene_name))+
  geom_point()
ggplot(df_lowmass, aes(x = chthSkew, y = Mass, color = gene_name))+
  geom_point()
names_v = unique(df_mtdna$species_name)
df_short = data.frame()
for (i in names_v)
{
  df1 = df_mtdna[df_mtdna$species_name == i,]
  a = sum(df1$ghahSkew)/12
  b = sum(df1$chthSkew)/12
  v = sum(df1$Mass)/12
  ab = c(i, a, b, v)
  df_short = rbind(df_short, ab)
}
names(df_short) = c('species_name', 'GhAhSkew', 'ThChSkew', 'Mass')
df_short$Mass = as.numeric(df_short$Mass)
df_short$GhAhSkew = as.numeric(df_short$GhAhSkew)
df_short$ThChSkew = as.numeric(df_short$ThChSkew)

ggplot(df_short, aes(x = Mass, y = GhAhSkew))+
  geom_point()+
  geom_smooth(method = lm)


ggplot(df_short, aes(x = Mass, y = ThChSkew))+
  geom_point()+
  ylim(0.0,3.0)+
  xlim(0,3000)

df_short$log_mass = log10(df_short$Mass)
df_short$log_ghskew = log10(df_short$GhAhSkew)

ggplot(df_short, aes(x = GhAhSkew, y = log_mass))+
  geom_point()+
  geom_smooth(method=lm)

ggplot(df_short, aes(x = log_mass, y = GhAhSkew))+
  geom_point()+
  geom_smooth(method=lm)

ggplot(df_short, aes(x = log_mass, y = ThChSkew))+
  geom_point()+
  geom_smooth(method=lm)

lm1 = lm(GhAhSkew ~ Mass, data = df_short)
summary(lm1)

lm2 = lm(ThChSkew ~ Mass, data = df_short)
summary(lm2)


df_short_cytb = df_mtdna[df_mtdna$gene_name == 'CYTB',]
df_short_cytb$log_mass = log10(df_short_cytb$Mass)
ggplot(df_short_cytb, aes(x = log_mass, y = ghahSkew))+
  geom_point()+
  geom_smooth(method = lm)

ggplot(df_short_cytb, aes(x = log_mass, y = chthSkew))+
  geom_point()+
  geom_smooth(method = lm)

df_short_cox1 = df_mtdna[df_mtdna$gene_name == 'COX1',]
df_short_cox1$log_mass = log10(df_short_cox1$Mass)
ggplot(df_short_cox1, aes(x = log_mass, y = ghahSkew))+
  geom_point()+
  geom_smooth(method = lm)
ggplot(df_short_cox1, aes(x = log_mass, y = chthSkew))+
  geom_point()+
  geom_smooth(method = lm)

