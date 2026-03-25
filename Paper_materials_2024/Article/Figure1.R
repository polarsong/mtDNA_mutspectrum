rm(list = ls(all=TRUE))
library(ggplot2)
library(ggpubr)
library(dplyr)
df_mtdna = read.csv('../Work_with_Andrey/Birds_dataset_paper.csv', header = TRUE, sep = ';')
df_nd6 = read.csv('../Birds_mtDNA_data.csv')


df_nd6$GhAhSkew = (df_nd6$neutral_c- df_nd6$neutral_T)/(df_nd6$neutral_c + df_nd6$neutral_T)
df_nd6$ThChSkew = (df_nd6$neutral_A - df_nd6$neutral_g)/(df_nd6$neutral_A + df_nd6$neutral_g)
df_nd6$fTn = df_nd6$neutral_A/df_nd6$neutral_amount
df_nd6$fAn = df_nd6$neutral_T/df_nd6$neutral_amount
df_nd6$fCn = df_nd6$neutral_g/df_nd6$neutral_amount
df_nd6$fGn = df_nd6$neutral_c/df_nd6$neutral_amount

SynNuc = read.table('AllGenesCodonUsageNoOverlap.txt', header = TRUE, sep = '\t')
SynNuc$ghahSkew = ((SynNuc$NeutralC - SynNuc$NeutralT))/((SynNuc$NeutralC + SynNuc$NeutralT))
SynNuc$chthSkew = ((SynNuc$NeutralA - SynNuc$NeutralG))/((SynNuc$NeutralA + SynNuc$NeutralG))
new_mam = SynNuc[, c('Species', 'Gene', "CCA",
                     "CCC", "CCG", "CCT", "TTA", "TTC", "TTG", "TTT", 'AAA', 'GAA', 'TAA', 'CAA', 'AGG', 'GGG', 'TGG', 'CGG', 'ghahSkew','chthSkew')]
new_mam$Сlass = 'Mammalia'
new_bird = df_nd6[, c('species_name', 'gene_name', "CCA",
                      "CCC", "CCG", "CCT", "TTA", "TTC", "TTG", "TTT", 'AAA', 'GAA', 'TAA', 'CAA', 'AGG', 'GGG', 'TGG', 'CGG', 'GhAhSkew','ThChSkew')]
new_bird$Сlass = 'Aves'
new_bird$species_name = gsub(' ', '_', new_bird$species_name)
new_mam$Gene[new_mam$Gene == 'CytB'] = 'CYTB'
names(new_mam) = c('species_name', 'gene_name', "CCA","CCC", "CCG", "CCT", "TTA", "TTC", "TTG", "TTT", 'AAA', 'GAA', 'TAA', 'CAA', 'AGG', 'GGG', 'TGG', 'CGG','GhAhSkew', 'ThChSkew', 'Class')
names(new_bird) = c('species_name', 'gene_name', "CCA","CCC", "CCG", "CCT", "TTA", "TTC", "TTG", "TTT", 'AAA', 'GAA', 'TAA', 'CAA', 'AGG', 'GGG', 'TGG', 'CGG','GhAhSkew', 'ThChSkew', 'Class')

new_big = rbind(new_mam, new_bird)
graph1 = ggplot(new_big, aes(x = gene_name, y = GhAhSkew, fill = Class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Mitochondrial genes')+
  ylab('GhAhSkew')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5",'CYTB',"ND6","ND1","ND2"))+
  ylim(-1,1)+
  annotate('text', x = 4.5, y = -0.75, label = 'N birds = 766')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")

graph2 = ggplot(new_big, aes(x = gene_name, y = ThChSkew, fill = Class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlab('Mitochondrial genes')+
  ylab('ThChSkew')+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5",'CYTB',"ND6","ND1","ND2"))+
  ylim(-1,1)+
  annotate('text', x = 4.5, y = -0.75, label = 'N mammals = 4356')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")

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
df_cytb$color = 'blue'
df_cytb[df_cytb$Mut =='C>G' | df_cytb$Mut == 'G>C',]$color = 'black'
df_cytb[df_cytb$Mut =='C>T' | df_cytb$Mut == 'G>A',]$color = 'red'
df_cytb[df_cytb$Mut =='T>A' | df_cytb$Mut == 'A>T',]$color = 'grey'
df_cytb[df_cytb$Mut =='T>C' | df_cytb$Mut == 'A>G',]$color = 'green'
df_cytb[df_cytb$Mut =='T>G' | df_cytb$Mut == 'A>C',]$color = 'pink'
graph3 = ggplot(df_cytb, aes(x = Mut, y = MutSpec))+
  geom_bar(stat = 'identity')+
  xlim(c("C>A","G>T","C>G","G>C","C>T", "G>A", "T>A","A>T","T>C",'A>G',"T>G","A>C"))+
  ylab('Mutspec for CytB')
graph3
graph_final = ggarrange(graph3, graph3, graph1, graph2,
                     ncol = 2, nrow = 2)
graph_final

#lm
new_big$classbinar = 0
new_big[new_big$Class == 'Aves',]$classbinar = 1
new_big$TBSS = 1
new_big[new_big$gene_name == 'COX2',]$TBSS = 2
new_big[new_big$gene_name == 'ATP8',]$TBSS = 3
new_big[new_big$gene_name == 'ATP6',]$TBSS = 4
new_big[new_big$gene_name == 'COX3',]$TBSS = 5
new_big[new_big$gene_name == 'ND3',]$TBSS = 6
new_big[new_big$gene_name == 'ND4L',]$TBSS = 7
new_big[new_big$gene_name == 'ND4',]$TBSS = 8
new_big[new_big$gene_name == 'ND5',]$TBSS = 9
new_big[new_big$gene_name == 'CYTB',]$TBSS = 10
new_big[new_big$gene_name == 'ND6',]$TBSS = 11
new_big[new_big$gene_name == 'ND1',]$TBSS = 12
new_big[new_big$gene_name == 'ND2',]$TBSS = 13
new_big$GhAhSkew_abs = new_big$GhAhSkew
new_big$ThChSkew_abs = new_big$ThChSkew
absgh = new_big[new_big$gene_name == 'ND6',]
absgh$GhAhSkew_abs = abs(absgh$GhAhSkew_abs)
absgh$ThChSkew_abs = abs(absgh$ThChSkew_abs)
new_big_cut = new_big[new_big$gene_name != 'ND6',]
new_big_lm = rbind(absgh, new_big_cut)
TBSS_lmGh = lm(GhAhSkew_abs ~ scale(classbinar) + scale(TBSS), data = new_big_lm)
summary(TBSS_lmGh)
summary(TBSS_lmGh)$coefficient
TBSS_lmGh_1 = lm(GhAhSkew_abs ~ scale(classbinar) * scale(TBSS), data = new_big_lm)
summary(TBSS_lmGh_1)


TBSS_lmTh = lm(ThChSkew_abs ~ scale(classbinar) + scale(TBSS), data = new_big_lm)
summary(TBSS_lmTh)
TBSS_lmTh_1 = lm(ThChSkew_abs ~ scale(classbinar) * scale(TBSS), data = new_big_lm)
summary(TBSS_lmTh_1)


sigma(TBSS_lm)/mean(new_big_lm$GhAhSkew_abs)

#AA shift
new_big$Pro_Phe = (new_big$CCT+new_big$CCA+new_big$CCG+new_big$CCC)/(new_big$TTT+new_big$TTC)
new_big$Pro_PheLeu = (new_big$CCT+new_big$CCA+new_big$CCG+new_big$CCC)/(new_big$TTT+new_big$TTC+new_big$TTA+new_big$TTG)
new_big12g = new_big[new_big$gene_name != 'ND6',]
ggplot(new_big, aes(x = Class, y = Pro_Phe))+
  geom_boxplot(outlier.alpha = FALSE)+
  ylim(0,3)
ggplot(new_big, aes(x = Class, y = Pro_PheLeu))+
  geom_boxplot(outlier.alpha = FALSE)+
  ylim(0,1.8)
ggplot(new_big12g, aes(x = Class, y = Pro_Phe))+
  geom_boxplot(outlier.alpha = FALSE)+
  ylim(0,3)
graph_3 = ggplot(new_big12g, aes(x = Class, y = Pro_PheLeu))+
  geom_boxplot(outlier.alpha = FALSE, notch = TRUE)+
  ylim(0,1.8)
#stats
wilcox.test(new_big12g[new_big12g$Class == 'Mammalia',]$Pro_PheLeu,new_big12g[new_big12g$Class == 'Aves',]$Pro_PheLeu)

#ND6
new_nd6 = new_big[new_big$gene_name == 'ND6',]
new_nd6$nd6ppl = (new_nd6$AGG+new_nd6$GGG+new_nd6$TGG+new_nd6$CGG)/(new_nd6$AAA+new_nd6$GAA+new_nd6$TAA+new_nd6$CAA)
new_nd6$nd6ppl_1 = (new_nd6$AAA+new_nd6$GAA+new_nd6$TAA+new_nd6$CAA)/(new_nd6$AGG+new_nd6$GGG+new_nd6$TGG+new_nd6$CGG)
ggplot(new_nd6, aes(x = Class, y = nd6ppl))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)
ggplot(new_nd6, aes(x = Class, y = nd6ppl_1))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  ylim(0,0.45)
wilcox.test(new_nd6[new_nd6$Class == 'Aves',]$nd6ppl_1, new_nd6[new_nd6$Class == 'Mammalia',]$nd6ppl_1)

#propheleu
df1 = new_big12g[,c('species_name', 'gene_name', 'Class', 'Pro_PheLeu')]
df2 = new_nd6[,c('species_name','gene_name', 'Class', 'nd6ppl_1')]
names(df2) = c('species_name', 'gene_name','Class', 'Pro_PheLeu')
df_3 = rbind(df1, df2)
ggplot(df_3, aes(x = gene_name, y = Pro_PheLeu, fill = Class))+
  geom_boxplot(notch = TRUE, outlier.alpha = FALSE)+
  xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5",'CYTB',"ND6","ND1","ND2"))+
  ylim(0,4.1)
new_big_exp = new_big[new_big$Pro_PheLeu != 'Inf',]
TBSS_ppl = lm(Pro_PheLeu ~ scale(classbinar) + scale(TBSS), data = new_big_exp)
summary(TBSS_ppl)


#mutspec 
mutspec = read.csv("MutSups.csv")
mutspec = mutspec[,c(1,4,6)]
names(mutspec) = mutspec[c(1),]
mutspec = mutspec[c(2:193),]
mutspec$Mut = substr(mutspec$Mut, 3,5)
mutspec$Aves = as.numeric(as.character(mutspec$Aves))
mutspec$Mammalia = as.numeric(as.character(mutspec$Mammalia))
ggplot(mutspec, aes(x = Mut, y = Aves))+
  geom_boxplot()
ggplot(mutspec, aes(x = Mut, y = Mammalia))+
  geom_boxplot()
graph_4 = ggplot(mutspec, aes(x = Mut, y = Aves))+
  geom_bar(stat = 'identity')
graph_5 = ggplot(mutspec, aes(x = Mut, y = Mammalia))+
  geom_bar(stat = 'identity')
graph_final = ggarrange(graph_4, graph_5, graph1, graph2,
                        ncol = 2, nrow = 2)
