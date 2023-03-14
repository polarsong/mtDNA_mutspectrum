rm(list=ls(all=TRUE))
library(ggplot2)
library(plotly)
library(ggpubr)
df_mtdna = read.csv('../../Head/2Scripts/Birds_dataset_paper.csv')
mut_data = read.table("C:/Users/User/Desktop/Birds mutspec results from Bogdan/mutspec12.tsv", header = TRUE, fill = TRUE)
#3 nucleotide against four-fold
df_mtdna$all_A_position_3 = df_mtdna$TTT + df_mtdna$TCT + df_mtdna$TAT + df_mtdna$TGT + df_mtdna$CTT + df_mtdna$CCT + df_mtdna$CAT + df_mtdna$CGT + df_mtdna$ATT + df_mtdna$ACT + df_mtdna$AAT + df_mtdna$AGT + df_mtdna$GTT + df_mtdna$GCT + df_mtdna$GAT + df_mtdna$GGT
df_mtdna$all_T_position_3 = df_mtdna$TTA + df_mtdna$TCA + df_mtdna$TAA + df_mtdna$TGA + df_mtdna$CTA + df_mtdna$CCA + df_mtdna$CAA + df_mtdna$CGA + df_mtdna$ATA + df_mtdna$ACA + df_mtdna$AAA + df_mtdna$AGA + df_mtdna$GTA + df_mtdna$GCA + df_mtdna$GAA + df_mtdna$GGA
df_mtdna$all_G_position_3 = df_mtdna$TTC + df_mtdna$TCC + df_mtdna$TAC + df_mtdna$TGC + df_mtdna$CTC + df_mtdna$CCC + df_mtdna$CAC + df_mtdna$CGC + df_mtdna$ATC + df_mtdna$ACC + df_mtdna$AAC + df_mtdna$AGC + df_mtdna$GTC + df_mtdna$GCC + df_mtdna$GAC + df_mtdna$GGC
df_mtdna$all_C_position_3 = df_mtdna$TTG + df_mtdna$TCG + df_mtdna$TAG + df_mtdna$TGG + df_mtdna$CTG + df_mtdna$CCG + df_mtdna$CAG + df_mtdna$CGG + df_mtdna$ATG + df_mtdna$ACG + df_mtdna$AAG + df_mtdna$AGG + df_mtdna$GTG + df_mtdna$GCG + df_mtdna$GAG + df_mtdna$GGG


df_mtdna$G_ratio_3_position = df_mtdna$neutral_c/df_mtdna$all_G_position_3
df_mtdna$C_ratio_3_position = df_mtdna$neutral_g/df_mtdna$all_C_position_3
df_mtdna$A_ratio_3_position = df_mtdna$neutral_T/df_mtdna$all_A_position_3
df_mtdna$T_ratio_3_position = df_mtdna$neutral_A/df_mtdna$all_T_position_3

f1 = ggplot(data = df_mtdna, aes(x = realm, y = T_ratio_3_position))+
  geom_boxplot(outlier.shape = NA)+
  xlab('realms')+
  ylab('Thymine ratio in third positions')
f1 = f1 +xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))
f1 = f1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
f2 = ggplot(data = df_mtdna, aes(x = realm, y = A_ratio_3_position))+
  geom_boxplot(outlier.shape = NA)+
  xlab('realms')+
  ylab('Adenine ratio in third positions')
f2 = f2 +xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))
f2 = f2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
f3 = ggplot(data = df_mtdna, aes(x = realm, y = G_ratio_3_position))+
  geom_boxplot(outlier.shape = NA)+
  xlab('realms')+
  ylab('Guanine ratio in third positions')
f3 = f3 +xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))
f3 = f3 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
f4 = ggplot(data = df_mtdna, aes(x = realm, y = C_ratio_3_position))+
  geom_boxplot(outlier.shape = NA)+
  xlab('realms')+
  ylab('Cytosine ratio in third positions')
f4 = f4 +xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))
f4 = f4 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
               
f5 = ggarrange(f1, f3, f4, f2, 
               ncol = 2, nrow = 2)
f5

df_look = df_mtdna[df_mtdna$realm == 'Antarctic',]
#Syn_mutations
mut_data_syn = mut_data[mut_data$Label == 'syn',]
mut_data_syn = mut_data_syn[,c(1,2,3,4,5,7,8)]
ecozone_data = df_mtdna[,c('species_name', 'realm', 'Trophic_niche')]
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

ggplot(data = mut_data_syn1, aes(x = Mut, y = MutSpec, fill = Trophic_niche))+
  geom_boxplot()

fig <- plot_ly(mut_data_syn1, x = ~Mut, y = ~MutSpec, color = ~Trophic_niche, type = "box")
fig <- fig %>% layout(boxmode = "group")

fig

#TR/TS
library(data.table)
pca_data = mut_data_syn1[,c(1,5,7,8,9)]

pca_data_shaped = dcast.data.table(setDT(pca_data), species_name+realm+Trophic_niche~Mut,
                                   value.var='MutSpec')
pca_data_shaped$TR_TS = (pca_data_shaped$`A>G`+pca_data_shaped$`C>T`+pca_data_shaped$`G>A`+pca_data_shaped$`T>C`)/(pca_data_shaped$`A>C`+pca_data_shaped$`A>T`+pca_data_shaped$`C>A`+pca_data_shaped$`C>G`+pca_data_shaped$`G>C`+pca_data_shaped$`G>T`+pca_data_shaped$`T>A`+pca_data_shaped$`T>G`)

realm_tr_ts = ggplot(data = pca_data_shaped, aes(x = realm, y = TR_TS))+
  geom_boxplot()+
  ylim(0,75)
realm_tr_ts = realm_tr_ts + xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))
realm_tr_ts = realm_tr_ts + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
realm_tr_ts

fig <- plot_ly(mut_data_syn1, x = ~realm, y = ~MutSpec, color = ~Trophic_niche, type = "box")
fig <- fig %>% layout(boxmode = "group")


fig
