rm(list=ls(all=T))
df_mtdna = read.csv('../../Body/3Results/Birds_mtDNA_data.csv')
df_mtdna$trophic_level = NA
df_mtdna$trophic_niche = NA
df_mtdna$foraging_niche = NA
df_mtdna = df_mtdna[, colSums(is.na(df_mtdna)) < nrow(df_mtdna)]
df_mtdna$ghahSkew = (df_mtdna$neutral_c - df_mtdna$neutral_T)/(df_mtdna$neutral_c + df_mtdna$neutral_T)
df_mtdna$fAn = df_mtdna$neutral_T/df_mtdna$neutral_amount
df_mtdna$fGn = df_mtdna$neutral_c/df_mtdna$neutral_amount
df_mtdna$fCn = df_mtdna$neutral_g/df_mtdna$neutral_amount
df_mtdna$fTn = df_mtdna$neutral_A/df_mtdna$neutral_amount
df_mtdna$Stg = df_mtdna$fTn + df_mtdna$fGn 
df_mtdna$Sac = df_mtdna$fAn + df_mtdna$fCn
df_mtdna$Stg_Sac = df_mtdna$Stg - df_mtdna$Sac
library(dplyr)
df_norm = df_mtdna[df_mtdna$gene_name != 'ND6',] #deleting ND6
vec_all = c('TTC','TTT','TCC','TCT','TAC','TAT','TGC','TGT',
            'TTA','TTG','TCA','TCG','TAA','TAG','TGA','TGG',
            'CTC','CTT','CCC','CCT','CAC','CAT','CGC','CGT',
            'CTA','CTG','CCA','CCG','CAA','CAG','CGA','CGG',
            'ATC','ATT','ACC','ACT','AAC','AAT','AGC','AGT',
            'ATA','ATG','ACA','ACG','AAA','AAG','AGA','AGG',
            'GTC','GTT','GCC','GCT','GAC','GAT','GGC','GGT',
            'GTA','GTG','GCA','GCG','GAA','GAG','GGA','GGG')

needed_codons = c('TTC','TCC','TAC','TGC',
                  'TTA','TCA','TAA','TGA',
                  'CTC','CCC','CAC','CGC',
                  'CTA','CCA','CAA','CGA',
                  'ATC','ACC','AAC','AGC',
                  'ATA','ACA','AAA','AGA',
                  'GTC','GCC','GAC','GGC',
                  'GTA','GCA','GAA','GGA')
dfn = c(names(df_norm))
sp_sum_gen = data.frame(unique(df_norm$species_name))

for ( codon in vec_all){
  
  sum_of_codon = aggregate(df_norm[ ,codon], by = list(df_norm$species_name), FUN = 'sum')[2]
  sp_sum_gen = cbind(sp_sum_gen, sum_of_codon)
  
}
names(sp_sum_gen) = c('Species', vec_all)

codon_norm = data.frame()

for (i in 1:nrow(sp_sum_gen)){
  org_gen = sp_sum_gen[i,]
  org_gen = as.vector(org_gen)
  df_out= data.frame(sp_sum_gen[i,]$Species) 
  for (codon in seq(from = 2, to = 65, by = 2)){
    if (org_gen[1,codon] == 0) {df_out = cbind(df_out, 0)}
    else {
      norm_cod = org_gen[1,codon] / (org_gen[1,codon+1] + org_gen[1,codon])
      df_out = cbind(df_out, norm_cod)
    }
  }
  names(df_out) = c('Species', needed_codons)
  codon_norm = rbind(codon_norm,df_out)
}
names(codon_norm) = c('species_name', needed_codons)
codon_norm = codon_norm %>% select(-c('TAA','AGA'))



final = data.frame()
for (org in 1:nrow(codon_norm)){
  sp_r = codon_norm[org,]
  
  vec_of_c = sp_r %>% select(TTC, TCC, TAC, TGC, CTC, CCC, CAC, CGC,
                             ATC, ACC, AAC, AGC, GTC, GCC, GAC, GGC)
  vec_of_a = sp_r %>% select(TTA, TCA, TGA, CTA, CCA, CAA, CGA,
                             ATA, ACA, AAA, GTA, GCA, GAA, GGA)
  
  med_c = median(as.numeric(vec_of_c), na.rm = TRUE)
  med_a = median(as.numeric(vec_of_a), na.rm = TRUE)
  sp_out = data.frame(sp_r$species_name, med_c, med_a) 
  final = rbind(final,sp_out)
}
names(final) = c('species_name', 'med_G', 'med_T')
df_mtdna = merge(df_mtdna, final, by = 'species_name')
df_eco = read.csv('../../Body/1Raw/Avonet_data.csv')
df_eco = df_eco[c(1,9,10,11,12,13,14,15,17,18,19,25,28,29,30)]
names(df_eco) = c('species_name','Beak_length_Culmen', 'Beak_length_Nares', 'Beak_width', 'Beak_depth', 'Tarsus_length', 'Wing_length','Kipps_distance', 'Hand_wing_index', 'Tail_length', 'Mass', 'Habitat', 'Trophic_level', 'Trophic_niche', 'Primary_lifestyle')
df_mtdna = merge(df_mtdna, df_eco, by = 'species_name')


library(ggbiplot)
library(ggplot2)
library(ggpubr)


f1 = ggplot(data = df_mtdna, aes(x = gene_name, y = fTn))+
  geom_boxplot()+
  xlab('Mitochondrial genes')+
  ylab('Thymine frequence')
f1 = f1 +xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))
f1 = f1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
f2 = ggplot(data = df_mtdna, aes(x = gene_name, y = fCn))+
  geom_boxplot()+
  xlab('Mitochondrial genes')+
  ylab('Cytosine frequence')
f2 = f2 +xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))
f2 = f2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
f3 = ggplot(data = df_mtdna, aes(x = gene_name, y = fGn))+
  geom_boxplot()+
  xlab('Mitochondrial genes')+
  ylab('Guanine frequence')
f3 = f3 +xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))
f3 = f3 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
f4 = ggplot(data = df_mtdna, aes(x = gene_name, y = fAn))+
  geom_boxplot()+
  xlab('Mitochondrial genes')+
  ylab('Adenine frequence')
f4 = f4 +xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))
f4 = f4 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
f5 = ggarrange(f1, f2, f3, f4, 
               ncol = 2, nrow = 2)
f5

cox1 = subset(df_mtdna, df_mtdna$gene_name == 'COX1',)
cox2 = subset(df_mtdna, df_mtdna$gene_name == 'COX2',)
cox3 = subset(df_mtdna, df_mtdna$gene_name == 'COX3',)
atp8 = subset(df_mtdna, df_mtdna$gene_name == 'ATP8',)
atp6 = subset(df_mtdna, df_mtdna$gene_name == 'ATP6',)
nd2 = subset(df_mtdna, df_mtdna$gene_name == 'ND2',)
nd3 = subset(df_mtdna, df_mtdna$gene_name == 'ND3',)
nd4l = subset(df_mtdna, df_mtdna$gene_name == 'ND4L',)
nd4 = subset(df_mtdna, df_mtdna$gene_name == 'ND4',)
nd5 = subset(df_mtdna, df_mtdna$gene_name == 'ND5',)
cytb = subset(df_mtdna, df_mtdna$gene_name == 'CYTB',)
nd6 = subset(df_mtdna, df_mtdna$gene_name == 'ND6',)
nd1 = subset(df_mtdna, df_mtdna$gene_name == 'ND1',)

skew_all = ggplot(data = df_mtdna, aes(x = gene_name, y = ghahSkew))+
  geom_boxplot()+
  xlab('Gene names')+
  ylab('GhAhSkew')
skew_all = skew_all + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))
skew_all

skew_eco = ggplot(data = df_mtdna, aes(x = realm, y = ghahSkew))+
  geom_boxplot()+
  xlab('Birds realms')+
  ylab('GhAhSkew')
skew_eco = skew_eco + xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))
skew_eco = skew_eco + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
skew_eco

stg_all = ggplot(data = df_mtdna, aes(x = gene_name, y = Stg_Sac))+
  geom_boxplot()+
  xlab('Gene names')+
  ylab('Stg-Sac')
stg_all = stg_all + xlim(c("COX1","COX2","ATP8","ATP6","COX3", "ND3", "ND4L","ND4","ND5","CYTB","ND6","ND1","ND2"))
stg_all

stg_eco = ggplot(data = df_mtdna, aes(x = realm, y = Stg_Sac))+
  geom_boxplot()+
  xlab('Birds realm')+
  ylab('Stg-Sac')
stg_eco = stg_eco + xlim(c('Antarctic', 'Nearctic', 'Palearctic', 'Indo_Malay', 'Afrotropic', 'Madagascar', 'Neotropic', 'Australian', 'Oceania'))
stg_eco = stg_eco + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
stg_eco

