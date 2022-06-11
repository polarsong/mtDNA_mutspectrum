#final math
rm(list = ls(all.names = TRUE))

library(dplyr)
library(ggplot2)

df = read.csv('../../Body/3Results/For_Bogdan.csv')
df = df[df$gene_name != 'ND6',] #deleting ND6
df_sgc = df[,c(1,2,3,4,5,8, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53,
               54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77,
               78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97)] #getting codon usage


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

df_codons_realm = df_sgc %>% select(species_name, gene_name, realm, trophic_niche, trophic_level, all_of(vec_all))

sp_sum_gen = data.frame(unique(df_codons_realm$species_name))

for ( codon in vec_all)
{
  
  sum_of_codon = aggregate(df_codons_realm[ ,codon], by = list(df_codons_realm$species_name), FUN = 'sum')[2]
  sp_sum_gen = cbind(sp_sum_gen, sum_of_codon)
  
}
names(sp_sum_gen) = c('Species', vec_all)


codon_norm = data.frame()

for (i in 1:nrow(sp_sum_gen))
{
  org_gen = sp_sum_gen[i,]
  org_gen = as.vector(org_gen)
  df_out= data.frame(sp_sum_gen[i,]$Species) 
  for (codon in seq(from = 2, to = 65, by = 2))
  {
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
df_eco = df_codons_realm[,c(1,3,4,5)]
codon_norm = merge(codon_norm, df_eco)

df_try = data.frame(unique(codon_norm))


final = data.frame()
for (org in 1:nrow(df_try))
{
  sp_r = df_try[org,]
  
  vec_of_c = sp_r %>% select(TTC, TCC, TAC, TGC, CTC, CCC, CAC, CGC,
                             ATC, ACC, AAC, AGC, GTC, GCC, GAC, GGC)
  vec_of_a = sp_r %>% select(TTA, TCA, TGA, CTA, CCA, CAA, CGA,
                             ATA, ACA, AAA, GTA, GCA, GAA, GGA)
  
  med_c = median(as.numeric(vec_of_c), na.rm = TRUE)
  med_a = median(as.numeric(vec_of_a), na.rm = TRUE)
  sp_out = data.frame(sp_r$species_name, med_c, med_a, sp_r$realm, sp_r$trophic_level, sp_r$trophic_niche) 
  final = rbind(final,sp_out)
}
names(final) = c('species_name', 'med_c', 'med_a', 'realm', 'trophic_level', 'trophic_niche')


df_antarctic = final[final$realm == 'Antarctic',]

df_afrotropic = final[final$realm == 'Afrotropic',]

sample(df_afrotropic$med_c, 8)


cor.test(df_antarctic$med_c,sample(df_afrotropic$med_c, 8), method = 'spearman')

#course work
g1 = ggplot(data = final, aes(x = realm, y = med_a))+
  geom_violin()+
  xlab("Ёкозона")+
  ylab('ћедиана јденина')
g1
g2 = ggplot(data = final, aes(x = realm, y = med_c))+
  geom_violin()+
  xlab("Ёкозона")+
  ylab('ћедиана ÷итозина')
g2

g3 = ggplot(data = final, aes(x = realm, y = med_a))+
  geom_boxplot()+
  xlab("Ёкозона")+
  ylab('ћедиана јденина')
g3
g4 = ggplot(data = final, aes(x = realm, y = med_c))+
  geom_boxplot()+
  xlab("Ёкозона")+
  ylab('ћедиана ÷итозина')
g4

g5 = ggplot(data = final, aes(x = trophic_level, y = med_a))+
  geom_violin()+
  xlab("“рофический уровень")+
  ylab('ћедиана јденина')
g5
#work here
g6 = ggplot(data = final, aes(x = trophic_level, y = med_c))+
  geom_violin()+
  xlab("“рофический уровень")+
  ylab('ћедиана ÷итозина')
g6

g7 = ggplot(data = final, aes(x = trophic_level, y = med_a))+
  geom_boxplot()+
  xlab("“рофический уровень")+
  ylab('ћедиана јденина')
g7
g8 = ggplot(data = final, aes(x = trophic_level, y = med_c))+
  geom_boxplot()+
  xlab("“рофический уровень")+
  ylab('ћедиана ÷итозина')
g8



#Alya's PGLS
library(ape)
library(geiger)
library(caper)

df = read.csv('../../Body/3Results/For_Bogdan.csv')


tree <- read.tree('../../Body/3Results/phylo.treefile')
row.names(MamGt) = MamGt$Species
tree_w = treedata(tree, df[, c('species_name', )],
                  sort=T, warnings=T)$phy
data<-as.data.frame(treedata(tree_w, MamGt[, c(СSpeciesТ, СTsTvТ, СT_CТ, СTC_TCGAТ,СG_AТ,СTC_TATGTCТ,СGenerationLength_dТ,СNumOfFourFoldMutInCytBТ)],
                             sort=T, warnings=T)$data)
nrow(data)
data$Species = as.character(data$Species)
data$TsTv = as.numeric(as.character(data$TsTv))
data$T_C = as.numeric(as.character(data$T_C))
data$G_A = as.numeric(as.character(data$G_A))
data$TC_TCGA = as.numeric(as.character(data$TC_TCGA))
data$TC_TATGTC = as.numeric(as.character(data$TC_TATGTC))
data$GenerationLength_d = as.numeric(as.character(data$GenerationLength_d))
data$NumOfFourFoldMutInCytB = as.numeric(as.character(data$NumOfFourFoldMutInCytB))
MutComp = comparative.data(tree_w, data, Species, vcv=TRUE)
summary(pgls(TsTv ~ log2(GenerationLength_d), MutComp, lambda=УMLФ)) # 1A
summary(pgls(TsTv ~ log2(GenerationLength_d) + log2(NumOfFourFoldMutInCytB), MutComp, lambda=УMLФ)) # 1B
summary(pgls(TsTv ~ 0 + log2(GenerationLength_d) + log2(NumOfFourFoldMutInCytB), MutComp, lambda=УMLФ)) # 1C
summary(pgls(T_C ~ log2(GenerationLength_d), MutComp, lambda=УMLФ)) # 2A
summary(pgls(T_C ~ 0 + log2(GenerationLength_d), MutComp, lambda=УMLФ)) # 2B
summary(pgls(T_C ~ 0 + log2(GenerationLength_d)  + log2(NumOfFourFoldMutInCytB), MutComp, lambda=УMLФ)) # 2C
summary(pgls(TC_TCGA ~ log2(GenerationLength_d), MutComp, lambda=УMLФ)) # 3A
summary(pgls(TC_TCGA ~ 0 + log2(GenerationLength_d), MutComp, lambda=УMLФ)) # 3B
summary(pgls(TC_TCGA ~ 0 + log2(GenerationLength_d) + log2(NumOfFourFoldMutInCytB), MutComp, lambda=УMLФ)) # 3C
summary(pgls(TC_TATGTC ~ log2(GenerationLength_d), MutComp, lambda=УMLФ)) # 4A
summary(pgls(TC_TATGTC ~ log2(GenerationLength_d) + log2(NumOfFourFoldMutInCytB), MutComp, lambda=УMLФ)) # 4B

