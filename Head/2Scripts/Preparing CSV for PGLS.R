rm(list = ls(all=TRUE))
df_mtdna = read.csv('../../Body/3Results/Birds_mtDNA_data.csv')
df_mtdna$trophic_level = NA
df_mtdna$trophic_niche = NA
df_mtdna$foraging_niche = NA
df_mtdna = df_mtdna[, colSums(is.na(df_mtdna)) < nrow(df_mtdna)]

df_mtdna$ghahSkew = (df_mtdna$neutral_c - df_mtdna$neutral_T)/(df_mtdna$neutral_c + df_mtdna$neutral_T)
df_mtdna$fAn = df_mtdna$neutral_A/df_mtdna$neutral_amount
df_mtdna$fGn = df_mtdna$neutral_g/df_mtdna$neutral_amount
df_mtdna$fCn = df_mtdna$neutral_c/df_mtdna$neutral_amount
df_mtdna$fTn = df_mtdna$neutral_T/df_mtdna$neutral_amount
df_mtdna$Stg = df_mtdna$fAn + df_mtdna$fCn 
df_mtdna$Sac = df_mtdna$fTn + df_mtdna$fGn
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
names(final) = c('species_name', 'med_c', 'med_a')
df_mtdna = merge(df_mtdna, final, by = 'species_name')
df_eco = read.csv('../../Body/1Raw/Avonet_data.csv')
df_eco = df_eco[c(1,9,10,11,12,13,14,15,17,18,19,25,28,29,30)]
names(df_eco) = c('species_name','Beak_length_Culmen', 'Beak_length_Nares', 'Beak_width', 'Beak_depth', 'Tarsus_length', 'Wing_length','Kipps_distance', 'Hand_wing_index', 'Tail_length', 'Mass', 'Habitat', 'Trophic_level', 'Trophic_niche', 'Primary_lifestyle')
df_mtdna = merge(df_mtdna, df_eco, by = 'species_name')


unique(df_mtdna$Habitat)
df_mtdna$forest1 = 0
df_mtdna$grassland1 = 0
df_mtdna$human1 = 0
df_mtdna$wetland1 = 0
df_mtdna$shrubland1 = 0
df_mtdna$woodland1 = 0
df_mtdna$marine1 = 0
df_mtdna$riverline1 = 0
df_mtdna$coastal1 = 0
df_mtdna$rock1 = 0
df_mtdna$desert1 = 0

df1 = subset(df_mtdna, df_mtdna$Habitat == 'Forest',)
df2 = subset(df_mtdna, df_mtdna$Habitat == 'Grassland',)
df3 = subset(df_mtdna, df_mtdna$Habitat == 'Human Modified',)
df4 = subset(df_mtdna, df_mtdna$Habitat == 'Wetland',)
df5 = subset(df_mtdna, df_mtdna$Habitat == 'Shrubland',)
df6 = subset(df_mtdna, df_mtdna$Habitat == 'Woodland',)
df7 = subset(df_mtdna, df_mtdna$Habitat == 'Marine',)
df8 = subset(df_mtdna, df_mtdna$Habitat == 'Riverine',)
df9 = subset(df_mtdna, df_mtdna$Habitat == 'Coastal',)
df10 = subset(df_mtdna, df_mtdna$Habitat == 'Rock',)
df11 = subset(df_mtdna, df_mtdna$Habitat == 'Desert',)

df1$forest1 = 1
df2$grassland1 = 1
df3$human1 = 1
df4$wetland1 = 1
df5$shrubland1 = 1
df6$woodland1 = 1
df7$marine1 = 1
df8$riverline1 = 1
df9$coastal1 = 1
df10$rock1 = 1
df11$desert1 = 1
df_mtdna = rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11)

unique(df_mtdna$Trophic_level)

df_mtdna$carnivore1 = 0
df_mtdna$herbivore1 = 0
df_mtdna$omnivore1 = 0
df_mtdna$scavenger1 = 0

df1 = subset(df_mtdna, df_mtdna$Trophic_level == 'Carnivore',)
df2 = subset(df_mtdna, df_mtdna$Trophic_level == 'Herbivore',)
df3 = subset(df_mtdna, df_mtdna$Trophic_level == 'Omnivore',)
df4 = subset(df_mtdna, df_mtdna$Trophic_level == 'Scavenger',)

df1$carnivore1 = 1
df2$herbivore1 = 1
df3$omnivore1 = 1
df4$scavenger1 = 1

df_mtdna = rbind(df1, df2, df3, df4)

unique(df_mtdna$realm)

df_mtdna$australian1 = 0
df_mtdna$palearctic1 = 0
df_mtdna$indo_malay1 = 0
df_mtdna$madagascar1 = 0
df_mtdna$neotropic1 = 0
df_mtdna$nearctic1 = 0
df_mtdna$afrotropic1 = 0
df_mtdna$oceania1 = 0
df_mtdna$antarctic1 = 0

df1 = subset(df_mtdna, df_mtdna$realm == 'Australian',)
df2 = subset(df_mtdna, df_mtdna$realm == 'Palearctic',)
df3 = subset(df_mtdna, df_mtdna$realm == 'Indo_Malay',)
df4 = subset(df_mtdna, df_mtdna$realm == 'Madagascar',)
df5 = subset(df_mtdna, df_mtdna$realm == 'Neotropic',)
df6 = subset(df_mtdna, df_mtdna$realm == 'Nearctic',)
df7 = subset(df_mtdna, df_mtdna$realm == 'Afrotropic',)
df8 = subset(df_mtdna, df_mtdna$realm == 'Oceania',)
df9 = subset(df_mtdna, df_mtdna$realm == 'Antarctic',)

df1$australian1 = 1
df2$palearctic1 = 1
df3$indo_malay1 = 1
df4$madagascar1 = 1
df5$neotropic1 = 1
df6$nearctic1 = 1
df7$afrotropic1 = 1
df8$oceania1 = 1
df9$antarctic1 = 1

df_mtdna = rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9)

unique(df_mtdna$Trophic_niche)
df_mtdna$invertivore1 = 0
df_mtdna$vertivore1 = 0
df_mtdna$aquatic_predator1 = 0
df_mtdna$omnivoretn1 = 0
df_mtdna$granivore1 = 0
df_mtdna$frugivore1 = 0
df_mtdna$nectarivore1 = 0
df_mtdna$herbivore_aquatic1 = 0
df_mtdna$herbivore_terrestrial1 = 0
df_mtdna$scavengertn1 = 0

df1 = subset(df_mtdna, df_mtdna$Trophic_niche == 'Invertivore',)
df2 = subset(df_mtdna, df_mtdna$Trophic_niche == 'Vertivore',)
df3 = subset(df_mtdna, df_mtdna$Trophic_niche == 'Aquatic predator',)
df4 = subset(df_mtdna, df_mtdna$Trophic_niche == 'Omnivore',)
df5 = subset(df_mtdna, df_mtdna$Trophic_niche == 'Granivore',)
df6 = subset(df_mtdna, df_mtdna$Trophic_niche == 'Frugivore',)
df7 = subset(df_mtdna, df_mtdna$Trophic_niche == 'Nectarivore',)
df8 = subset(df_mtdna, df_mtdna$Trophic_niche == 'Herbivore aquatic',)
df9 = subset(df_mtdna, df_mtdna$Trophic_niche == 'Herbivore terrestrial',)
df10 = subset(df_mtdna, df_mtdna$Trophic_niche == 'Scavenger',)

df1$invertivore1 = 1
df2$vertivore1 = 1
df3$aquatic_predator1 = 1
df4$omnivoretn1 = 1
df5$granivore1 = 1
df6$frugivore1 = 1
df7$nectarivore1 = 1
df8$herbivore_aquatic1 = 1
df9$herbivore_terrestrial1 = 1
df10$scavengertn1 = 1

df_mtdna = rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10)

unique(df_mtdna$Primary_lifestyle)
df_mtdna$insessorial1 = 0
df_mtdna$terrestrial1 = 0
df_mtdna$generalist1 = 0
df_mtdna$aerial1 = 0
df_mtdna$aquatic1 = 0


df1 = subset(df_mtdna, df_mtdna$Primary_lifestyle == 'Insessorial',)
df2 = subset(df_mtdna, df_mtdna$Primary_lifestyle == 'Terrestrial',)
df3 = subset(df_mtdna, df_mtdna$Primary_lifestyle == 'Generalist',)
df4 = subset(df_mtdna, df_mtdna$Primary_lifestyle == 'Aerial',)
df5 = subset(df_mtdna, df_mtdna$Primary_lifestyle == 'Aquatic',)

df1$insessorial1 = 1
df2$terrestrial1 = 1
df3$generalist1 = 1
df4$aerial1 = 1
df5$aquatic1 = 1

df_mtdna = rbind(df1, df2, df3, df4, df5)

write.csv(df_mtdna, file = 'Table_for_PGLS')

df_pgls = df_mtdna[c(1,2,5,26,27,28,29,30)]
df_pgls = df_pgls[df_pgls$gene_name != 'ND6',]
vecn = c('neutral_A', 'neutral_g', 'neutral_c', 'neutral_T', 'neutral_amount')
sp_sum_nucl = data.frame(unique(df_pgls$species_name))
for (i in vecn){
  
  sum_of_nucl = aggregate(df_pgls[ ,i], by = list(df_pgls$species_name), FUN = 'sum')[2]
  sp_sum_nucl = cbind(sp_sum_nucl, sum_of_nucl)
  
}
names(sp_sum_nucl) = c('species_name', vecn)
sp_sum_nucl$fTn = (sp_sum_nucl$neutral_A)/(sp_sum_nucl$neutral_amount)
sp_sum_nucl$fCn = (sp_sum_nucl$neutral_g)/(sp_sum_nucl$neutral_amount)
sp_sum_nucl$fAn = (sp_sum_nucl$neutral_T)/(sp_sum_nucl$neutral_amount)
sp_sum_nucl$fGn = (sp_sum_nucl$neutral_c)/(sp_sum_nucl$neutral_amount)
sp_sum_nucl$ghahSkew = (sp_sum_nucl$neutral_c - sp_sum_nucl$neutral_T)/(sp_sum_nucl$neutral_c + sp_sum_nucl$neutral_T)
sp_sum_nucl$Stg_Sac = (sp_sum_nucl$fTn+sp_sum_nucl$fGn)-(sp_sum_nucl$fAn+sp_sum_nucl$fCn)

df_pgls1 = unique(df_mtdna[c(1,2,104,105,106,107,108,109,110,111,112,113,114,115,120:158)])
df_finale = merge(df_pgls1, sp_sum_nucl)
df_finale$realm = NA
df_finale$neutral_A = NA
df_finale$neutral_c = NA
df_finale$neutral_g = NA
df_finale$neutral_T = NA
df_finale$neutral_amount = NA
df_finale = df_finale[, colSums(is.na(df_finale)) < nrow(df_finale)]

df_finale$species_name = gsub(" ", "_", df_finale$species_name)
row.names(df_finale) = df_finale$species_name

write.csv(df_finale, file = 'Table_for_PGLS.csv', row.names = FALSE)

a = read.csv("../../Head/2Scripts/Table_for_PGLS.csv")


