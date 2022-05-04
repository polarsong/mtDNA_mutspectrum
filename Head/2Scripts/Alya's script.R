#changing data and PGLS analysis
#1 - antarctic 0 - afrotropic

rm(list=ls(all=T))

#Dima's analysis

library(dplyr)
library(ggplot2)

df = read.csv('../../Body/3Results/For_Bogdan.csv') #reading file
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

df_codons_realm = df_sgc %>% select(species_name, gene_name, realm,all_of(vec_all))

sp_sum_gen = data.frame(unique(df_codons_realm$species_name))

for ( codon in vec_all){
  
  sum_of_codon = aggregate(df_codons_realm[ ,codon], by = list(df_codons_realm$species_name), FUN = 'sum')[2]
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
df_eco = df_codons_realm[,c(1,3)]
codon_norm = merge(codon_norm, df_eco)

df_try = data.frame(unique(codon_norm))


final = data.frame()
for (org in 1:nrow(df_try)){
  sp_r = df_try[org,]
  
  vec_of_c = sp_r %>% select(TTC, TCC, TAC, TGC, CTC, CCC, CAC, CGC,
                             ATC, ACC, AAC, AGC, GTC, GCC, GAC, GGC)
  vec_of_a = sp_r %>% select(TTA, TCA, TGA, CTA, CCA, CAA, CGA,
                             ATA, ACA, AAA, GTA, GCA, GAA, GGA)
  
  med_c = median(as.numeric(vec_of_c), na.rm = TRUE)
  med_a = median(as.numeric(vec_of_a), na.rm = TRUE)
  sp_out = data.frame(sp_r$species_name, med_c, med_a, sp_r$realm) 
  final = rbind(final,sp_out)
}
names(final) = c('Species', 'med_c', 'med_a', 'realm')

df_antarctic = final[final$realm == 'Antarctic',]
df_every = final[final$realm != 'Antarctic', ]
df_antarctic$point = 1
df_every$point = 0
df_all = data.frame()
df_all = rbind(df_antarctic, df_every)
df_all$Species = gsub(' ','_', df_all$Species)
df_all$Species = as.character(df_all$Species)


#Alya's script
library(ape)
library(geiger)
library(caper)
tree <- read.tree("../../Body/1Raw/mtalign.aln.treefile.rooted")
row.names(df_all) = df_all$Species
tree_w = treedata(tree, df_all, sort=T, warnings=T)$phy

data<-as.data.frame(treedata(tree_w, df_all, sort=T, warnings=T)$data)

data$med_c = as.numeric(as.character(data$med_c))
data$med_a = as.numeric(as.character(data$med_a))
data$point= as.numeric(as.character(data$point))
data$realm= as.character(data$realm)
table(data$point)

MutComp = comparative.data(tree_w, data, Species, vcv=TRUE)

model = pgls(point ~ med_a + med_c, MutComp, lambda="ML")
summary(model)

#new birds data
avonet = read.csv2("../../Body/1Raw/Avonet_data.csv", header = TRUE, sep = ",")
names(avonet) = c('Species', 'Family', 'Order', 'Total_individuals', 'Female', 'Male', 'Unknown', 'Complete_measures',
                  'Beak_length_Culmen', 'Beak_length_Nares', 'Beak_width', 'Beak_depth', 'Tarsus_length', 'Wing_length',
                  'Kipps_distance', 'Secondary1', 'Hand_wing_index', 'Tail_length', 'Mass', 'Mass_source', 'Mass_refs_other',
                  'Inference', 'Treits_inferred', 'Reference_species', 'Habitat', 'Habitat_density', 'Migration', 'Trophic_level',
                  'Trophic_niche', 'Primary_lifestyle', 'Min_latitude', 'Max_latitude', 'Centroid_latitude', 'Centroid_longitude',
                  'Range_size', 'Species_status')
avonet$Species = gsub(' ','_', avonet$Species)
df_all1 = merge(df_all, avonet, by = 'Species')
#new PGLS
tree <- read.tree("../../Body/1Raw/mtalign.aln.treefile.rooted")
row.names(df_all1) = df_all1$Species
tree_w1 = treedata(tree, df_all1, sort=T, warnings=T)$phy
data1<-as.data.frame(treedata(tree_w1, df_all1, sort=T, warnings=T)$data)
data1$Beak_width = as.numeric(as.character(data1$Beak_width))
