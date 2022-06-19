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
#End of analysis

#Making points with realms
df_antarctic = final[final$realm == 'Antarctic' | final$realm == 'Nearctic' | final$realm == 'Palearctic',]
df_every = final[final$realm != 'Antarctic' & final$realm != 'Nearctic' & final$realm != 'Palearctic',]
df_antarctic$point = 1
df_every$point = 0
df_all = data.frame()
df_all = rbind(df_antarctic, df_every)
df_all$Species = gsub(' ','_', df_all$Species)
df_all$Species = as.character(df_all$Species)
#end




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
#end of Alya's script

#Bogdan tree
tree1 = read.tree("../../Body/3Results/phylo.treefile")
tree_w3 = treedata(tree1, df_all, sort=T, warnings=T)$phy
data3 = as.data.frame(treedata(tree_w3, df_all, sort=T, warnings=T)$data)
MutComp3 = comparative.data(tree_w3, data3, Species, vcv = TRUE)
#end


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
#end


#all ecology PGLS analysis
tree <- read.tree("../../Body/1Raw/mtalign.aln.treefile.rooted")
phy=multi2di(tree1)
row.names(df_all1) = df_all1$Species
tree_all = treedata(phy, df_all1, sort=T, warnings=T)$phy
data_all<-as.data.frame(treedata(tree_all, df_all1, sort=T, warnings=T)$data)



data_all[is.na(data_all)] <- 0
data_all$Beak_width = as.numeric(as.character(data_all$Beak_width))
data_all$med_c = as.numeric(as.character(data_all$med_c))
data_all$med_a = as.numeric(as.character(data_all$med_a))
data_all$point= as.numeric(as.character(data_all$point))
data_all$realm= as.character(data_all$realm)
data_all$Tail_length = as.numeric(as.character(data_all$Tail_length))
MutComp_all = comparative.data(tree_all, data_all, Species, vcv = TRUE)
model_all = pgls(med_c ~ Beak_width + Tail_length, MutComp_all, lambda = "ML")

summary(model_all)

#med_c against med_a
ggplot(df_all1, aes(x = med_c, y = med_a))+
  geom_point()+
  geom_smooth(method = lm)+
  geom_text(label=rownames(df_all1))
plot(df_all1$med_c, df_all1$med_a)
ggplot(df_all1, aes(x = med_c, y = med_a))+
  geom_point()




#new work
tree1 = read.tree("../../Body/3Results/phylo.treefile")
phy=multi2di(tree1)
row.names(df_all1) = df_all1$Species
tree_all = treedata(phy, df_all1, sort=T, warnings=T)$phy
data_all<-as.data.frame(treedata(tree_all, df_all1, sort=T, warnings=T)$data)

data_all = data_all[,c(1,2,3,4,6,7,13,14,15,16,17,18,19,21,22,23,29,32,33,34,40)]
data_all$Beak_width = as.numeric(as.character(data_all$Beak_width))
data_all$med_c = as.numeric(as.character(data_all$med_c))
data_all$med_a = as.numeric(as.character(data_all$med_a))
data_all$realm= as.character(data_all$realm)
data_all$Tail_length = as.numeric(as.character(data_all$Tail_length))
data_all$Family = as.character(data_all$realm)
data_all$Order = as.character(data_all$Order)
data_all$Beak_length_Culmen = as.numeric(as.character(data_all$Beak_length_Culmen))
data_all$Beak_length_Nares = as.numeric(as.character(data_all$Beak_length_Nares))
data_all$Beak_depth = as.numeric(as.character(data_all$Beak_depth))
data_all$Tarsus_length = as.numeric(as.character(data_all$Tarsus_length))
data_all$Wing_length = as.numeric(as.character(data_all$Wing_length))
data_all$Kipps_distance = as.numeric(as.character(data_all$Kipps_distance))
data_all$Hand_wing_index = as.numeric(as.character(data_all$Hand_wing_index))
data_all$Mass = as.numeric(as.character(data_all$Mass))
data_all$Habitat = as.character(data_all$Habitat)
data_all$Trophic_level = as.character(data_all$Trophic_level)
data_all$Trophic_niche = as.character(data_all$Trophic_niche)
data_all$Primary_lifestyle = as.character(data_all$Primary_lifestyle)
data_all$Species_status = as.character(data_all$Species_status)

MutComp_all = comparative.data(tree_all, data_all, Species, vcv = TRUE)
model_all = pgls(med_a ~ Beak_width, MutComp_all, lambda = "ML")
summary(model_all)
summary(model_all)$r.squared #to get r squared
a[2]
#end



#PGLS table 
eco_vec = c('Beak_length_Culmen', 'Beak_length_Nares', 'Beak_width', 'Beak_depth', 'Tarsus_length',
            'Wing_length', 'Kipps_distance', 'Hand_wing_index', 'Tail_length', 'Mass')
gen_vec = c('med_a', 'med_c')
megapgls = data.frame()
for (i in gen_vec)
{
  for (j in eco_vec)
  {
    model_meda = pgls(data_all[,i] ~ data_all[,j], MutComp_all, lambda = 'ML')
    temp = summary(model_meda)$coefficients[,4]
    megapgls = rbind(megapgls, c(j, temp[2]))
  }
  names(megapgls) = c("ecology", "p-value")
  write.table(megapgls,
              file = 'all_pgls.txt',
              append = F,
              sep = ",",
              row.names = F,
              col.names = T,
              na="",
              quote = F)
  write( "\n Next Dataframe \n", file = "all_pgls.txt", append = T)
}

effect_size() - ?

#names(medapgls) = c("ecology", "p-value")
#write.csv(medapgls, "med_a_pgls.csv")



#one file 
write.table( df,
             file = "dataframes.txt", 
             append = F,
             sep = ",",
             row.names = F,
             col.names = T,
             na="",
             quote = F)
write( "\n Next Dataframe \n", file = "dataframes.txt", append = T)
write.table( df1,
             file = "dataframes.txt", 
             append = T,
             sep = ",",
             row.names = F,
             col.names = T,
             na="",
             quote = F)


