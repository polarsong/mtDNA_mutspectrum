rm(list = ls(all=TRUE))
library(ggplot2)
df_midori = read.csv('Midori2_new_mutspec.csv')
df_midori = df_midori[,c(2,5,6,7)]
df_midori$gene_and_species = paste(df_midori$Gene, df_midori$Species)
df_midori$gene_and_species = gsub(' ', '_', df_midori$gene_and_species)
df_midori = df_midori[,c(1,2,5)]
df_mtdna = read.csv('../../Paper_materials_2024/Birds_dataset_paper.csv')
df_midori1 = reshape(data = df_midori, idvar = 'gene_and_species',
                     timevar = 'Mut',
                     direction = 'wide')
library(stringr)
df_midori1[,c('gene', 'species_name')] = str_split_fixed(df_midori1$gene_and_species, "_", 2)
df_midori1 = df_midori1[,c(2:15)]
names(df_midori1) = c('Mutation_TC_midori', 'Mutation_GA_midori', 'Mutation_CT_midori', 'Mutation_AG_midori',
                      'Mutation_AT_midori', 'Mutation_GT_midori', 'Mutation_TA_midori', 'Mutation_CG_midori',
                      'Mutation_TG_midori', 'Mutation_CA_midori', 'Mutation_GC_midori', 'Mutation_AC_midori',
                      'gene', 'species_name')
nd1 = df_midori1[df_midori1$gene == "ND1",]
nd2 = df_midori1[df_midori1$gene == "ND2",]
nd6 = df_midori1[df_midori1$gene == "ND6",]
cox1 = df_midori1[df_midori1$gene == "CO1",]
cytb = df_midori1[df_midori1$gene == "Cytb",]
atp8 = df_midori1[df_midori1$gene == "A8",]
nd3 = df_midori1[df_midori1$gene == "ND3",]
cox3 = df_midori1[df_midori1$gene == "CO3",]
cox2 = df_midori1[df_midori1$gene == "CO2",]
nd4 = df_midori1[df_midori1$gene == "ND4",]
atp6 = df_midori1[df_midori1$gene == "A6",]
nd5 = df_midori1[df_midori1$gene == "ND5",]
nd4l = df_midori1[df_midori1$gene == "ND4L",]

df_need = data.frame()
for (i in unique(df_mtdna$species_name))
{
  a = df_mtdna[df_mtdna$species_name == i,]
  b = sum(a$ghahSkew)/12
  ab = c(i,b)
  df_need = rbind(df_need, ab)
}

names(df_need) = c('species_name', 'GhAhSkew')
df_need$species_name = gsub(' ', '_', df_need$species_name)
df_b_cytb = merge(df_need, cytb, by = 'species_name')
df_b_nd2 = merge(df_need, nd2, by = 'species_name')
df_b_all = merge(df_need, df_midori1, by = 'species_name')
unique(df_b_all$species_name)
