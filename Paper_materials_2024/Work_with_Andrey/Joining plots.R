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
                      'gene', 'Species3')
df_midori1$Species3 = gsub('_', ' ', df_midori1$Species3)
df_fly = read.csv('../../Paper_materials_2024/flying_birds.csv')
df_fly = df_fly[,c(2,3,4)]
names(df_fly) = c('Species3', 'ability_to_fly', 'ability_to_dive')
avo_data = read.csv('../../Body/1Raw/Avonet_data.csv')
names_v = unique(df_mtdna$species_name)
df_short = data.frame()
for (i in names_v)
{
  df1 = df_mtdna[df_mtdna$species_name == i,]
  a = sum(df1$ghahSkew)/12
  b = sum(df1$chthSkew)/12
    ab = c(i, a, b)
  df_short = rbind(df_short, ab)
}
names(df_short) = c('Species3', 'GhAhSkew', 'ThChSkew')
df_fly = na.omit(df_fly)
df_midori_birds = merge(df_midori1, avo_data)
names_v1 = unique(df_midori_birds$Species3)
df_short1 = data.frame()
for (i in names_v1)
{
  df1 = df_midori_birds[df_midori_birds$Species3 == i,]
  a = (sum(df1$Mutation_AG_midori))/(nrow(df1))
  b = (sum(df1$Mutation_CT_midori))/(nrow(df1))
  ab = c(i, a, b)
  df_short1 = rbind(df_short1, ab)
}
names(df_short1) = c('Species3', 'Mutation_AG_midori', 'Mutation_CT_midori')
df_cytb = read.csv('Midori2_birds_cytb_ghahskew_better.csv')
library(dplyr)
df_big = full_join(df_short, df_fly)
df_cytb = df_cytb[,c(2,3,4)]
names(df_cytb) = c('Species3', 'GhAhSkew_seq_beg', 'GhAhSkew_start_codon')
#df_midori_birds = merge(df_midori1, avo_data, by = 'Species3')
df_big1 = full_join(df_big, df_short1)
df_big11 = full_join(df_big1, df_cytb)
df_big2 = full_join(df_big11, avo_data)
write.csv(df_big2, file = 'big_birds_new_data.csv')
