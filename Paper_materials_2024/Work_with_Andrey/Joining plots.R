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
df_midori1$species_name = gsub('_', ' ', df_midori1$species_name)
df_fly = read.csv('../../Paper_materials_2024/flying_birds.csv')
df_fly = df_fly[,c(2,3,4)]
names(df_fly) = c('species_name', 'ability_to_fly', 'ability_to_dive')
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
names(avo_data)[names(avo_data) == 'Species3'] <- 'species_name'
names(avo_data)[names(avo_data) == 'Family3'] <- 'family'
names(avo_data)[names(avo_data) == 'Order3'] <- 'order'
names(df_short) = c('species_name', 'GhAhSkew_refseq', 'ThChSkew_refseq')
df_fly = na.omit(df_fly)
df_midori_birds = merge(df_midori1, avo_data)
names_v1 = unique(df_midori_birds$species_name)
df_short1 = data.frame()
for (i in names_v1)
{
  df1 = df_midori_birds[df_midori_birds$species_name == i,]
  a = (sum(df1$Mutation_AG_midori))/(nrow(df1))
  b = (sum(df1$Mutation_CT_midori))/(nrow(df1))
  ab = c(i, a, b)
  df_short1 = rbind(df_short1, ab)
}
names(df_short1) = c('species_name', 'Mutation_AG_midori', 'Mutation_CT_midori')
df_cytb = read.csv('Midori2_birds_cytb_ghahskew_better.csv')
df_syst = read.csv('TaxonomyIOC14_1.csv', sep = ';')
df_midori_eco = read.csv('Midori_birds_eco.csv')
df_midori_eco = df_midori_eco[,c(2,3,4)]
df_midori_eco$species_name = gsub('_',' ', df_midori_eco$species_name)
names(df_midori_eco) = c('species_name', 'ability_to_fly_midori', 'ability_to_dive_midori')
names(df_syst)[names(df_syst) == 'species_name'] <- 'species_name_IOC'
names(df_syst) = c('order', 'family', 'genus', 'species_name_IOC')
df_syst = df_syst[,c(1,4)]
library(dplyr)
df_big = full_join(df_short, df_fly)
df_cytb = df_cytb[,c(2,3,4)]
names(df_cytb) = c('species_name', 'GhAhSkew_seq_beg', 'GhAhSkew_start_codon')
#df_midori_birds = merge(df_midori1, avo_data, by = 'Species3')
df_big1 = full_join(df_big, df_short1)
df_big11 = full_join(df_big1, df_cytb)
df_big2 = full_join(df_big11, avo_data)
df_big3 = full_join(df_big2, df_midori_eco)
write.csv(df_big3, file = 'big_birds_new_data.csv')

#Vienn diagram
library(devtools)
devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)
df_venn = df_big3[,c(1,2,3,4,5,6,7,8,9,45,46)]
df_venn = df_venn[,c(2:11)]
ggvenn(
  df_venn
)
