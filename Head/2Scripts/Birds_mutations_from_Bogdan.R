rm(list=ls(all=TRUE))
library(ggplot2)
library(plotly)
#four-fold mutations
df_mtdna = read.csv('../../Head/2Scripts/Birds_dataset_paper.csv')
mut_data = read.table("C:/Users/User/Desktop/Birds mutspec results from Bogdan/mutspec12.tsv", header = TRUE, fill = TRUE)
unique(mut_data$AltNode)
mut_data_ff = mut_data[mut_data$Label == 'ff',]
mut_data_ff = mut_data_ff[,c(1,2,3,4,5,7,8)]
ecozone_data = df_mtdna[,c('species_name', 'realm', 'Trophic_niche')]
ecozone_data = unique(ecozone_data)
ecozone_data$species_name = gsub(' ', '_', ecozone_data$species_name)
mut_data_ff = mut_data_ff[!grepl('Node', mut_data_ff$AltNode),]
names(mut_data_ff) = c('Mut', 'ObsNum', 'ExpNum', 'RawMutSpec', 'MutSpec', 'species_name', 'Label') 
mut_data_ff = merge(mut_data_ff, ecozone_data, by = 'species_name')
mut_data_ff$MutSpec =  as.numeric(mut_data_ff$MutSpec)

AC = mut_data_ff[mut_data_ff$Mut == 'A>C',] 
AG = mut_data_ff[mut_data_ff$Mut == 'A>G',]
AT = mut_data_ff[mut_data_ff$Mut == 'A>T',]
GC = mut_data_ff[mut_data_ff$Mut == 'G>C',]
GT = mut_data_ff[mut_data_ff$Mut == 'G>T',]
GA = mut_data_ff[mut_data_ff$Mut == 'G>A',]
CG = mut_data_ff[mut_data_ff$Mut == 'C>G',]
CT = mut_data_ff[mut_data_ff$Mut == 'C>T',]
CA = mut_data_ff[mut_data_ff$Mut == 'C>A',]
TG = mut_data_ff[mut_data_ff$Mut == 'T>G',]
TC = mut_data_ff[mut_data_ff$Mut == 'T>C',]
TA = mut_data_ff[mut_data_ff$Mut == 'T>A',]

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

mut_data_ff1 = rbind(AC, AG, AT, GC, GT, GA, CT, CA, CG, TA, TG, TC)


ggplot(data = mut_data_ff, aes(x = Mut, y = MutSpec, fill = Trophic_niche))+
  geom_boxplot()

fig <- plot_ly(mut_data_ff1, x = ~Mut, y = ~MutSpec, color = ~realm, type = "box")
fig <- fig %>% layout(boxmode = "group")

fig

#PCA work
library(ggbiplot)
pca_data = mut_data_ff1[,c(1,5,7,8,9)]
stats_pca = prcomp(pca_data[c(2, 5)], center = TRUE, scale. = TRUE)
summary(stats_pca)


#192 mutspec
mut_data_192 = read.table("C:/Users/User/Desktop/Birds mutspec results from Bogdan/mutspec192.tsv", header = TRUE, fill = TRUE)
mut_data_192_ff = mut_data_192[mut_data_192$Label == 'ff',]
mut_data_192_ff = mut_data_192_ff[,c(1,2,3,4,5,7,8)]

mut_data_192_ff = mut_data_192_ff[!grepl('Node', mut_data_192_ff$AltNode),]
names(mut_data_192_ff) = c('Mut', 'ObsNum', 'ExpNum', 'RawMutSpec', 'MutSpec', 'species_name', 'Label') 
mut_data_192_ff = merge(mut_data_192_ff, ecozone_data, by = 'species_name')
mut_data_192_ff$MutSpec =  as.numeric(mut_data_192_ff$MutSpec)
fig <- plot_ly(mut_data_192_ff, x = ~Mut, y = ~MutSpec, color = ~realm, type = "box")
fig <- fig %>% layout(boxmode = "group")

fig
