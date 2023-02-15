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
library(data.table)
pca_data = mut_data_ff1[,c(1,5,7,8,9)]

pca_data_shaped = dcast.data.table(setDT(pca_data), species_name+realm+Trophic_niche~Mut,
                 value.var='MutSpec')
pca_data_shaped$`A>C` = as.numeric(pca_data_shaped$`A>C`)
pca_data_shaped$`A>G` = as.numeric(pca_data_shaped$`A>G`)
pca_data_shaped$`A>T` = as.numeric(pca_data_shaped$`A>T`)
pca_data_shaped$`C>A` = as.numeric(pca_data_shaped$`C>A`)
pca_data_shaped$`C>G` = as.numeric(pca_data_shaped$`C>G`)
pca_data_shaped$`C>T` = as.numeric(pca_data_shaped$`C>T`)
pca_data_shaped$`G>A` = as.numeric(pca_data_shaped$`G>A`)
pca_data_shaped$`G>C` = as.numeric(pca_data_shaped$`G>C`)
pca_data_shaped$`G>T` = as.numeric(pca_data_shaped$`G>T`)
pca_data_shaped$`T>A` = as.numeric(pca_data_shaped$`T>A`)
pca_data_shaped$`T>C` = as.numeric(pca_data_shaped$`T>C`)
pca_data_shaped$`T>G` = as.numeric(pca_data_shaped$`T>G`)

stats_pca = prcomp(pca_data_shaped[,c(4,5,6,7,8,9,10,11,12,13,14,15)], center = TRUE, scale. = TRUE)
summary(stats_pca)

bipl = ggbiplot(stats_pca, groups = pca_data_shaped$realm, labels = pca_data_shaped$species_name, labels.size = 2)
bipl
ggplotly(bipl)

bipl1 = ggbiplot(stats_pca, groups = pca_data_shaped$Trophic_niche, labels = pca_data_shaped$species_name, labels.size = 2)
bipl1
ggplotly(bipl1)

bipl2 = ggbiplot(stats_pca, groups = pca_data_shaped$realm, labels.size = 2)
bipl2
ggplotly(bipl2)

bipl3 = ggbiplot(stats_pca, groups = pca_data_shaped$Trophic_niche, labels.size = 2)
bipl3
ggplotly(bipl3)


#Valya's data
valya_data = read.csv('../../Head/2Scripts/valyadata_final.csv')
valya_data = na.omit(valya_data)
valya_data = valya_data[,c(1,3,4,5,6)]
valya_data["species_name"][valya_data["species_name"] == "Strigops_habroptilus"] = "Strigops_habroptila"
valya_gene = merge(pca_data_shaped, valya_data, by = 'species_name')
stats_pca1 = prcomp(valya_gene[,c(4,5,6,7,8,9,10,11,12,13,14,15)], center = TRUE, scale. = TRUE)
summary(stats_pca1)
bipl_valya1 = ggbiplot(stats_pca1, groups = valya_gene$far_migration, labels = valya_gene$species_name, labels.size = 2)
bipl_valya1
ggplotly(bipl_valya1)

bipl_valya2 = ggbiplot(stats_pca1, groups = valya_gene$far_migration, labels.size = 2)
bipl_valya2
ggplotly(bipl_valya2)

bipl_valya3 = ggbiplot(stats_pca1, groups = valya_gene$wintering, labels = valya_gene$species_name, labels.size = 2)
bipl_valya3
ggplotly(bipl_valya3)

bipl_valya4 = ggbiplot(stats_pca1, groups = valya_gene$wintering, labels.size = 2)
bipl_valya4
ggplotly(bipl_valya4)

bipl_valya5 = ggbiplot(stats_pca1, groups = valya_gene$diving, labels = valya_gene$species_name, labels.size = 2)
bipl_valya5
ggplotly(bipl_valya5)

bipl_valya6 = ggbiplot(stats_pca1, groups = valya_gene$diving, labels.size = 2)
bipl_valya6
ggplotly(bipl_valya6)

bipl_valya7 = ggbiplot(stats_pca1, groups = valya_gene$flying, labels = valya_gene$species_name, labels.size = 2)
bipl_valya7
ggplotly(bipl_valya7)

bipl_valya8 = ggbiplot(stats_pca1, groups = valya_gene$flying, labels.size = 2)
bipl_valya8
ggplotly(bipl_valya8)
