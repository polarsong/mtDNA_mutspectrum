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

ggplot(data = mut_data_ff, aes(x = Mut, y = MutSpec, fill = Trophic_niche))+
  geom_boxplot()

fig <- plot_ly(mut_data_ff, x = ~Mut, y = ~MutSpec, color = ~Trophic_niche, type = "box")
fig <- fig %>% layout(boxmode = "group")

fig
