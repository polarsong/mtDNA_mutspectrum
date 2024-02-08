rm(list = ls(all=TRUE))
df_flighless = read.csv('Flightless_birds_from_wiki.csv', header = F)
names(df_flighless) = c('species_name')
DNA <- read.csv2("Birds_dataset_paper.csv", header=T)
DNA <- DNA[c('Species', 'ghahSkew', 'chthSkew')]
DNA <- aggregate (DNA[,-1], by=list(DNA$Species), mean)
names(DNA) = c('species_name', 'GhAhSkew', 'ThChSkew')
df_check = merge(df_flighless,DNA)
df_fly = read.csv('../../Paper_materials_2024/flying_birds.csv')
df_fly = na.omit(df_fly)
df_fly_check = df_fly[df_fly$Flightless != 0,]
df_fly = df_fly[,c(2,3)]
names(df_fly) = c('species_name', 'flightless')
df_check1 = merge(DNA, df_fly)
df_check2 = merge(df_check1, df_flighless)
