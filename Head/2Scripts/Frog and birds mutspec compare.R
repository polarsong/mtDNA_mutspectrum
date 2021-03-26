library('ggplot2')
rm(list=ls(all=TRUE))
library('heatmaply')
#frogs
mutspec = read.table('VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt')
classes = read.table('TaxaWithClasses.txt')
only_anura = data.frame(classes$Species[classes$Class == 'Amphibia'])
names(only_anura) = c('Anura_species')
UV_names = read.table('Names.txt')
UV_clutch = read.table('Clutch.txt')
UV_characteristic = read.table('UV.txt')
unfull = cbind(UV_names, UV_clutch)
anura_with_UV = cbind(unfull, UV_characteristic)
names(anura_with_UV) = c('Name', 'Type_of_clutch', 'UV_influence')
names(mutspec) = c('Name', "A_T", 'A_G', "A_C", 'T_A', "T_G", 'T_C', 'G_A', 'G_T', 'G_C', 'C_A', 'C_T', 'C_G')
anura_with_UV_and_mutspec = merge(anura_with_UV, mutspec, by = c('Name'))
MutSpecAnura = mutspec[mutspec$Name %in% only_anura$Anura_species,]
AmphibiaWithoutAnnotation = setdiff(MutSpecAnura$Name,anura_with_UV$Name) #найти что-то интересное!!!!!
MutSpecAnura = MutSpecAnura[order(MutSpecAnura$T_C),] #10-15 сверху и снизу
3168/12

#1. ƒонайти экологию
#birds
#389 просмотреть глазами сравнить названи€!!!!!
#7 признаков, 12 фенотипов, 
#¬сех птиц из мутспека сделать pdf!!!!!
#сортировать с_т, A_G - heavy chain!!!!!
#цветные  сосиски!!!!!
#heatmap - изучить!!!! - высокий уровень приоритета
#https://www.datanovia.com/en/blog/how-to-create-a-beautiful-interactive-heatmap-in-r/
birds_pheno = read.table('ALL_PHENOTYPES.txt')
names(birds_pheno) = c('Name', 'Phenotype')
birds_with_pheno_and_mutspec = merge(mutspec, birds_pheno, by = c("Name"))
table(birds_with_pheno_and_mutspec$Phenotype)
ggplot(data = birds_with_pheno_and_mutspec, aes(colour = Phenotype, x = Name, y = G_A)) +
  geom_point()+
  theme(axis.text.x = element_text(angle = 90))
ggplot(data = birds_with_pheno_and_mutspec, aes(colour = Phenotype, x = Name, y = G_C)) +
  geom_point()+
  theme(axis.text.x = element_text(angle = 90))
ggplot(data = birds_with_pheno_and_mutspec, aes(colour = Phenotype, x = Name, y = G_T)) +
  geom_point()+
  theme(axis.text.x = element_text(angle = 90))
ggplot(data = birds_with_pheno_and_mutspec, aes(colour = Phenotype, x = Name, y = A_G)) +
  geom_point()+
  theme(axis.text.x = element_text(angle = 90))             
ggplot(data = birds_with_pheno_and_mutspec, aes(colour = Phenotype, x = Name, y = A_C)) +
  geom_point()+
  theme(axis.text.x = element_text(angle = 90))
ggplot(data = birds_with_pheno_and_mutspec, aes(colour = Phenotype, x = Name, y = A_T)) +
  geom_point()+
  theme(axis.text.x = element_text(angle = 90))
ggplot(data = birds_with_pheno_and_mutspec, aes(colour = Phenotype, x = Name, y = T_G)) +
  geom_point()+
  theme(axis.text.x = element_text(angle = 90))
ggplot(data = birds_with_pheno_and_mutspec, aes(colour = Phenotype, x = Name, y = T_C)) +
  geom_point()+
  theme(axis.text.x = element_text(angle = 90))
ggplot(data = birds_with_pheno_and_mutspec, aes(colour = Phenotype, x = Name, y = T_A)) +
  geom_point()+
  theme(axis.text.x = element_text(angle = 90))
ggplot(data = birds_with_pheno_and_mutspec, aes(colour = Phenotype, x = Name, y = C_G)) +
  geom_point()+
  theme(axis.text.x = element_text(angle = 90))
ggplot(data = birds_with_pheno_and_mutspec, aes(colour = Phenotype, x = Name, y = C_T)) +
  geom_point()+
  theme(axis.text.x = element_text(angle = 90))
ggplot(data = birds_with_pheno_and_mutspec, aes(colour = Phenotype, x = Name, y = C_A)) +
  geom_point()+
  theme(axis.text.x = element_text(angle = 90))
ggplot(birds_with_pheno_and_mutspec, aes(fill=Phenotype, y=A_T, x=Name)) + 
  geom_bar(position="fill", stat="identity")+
  theme(axis.text.x = element_text(angle = 90))
heatmap(one_bird)
normalize(birds_with_pheno_and_mutspec)

birds_with_pheno_and_mutspec = transform(birds_with_pheno_and_mutspec, A_T = as.numeric(A_T))
birds_with_pheno_and_mutspec = transform(birds_with_pheno_and_mutspec, A_G = as.numeric(A_G))
birds_with_pheno_and_mutspec = transform(birds_with_pheno_and_mutspec, A_C = as.numeric(A_C))
birds_with_pheno_and_mutspec = transform(birds_with_pheno_and_mutspec, G_T = as.numeric(G_T))
birds_with_pheno_and_mutspec = transform(birds_with_pheno_and_mutspec, G_A = as.numeric(G_A))
birds_with_pheno_and_mutspec = transform(birds_with_pheno_and_mutspec, G_C = as.numeric(G_C))
birds_with_pheno_and_mutspec = transform(birds_with_pheno_and_mutspec, C_T = as.numeric(C_T))
birds_with_pheno_and_mutspec = transform(birds_with_pheno_and_mutspec, C_G = as.numeric(C_G))
birds_with_pheno_and_mutspec = transform(birds_with_pheno_and_mutspec, C_A = as.numeric(C_A))
birds_with_pheno_and_mutspec = transform(birds_with_pheno_and_mutspec, T_A = as.numeric(T_A))
birds_with_pheno_and_mutspec = transform(birds_with_pheno_and_mutspec, T_C = as.numeric(T_C))
birds_with_pheno_and_mutspec = transform(birds_with_pheno_and_mutspec, T_G = as.numeric(T_G))


to_change = c(birds_with_pheno_and_mutspec$Name)
to_change[31] = 'Pygoscelis_papua1'
to_change[38] = 'Tetraogallus_himalayensis1'
row.names(birds_draw, ) = to_change
birds_draw$Name = NULL



birds_draw = normalize(birds_with_pheno_and_mutspec)
heatmaply(birds_draw, file = "birds heatmaps/bird mutspec.html")
phe_vec = birds_draw$Phenotype
birds_draw$Phenotype = NULL
birds_draw = round(birds_draw, digits = 2)
birds_draw$Phenotype = phe_vec

#frog drawing
anura = anura_with_UV_and_mutspec
row.names(anura) = anura$Name
anura$Name = NULL
anura$Type_of_clutch = NULL
anura1 = normalize(anura)
typeof(anura1$A_T)
heatmaply(anura1)
heatmaply(ex1, file = "folder/f.html")
anura$A_T = type.convert(anura1$A_T)
typeof(f1)
anura = transform(anura, A_T = as.numeric(A_T))
anura = transform(anura, A_G = as.numeric(A_G))
anura = transform(anura, A_C = as.numeric(A_C))
anura = transform(anura, G_T = as.numeric(G_T))
anura = transform(anura, G_A = as.numeric(G_A))
anura = transform(anura, G_C = as.numeric(G_C))
anura = transform(anura, C_T = as.numeric(C_T))
anura = transform(anura, C_G = as.numeric(C_G))
anura = transform(anura, C_A = as.numeric(C_A))
anura = transform(anura, T_A = as.numeric(T_A))
anura = transform(anura, T_C = as.numeric(T_C))
anura = transform(anura, T_G = as.numeric(T_G))
row.names(anura) = anura$Name
anura$Name = NULL
anura$Type_of_clutch = NULL
anura1 = normalize(anura)
dir.create("frog heatmaps")
heatmaply(anura1, file = "frog heatmaps/frog mutspec.html")
