library('ggplot2')
install.packages("heatmaply")
library('heatmaply')
rm(list=ls(all=TRUE))
#frog
codus = read.table('AllGenesCodonUsageNoOverlap.txt', header = TRUE, sep = '\t')
codus$ratioA = codus$NeutralA/(codus$NeutralA+codus$NeutralG+codus$NeutralC+codus$NeutralT)
codus$ratioT = codus$NeutralT/(codus$NeutralA+codus$NeutralG+codus$NeutralC+codus$NeutralT)
codus$ratioC = codus$NeutralC/(codus$NeutralA+codus$NeutralG+codus$NeutralC+codus$NeutralT)
codus$ratioG = codus$NeutralG/(codus$NeutralA+codus$NeutralG+codus$NeutralC+codus$NeutralT)
species_and_ratio = data.frame(codus$Species)
species_and_ratio$ratioA = codus$ratioA
species_and_ratio$ratioG = codus$ratioG
species_and_ratio$ratioC = codus$ratioC
species_and_ratio$ratioT = codus$ratioT
names(species_and_ratio) = c('Name', 'RatioA', 'RatioG', 'RatioC', 'RatioT')
species_and_ratio$Gene = codus$Gene
species_and_ratio$GeneStart = codus$GeneStart
species_and_ratio$GeneEnd = codus$GeneEnd
species_and_ratio$Taxonomy = codus$Taxonomy
species_and_ratio$Class = codus$Class
frog = read.table('Names.txt')
clutch = read.table('Clutch.txt')
uv = read.table('UV.txt')
frog = cbind(frog, clutch)
frog = cbind(frog, uv)
names(frog) = c('Name', 'Clutch', 'UV')
frog_with_ratio = merge(species_and_ratio, frog, by = 'Name' )
frog_with_ratio = frog_with_ratio[frog_with_ratio$Gene != 'ND6',]
ggplot(data = frog_with_ratio, aes(col = UV, x = Name, y = RatioA))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))
ggplot(data = frog_with_ratio, aes(col = UV, x = Name, y = RatioG))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))
ggplot(data = frog_with_ratio, aes(col = UV, x = Name, y = RatioC))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))
ggplot(data = frog_with_ratio, aes(col = UV, x = Name, y = RatioT))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))
#bird
#A_t, G_C сосиски, heatmap, интерактивный HTML-файл - почитать, как сохранить!!!!!
pheno = read.table("ALL_PHENOTYPES.txt")
names(pheno) = c('Name', 'Phenotype')
birds_with_ratio = merge(species_and_ratio, pheno, by = 'Name')
birds_with_ratio = birds_with_ratio[birds_with_ratio$Gene != 'ND6',]
one_bird = birds_with_ratio[c(1,2,3,4,5,6,7,8,9,10,11,12),]
one_bird$count = one_bird$RatioA + one_bird$RatioT + one_bird$RatioC + one_bird$RatioG
ggplot(data = one_bird, aes(x = Gene, y = count, fill = RatioA + RatioC))+
  geom_bar(position = 'fill', stat = 'identity')
ggplot(data = birds_with_ratio, aes(col = Phenotype, x = Name, y = RatioA))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))
ggplot(data = birds_with_ratio, aes(col = Phenotype, x = Name, y = RatioG))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))
ggplot(data = birds_with_ratio, aes(col = Phenotype, x = Name, y = RatioC))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))
ggplot(data = birds_with_ratio, aes(col = Phenotype, x = Name, y = RatioT))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))

ggplot(birds_with_ratio, aes(fill=Phenotype, y=RatioA, x=Name)) + 
  geom_bar(position="fill", stat="identity")+
  theme(axis.text.x = element_text(angle = 90))
df = normalize(one_bird)
heatmaply(df)
birds_with_ratio_for_heatmap = birds_with_ratio

birds_with_ratio_for_heatmap$GeneStart = NULL
birds_with_ratio_for_heatmap$GeneEnd = NULL
birds_with_ratio_for_heatmap$Taxonomy = NULL
birds_with_ratio_for_heatmap$Class = NULL
ex = normalize(birds_with_ratio_for_heatmap)
ex = birds_with_ratio[c(1,2,3,4,5,6,7,8,9,10,11,12),]
birds_with_ratio$Gene=NULL
birds_with_ratio$GeneStart=NULL
birds_with_ratio$GeneEnd=NULL
birds_with_ratio$Taxonomy=NULL
birds_with_ratio$Class=NULL


#drawing heatmaps
dir.create("birds heatmaps")
row_vector = c(1,2,3,4,5,6,7,8,9,10,11,12)
row_vector = row_vector + 12
ex = birds_with_ratio[row_vector,]
row.names(ex) = ex$Gene
ex$Gene = NULL
ex$Name
ex$Phenotype
ex1 = normalize(ex)
heatmaply(ex1)
heatmaply(ex1, file = "folder/f.html")
browseURL("folder/heatmaply_plot.html")
