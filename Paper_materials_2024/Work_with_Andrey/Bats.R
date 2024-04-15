rm(list = ls(all=TRUE))
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip", exdir = "../../Body/3Results/")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")) file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")
names(SynNuc)
SynNuc$ghahSkew = ((SynNuc$NeutralC - SynNuc$NeutralT))/((SynNuc$NeutralC + SynNuc$NeutralT))
SynNuc$chthSkew = ((SynNuc$NeutralA - SynNuc$NeutralG))/((SynNuc$NeutralA + SynNuc$NeutralG))
new_mam = SynNuc[, c(1, 2, 79, 80)]
SynNuc$Taxonomy
library(data.table)
df1 = SynNuc[SynNuc$Taxonomy %like% "Mammalia", ]
df_bats = SynNuc[SynNuc$Taxonomy %like% "Chiroptera", ]
names(df_bats)
df_bats = df_bats[,c(1,79)]
df_bats = aggregate(.~Species, FUN = mean, data = df_bats)
df_bats$animal = 'Bats'
df_mammals = df1[,c(1,79)]
df_mammals = aggregate(.~Species, FUN = mean, data = df_mammals)
df_mammals = df_mammals[!df_mammals$Species %in% df_bats$Species,]
df_mammals$animal = 'Mammals'
df_mtdna = read.csv('../../Paper_materials_2024/Birds_dataset_paper.csv')
df_ex = df_mtdna[,c(2,97)]
df_ex1 = aggregate(. ~ species_name, FUN = mean, data = df_ex)
df_ex1$animal = 'Birds'
library(ggplot2)
a = ggplot(data=df_bats, aes(x = animal, y = ghahSkew))+
  geom_boxplot()+
  ylim(-0.4,1)
b = ggplot(data=df_mammals, aes(x = animal, y = ghahSkew))+
  geom_boxplot()+
  ylim(-0.4,1)
c = ggplot(data=df_ex1, aes(x = animal, y = ghahSkew))+
  geom_boxplot()+
  ylim(-0.4,1)
library(ggpubr)
animals = ggarrange(a, b, c,
                    ncol = 3, nrow = 1)
animals
