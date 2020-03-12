###################################
###### 
###################################

rm(list=ls(all=TRUE))

############ list of hibernating animals
Hib = read.table("../../Body/1Raw/HibernatingDailytorporMammals.txt")
Hib$V1 = paste(Hib$V1, Hib$V2, sep="_")
ListOfHibSpecies = Hib[Hib$V3 == "HIB",]$V1; length(ListOfHibSpecies)
ListOfNHibSpecies = Hib[!Hib$V3 == "HIB",]$V1; length(ListOfNHibSpecies)
############ Syn mut
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")) {file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")}


### make ND6 complementary:
NotND6 = SynNuc[SynNuc$Gene != 'ND6',]
NotND6$FrA = NotND6$NeutralA / (NotND6$NeutralA + NotND6$NeutralT + NotND6$NeutralG + NotND6$NeutralC)
NotND6$FrT = NotND6$NeutralT / (NotND6$NeutralA + NotND6$NeutralT + NotND6$NeutralG + NotND6$NeutralC) 
NotND6$FrG = NotND6$NeutralG / (NotND6$NeutralA + NotND6$NeutralT + NotND6$NeutralG + NotND6$NeutralC) 
NotND6$FrC = NotND6$NeutralC / (NotND6$NeutralA + NotND6$NeutralT + NotND6$NeutralG + NotND6$NeutralC) 

ND6 = SynNuc[SynNuc$Gene == 'ND6',]
ND6$FrA = ND6$NeutralT / (ND6$NeutralA + ND6$NeutralT + ND6$NeutralG + ND6$NeutralC)
ND6$FrT = ND6$NeutralA / (ND6$NeutralA + ND6$NeutralT + ND6$NeutralG + ND6$NeutralC) 
ND6$FrG = ND6$NeutralC / (ND6$NeutralA + ND6$NeutralT + ND6$NeutralG + ND6$NeutralC) 
ND6$FrC = ND6$NeutralG / (ND6$NeutralA + ND6$NeutralT + ND6$NeutralG + ND6$NeutralC) 

SynNuc = rbind(NotND6,ND6)

VecOfTaxa = unique(SynNuc$Class)

############# Hibernation in different classes
#### who is in our list? only mammals!
for (i in 1:length(VecOfTaxa)){ # i = 1
  Species = unique(SynNuc[SynNuc$Class == VecOfTaxa[i],]$Species); length(Species)
  L = length(intersect(Species,ListOfHibSpecies))
  print(c(as.character(VecOfTaxa[i]), L))
}

########## compare nucleotide frequencies between hibernating and other mammals
SynNuc = SynNuc[SynNuc$Class == 'Mammalia',]; length(unique(SynNuc$Species))
AGG = aggregate(list(SynNuc$FrA,SynNuc$FrT,SynNuc$FrG,SynNuc$FrC), by = list(SynNuc$Species), FUN = mean)
names(AGG)=c('Species','FrA','FrT','FrG','FrC')


#pdf("../../Body/4Figures/WholeGenomeAnalyses.HibernationWithinMammals.NoOverlap.R.01.pdf")

par(mfrow=c(2,2))
boxplot(AGG[AGG$Species %in% ListOfHibSpecies,]$FrA,AGG[!AGG$Species %in% ListOfHibSpecies,]$FrA, notch = TRUE, main = 'A', names = c('Hib','Others'))
boxplot(AGG[AGG$Species %in% ListOfHibSpecies,]$FrT,AGG[!AGG$Species %in% ListOfHibSpecies,]$FrT, notch = TRUE, main = 'T', names = c('Hib','Others'))
boxplot(AGG[AGG$Species %in% ListOfHibSpecies,]$FrG,AGG[!AGG$Species %in% ListOfHibSpecies,]$FrG, notch = TRUE, main = 'G', names = c('Hib','Others'))
boxplot(AGG[AGG$Species %in% ListOfHibSpecies,]$FrC,AGG[!AGG$Species %in% ListOfHibSpecies,]$FrC, notch = TRUE, main = 'C', names = c('Hib','Others'))

dev.off()

### TO DO
# hibernation vs torpor
# effect of generation time (hib are short or long lived...)
# temp of envir of cold-blooded

############# Longevity 
AA = read.table("../../Body/1Raw/anage_data.txt", header = TRUE, sep = '\t')
AA$Species = paste(AA$Genus,AA$Species,sep = '_')

GL=read.table("../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt", header = TRUE, sep = '\t')
GL$Species = gsub(" ", "_", GL$Scientific_name)

vec_of_Marsupials_orders = c("Dasyuromorphia", 'Didelphimorphia', "Diprotodontia", "Microbiotheria", "Notoryctemorphia", "Paucituberculata", "Peramelemorphia")
Vec_of_Monotremata_genus = c("Tachyglossus","Zaglossus", "Ornithorhynchus") #Ornithorhynchus
Vec_of_Marsupials_genus = AA[AA$Order %in% vec_of_Marsupials_orders,]$Genus
vec_of_NonWarm = c(GL[GL$Genus %in% Vec_of_Monotremata_genus,]$Species, GL[GL$Genus %in% Vec_of_Marsupials_genus,]$Species, "Heterocephalus_glaber")
Vec_of_HibDai = c(ListOfNHibSpecies, ListOfHibSpecies)
GL$Hib = 1
GL[!GL$Species %in% ListOfHibSpecies,]$Hib = 0
GL$Daily = 1
GL[!GL$Species %in% ListOfNHibSpecies,]$Daily = 0
GL$Mono= 1
GL[!GL$Genus %in% Vec_of_Monotremata_genus,]$Mono = 0
GL$Mars = 1
GL[!GL$Genus %in% Vec_of_Marsupials_genus,]$Mars = 0
GL$NonWarm = 1
GL[!GL$Species %in% vec_of_NonWarm,]$NonWarm = 0
GL$HibDai = 1
GL[!GL$Species %in% Vec_of_HibDai,]$HibDai = 0


allparameters = merge(GL, AGG)
table(allparameters$Hib)
table(allparameters$Daily)
table(allparameters$Mono)
table(allparameters$Mars)
table(allparameters$NonWarm)
table(allparameters$HibDai)

allparameters$Temper = "Warm"

allparameters[allparameters$Species %in% Vec_of_HibDai, ]$Temper= "HibDay"
allparameters[allparameters$Species %in% vec_of_NonWarm, ]$Temper= "MonoMarsHG"


ltest = lm(formula = FrT ~ scale(GenerationLength_d)*scale(HibDai), data = allparameters)
summary(ltest)


ltest = lm(formula = FrT ~ scale(GenerationLength_d)+scale(HibDai), data = allparameters)
summary(ltest)


ltest = lm(formula = scale(FrT) ~ 0 + scale(GenerationLength_d)+scale(HibDai), data = allparameters)
summary(ltest)


library(ggpubr)
ggscatter(allparameters, x = "GenerationLength_d", y = "FrT", color = "Temper",
          palette = c("#00AFBB", "#756bb1", "#FC4E07"), add = "reg.line",  xscale = "log2", cor.coeff.args = list(method = "spearman"))
ggscatter(allparameters, x = "GenerationLength_d", y = "FrT",
          color = "Temper", shape = "Temper",
          palette = c("#00AFBB", "#756bb1", "#FC4E07"),
          ellipse = TRUE, mean.point = TRUE, add = "reg.line",  xscale = "log2")


gghistogram(allparameters, x = "GenerationLength_d",
            add = "mean", rug = TRUE,
            color = "Temper", fill = "Temper",
            palette = c("#00AFBB", "#756bb1", "#FC4E07"), bins = 50)

alleqgenl = allparameters[allparameters$GenerationLength_d < 2000,]
alleqgenl = allparameters[allparameters$GenerationLength_d < 2500,]
alleqgenl = allparameters[allparameters$GenerationLength_d > 1250,]

alleqgenl[alleqgenl$Temper == "MonoMarsHG",]

ggboxplot(alleqgenl, x = "Temper", y = "FrT",
               color = "Temper", palette =c("#FC4E07", "#00AFBB", "#756bb1"),
               add = "jitter", shape = "Temper")
nrow(alleqgenl[alleqgenl$Temper == "Warm",])
