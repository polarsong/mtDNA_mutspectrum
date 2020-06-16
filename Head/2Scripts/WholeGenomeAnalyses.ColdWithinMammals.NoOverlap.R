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

kuptsovtable = read.table("../../Body/2Derived/EcologyMammalianTable01_KuptsovA_ver2_Full.txt", sep='\t', header=TRUE)
  

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


table(allparameters$Hib)
table(allparameters$Daily)
table(allparameters$Mono)
table(allparameters$Mars)
table(allparameters$NonWarm)
table(allparameters$HibDai)

allparameters$Temperature = "Warm"

allparameters[allparameters$Species %in% Vec_of_HibDai, ]$Temperature= "HibDay"
allparameters[allparameters$Species %in% vec_of_NonWarm, ]$Temperature= "MonoMarsHG"


ltest = lm(formula = FrT ~ scale(GenerationLength_d), data = allparameters)
summary(ltest)

allparameters$residuals = ltest$residuals


ltest = lm(formula = FrT ~ scale(GenerationLength_d)+scale(HibDai), data = allparameters)
summary(ltest)


ltest = lm(formula = scale(FrT) ~ 0 + scale(GenerationLength_d)+scale(HibDai), data = allparameters)
summary(ltest)


library(ggpubr)
ggscatter(allparameters, x = "GenerationLength_d", y = "FrT", color = "Temperature",
          palette = c("#00AFBB", "#756bb1", "#FC4E07"), add = "reg.line",  xscale = "log2", cor.coeff.args = list(method = "spearman"))
ggscatter(allparameters, x = "GenerationLength_d", y = "FrT",
          color = "Temperature", shape = "Temperature",
          palette = c("#08519c", "#3182bd", "#fc9272"),
          ellipse = TRUE, mean.point = TRUE, add = "reg.line",  xscale = "log2", xlab="Generation Length, log2", ylab="Fraction of A")


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


allparameters = merge(allparameters, AA, all.x = TRUE)
EcologyMammalianTable= data.frame(allparameters$Kingdom, allparameters$Phylum, allparameters$Class, allparameters$Order, allparameters$Family,allparameters$Genus, allparameters$Species, allparameters$GenerationLength_d, allparameters$FrT, allparameters$Hib, allparameters$Daily, allparameters$Mars, allparameters$Mono, allparameters$residuals)
EcologyMammalianTable= EcologyMammalianTable[order(EcologyMammalianTable$allparameters.residuals),]
write.table(EcologyMammalianTable, file="../../Body/2Derived/EcologyMammalianTable01.txt", quote = FALSE, row.names = FALSE)
write.csv(EcologyMammalianTable, file="../../Body/2Derived/EcologyMammalianTable01.csv", quote = FALSE, row.names = FALSE)


#############kuptsov species

kuptsovtable = read.table("../../Body/2Derived/EcologyMammalianTable01_KuptsovA_ver2_Full.txt", sep='\t', header=TRUE)

allparameters$FrT= NULL
allparameters = merge(kuptsovtable, AGG, by="Species")
allparameters$Temper = as.numeric(gsub(",", ".", allparameters$Temperature.C._White2003.2006.other.close.species))
allparameters$GenerationLength_d = as.numeric(gsub(",", ".", allparameters$GenerationLength_d))

summary(allparameters$Temper)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#30.70   35.90   37.00   36.80   38.08   40.10     425 


#######Lm ~ hiber

ltest = lm(formula = FrT.y ~ scale(GenerationLength_d)+scale(Hib), data = allparameters)
summary(ltest)
#######Lm ~ Hib.unconfirmedHib

ltest = lm(formula = FrT.y ~ scale(GenerationLength_d)+scale(Hib.unconfirmedHib), data = allparameters)
summary(ltest)

#######Lm ~ Daily

ltest = lm(formula = FrT.y ~ scale(GenerationLength_d)+scale(Daily), data = allparameters)
summary(ltest)

#######Lm ~ Daily
ltest = lm(formula = FrT.y ~ scale(GenerationLength_d)+scale(Daily.unconfirmedDaily), data = allparameters)
summary(ltest)

#######Lm ~ MarsMono

allparameters$MarsMono = allparameters$Mars + allparameters$Mono
table(allparameters$MarsMono)

ltest = lm(formula = FrT.y ~ scale(GenerationLength_d)+scale(MarsMono), data = allparameters)
summary(ltest)

#######Lm ~ Temperature

ltest = lm(formula = FrT.y ~ scale(GenerationLength_d)+scale(Temper), data = allparameters)
summary(ltest)

#######Lm ~ cold

coldspeciesnames = allparameters[allparameters$Temper <= 35.90 & !is.na(allparameters$Temper),]$Species; length(coldspecies)
allparameters$colddummy = 0
allparameters[allparameters$Species %in% coldspeciesnames,]$colddummy = 1

ltest = lm(formula = FrT.y ~ scale(GenerationLength_d)+scale(colddummy), data = allparameters)
summary(ltest)

#######Lm ~ allcold

allparameters$allcolddummy = allparameters$Hib.unconfirmedHib + allparameters$Daily.unconfirmedDaily + allparameters$MarsMono + allparameters$colddummy
table(allparameters$allcolddummy)
allparameters[allparameters$allcolddummy > 0,]$allcolddummy = 1

ltest = lm(formula = FrT.y ~ log2(GenerationLength_d)+scale(allcolddummy), data = allparameters)
summary(ltest)


ltest = lm(formula = FrT.y ~ log2(GenerationLength_d), data = allparameters)
summary(ltest)
allparameters$residuals = ltest$residuals

ltest = lm(formula = FrT.y ~ scale(GenerationLength_d)+scale(allcolddummy), data = allparameters)
summary(ltest)

ltest = lm(formula = FrT.y ~ scale(GenerationLength_d)*scale(allcolddummy), data = allparameters)
summary(ltest)

#################Lm ~ Xen
allparameters$Xen = 0
allparameters[allparameters$Superorder == "Xenarthra", ]$Xen = 1
ltest = lm(formula = FrT.y ~ scale(GenerationLength_d)+scale(Xen), data = allparameters)
summary(ltest)

table(allparameters$Xen)

ltest = lm(formula = FrT.y ~ scale(GenerationLength_d)+scale(Xen)+scale(allcolddummy), data = allparameters)
summary(ltest)


##############################
####Fraction G
##############################

ltest = lm(formula = FrG ~ scale(GenerationLength_d), data = allparameters)
summary(ltest)

#######Lm ~ hiber

ltest = lm(formula = FrG ~ scale(GenerationLength_d)+scale(Hib), data = allparameters)
summary(ltest)

#########Lm ~ allcold

ltest = lm(formula = FrG ~ scale(allcolddummy), data = allparameters)
summary(ltest)

ltest = lm(formula = FrG ~ scale(GenerationLength_d)*scale(allcolddummy), data = allparameters)
summary(ltest)


##########mutspec
MUT = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt', header = TRUE)

HG = MUT[MUT$Species == "Heterocephalus_glaber",]
MM = MUT[MUT$Species == "Mus_musculus",]
CH = MUT[MUT$Species == "Cryptomys_hottentotus",]
	

HG$Temperature = allparameters[allparameters$Species == "Heterocephalus_glaber",]$Temper
MM$Temperature = allparameters[allparameters$Species == "Mus_musculus",]$Temper
HG$GL = allparameters[allparameters$Species == "Heterocephalus_glaber",]$GenerationLength_d
MM$GL = allparameters[allparameters$Species == "Mus_musculus",]$GenerationLength_d
CH$Temperature = allparameters[allparameters$Species == "Cryptomys_hottentotus",]$Temper
CH$GL = allparameters[allparameters$Species == "Cryptomys_hottentotus",]$GenerationLength_d

HGMM = rbind(HG, CH)

##############MUT + WG

MUTWG = merge(MUT, allparameters, by="Species")
ltest = lm(formula = T_C ~ scale(GenerationLength_d)+scale(allcolddummy), data = MUTWG)
summary(ltest)

ltest = lm(formula = T_C ~ scale(Xen)+scale(allcolddummy)+scale(GenerationLength_d), data = MUTWG)
summary(ltest)

ltest = lm(formula = T_C ~ scale(allcolddummy)+log2(GenerationLength_d), data = MUTWG)
summary(ltest)

ltest = lm(formula = A_G ~ scale(allcolddummy)+log2(GenerationLength_d), data = MUTWG)
summary(ltest)

#############residuals
summary(allparameters$residuals)
table(allparameters$Superorder)
table(allparameters$Order..or.infraorder.for.Cetacea.)

table(allparameters[allparameters$Superorder == "Laurasiatheria",]$Order..or.infraorder.for.Cetacea.)
table(allparameters[allparameters$Superorder == "Euarchontoglires",]$Order..or.infraorder.for.Cetacea.)

boxplot(allparameters[allparameters$Hib.unconfirmedHib == 1,]$residuals, allparameters[allparameters$Daily.unconfirmedDaily == 1,]$residuals, allparameters[allparameters$MarsMono == 1,]$residuals, allparameters[allparameters$colddummy == 1,]$residuals, allparameters[allparameters$allcolddummy == 1,]$residuals, allparameters[allparameters$allcolddummy == 0 & allparameters$Superorder == "Xenarthra",]$residuals, allparameters[allparameters$allcolddummy == 0 & allparameters$Superorder == "Euarchontoglires",]$residuals, allparameters[allparameters$allcolddummy == 0 & allparameters$Superorder == "Laurasiatheria",]$residuals, 
        allparameters[allparameters$allcolddummy == 0 & allparameters$Order..or.infraorder.for.Cetacea. == "Primates",]$residuals,
        allparameters[allparameters$allcolddummy == 0 & allparameters$Order..or.infraorder.for.Cetacea. == "Rodentia",]$residuals, 
        allparameters[allparameters$allcolddummy == 0 & allparameters$Order..or.infraorder.for.Cetacea. == "Artiodactyla",]$residuals,
        allparameters[allparameters$allcolddummy == 0 & allparameters$Order..or.infraorder.for.Cetacea. == "Carnivora",]$residuals,
        allparameters[allparameters$allcolddummy == 0 & allparameters$Order..or.infraorder.for.Cetacea. == "Chiroptera",]$residuals,
        allparameters[allparameters$allcolddummy == 0 & allparameters$Order..or.infraorder.for.Cetacea. == "Cetacea",]$residuals,
        allparameters[allparameters$allcolddummy == 0 & allparameters$Order..or.infraorder.for.Cetacea. == "Eulipotyphla",]$residuals,
        names = c("Hib", "Daily", "MarsMono", "Cold", "Allcold", "Xenarthra","Euarchontoglires", "Laurasiatheria", "Primates", "Rodentia", "Artiodactyla", "Carnivora", "Chiroptera", "Cetacea", "Eulipotyphla"), varwidth = TRUE, notch = TRUE)
abline(h = 0, col = "red")

summary(allparameters[allparameters$Order..or.infraorder.for.Cetacea. == "Cetacea",]$FrT.y)
summary(allparameters[allparameters$Order..or.infraorder.for.Cetacea. == "Artiodactyla",]$FrT.y)
