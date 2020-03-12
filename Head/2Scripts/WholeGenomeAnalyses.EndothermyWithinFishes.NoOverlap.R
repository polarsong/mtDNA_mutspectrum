###################################
###### 
###################################



rm(list=ls(all=TRUE))

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


########## compare nucleotide frequencies 
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


###########create fishes species table with Temp
Eu = read.table("../../Body/1Raw/EdothermicFishes.txt", sep = '\t')
ListOfEuSpecies = gsub(' ','_',Eu$V1); length(ListOfEuSpecies)
TP=read.table("../../Body/1Raw/FishBaseTemperature.txt", header = TRUE)
TP = aggregate(Temperature ~ ., mean, data = TP)
TP$EuFish = 1
TP[!TP$Species %in% ListOfEuSpecies,]$EuFish = 0
allparameters = merge(TP, AGG)
table(allparameters$EuFish)
allparameters = merge(allparameters, AA, all.x = TRUE)
EcologyFishTable= data.frame(allparameters$Kingdom, allparameters$Phylum, allparameters$Class, allparameters$Order, allparameters$Family,allparameters$Genus, allparameters$Species, allparameters$Temperature, allparameters$FrT, allparameters$EuFish)
write.table(EcologyFishTable, file="../../Body/2Derived/EcologyFishTable.txt", quote = FALSE, row.names = FALSE)

ltest = lm(formula = FrT ~ scale(Temperature)*scale(EuFish), data = allparameters)
summary(ltest)
ltest = lm(formula = FrT ~ scale(Temperature)+scale(EuFish), data = allparameters)
summary(ltest)
ltest = lm(formula = scale(FrT) ~ 0 + scale(Temperature)+scale(EuFish), data = allparameters)
summary(ltest)

ggscatter(allparameters, x = "Temperature", y = "FrT",
          palette = c("#00AFBB"), add = "reg.line", xlab="Temperature, C")


ggscatter(allparameters, x = "Temperature", y = "FrT",
          color = "blue", # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x = 3, label.sep = "\n")
)
