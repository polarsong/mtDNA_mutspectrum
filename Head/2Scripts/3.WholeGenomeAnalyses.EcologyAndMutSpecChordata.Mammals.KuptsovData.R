###################################
###### 
###################################
rm(list=ls(all=TRUE))
library(ggpubr)

############ Syn mut
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")) {file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")}

kuptsovtable = read.table("../../Body/2Derived/EcologyMammalianTable01_KuptsovA_ver2_Full.txt", sep='\t', header=TRUE)



  
SynNuc = SynNuc[SynNuc$Gene != 'ND6',]
SynNuc = aggregate(list(SynNuc$NeutralA,SynNuc$NeutralT,SynNuc$NeutralG,SynNuc$NeutralC), by = list(SynNuc$Species), FUN = sum)
names(SynNuc) = c('Species','NeutralA','NeutralT','NeutralG','NeutralC')
SynNuc$FrA = SynNuc$NeutralA / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC)
SynNuc$FrT = SynNuc$NeutralT / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 
SynNuc$FrG = SynNuc$NeutralG / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 
SynNuc$FrC = SynNuc$NeutralC / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 


kuptsovtable$FrT= NULL
allparameters = merge(kuptsovtable, SynNuc, by="Species")
allparameters$Temper = as.numeric(gsub(",", ".", allparameters$Temperature.C._White2003.2006.other.close.species))
allparameters$GenerationLength_d = as.numeric(gsub(",", ".", allparameters$GenerationLength_d))

allparameters$TG = allparameters$FrT+allparameters$FrG
allparameters$AC = allparameters$FrA+allparameters$FrC
allparameters$TG_ACSkew = (allparameters$TG-allparameters$AC)/(allparameters$TG+allparameters$AC); summary(allparameters$TG_ACSkew)
allparameters$TtoCSkew = (allparameters$FrT-allparameters$FrC)/(allparameters$FrT+allparameters$FrC); summary(allparameters$TtoCSkew)
 


summary(lm(TG_ACSkew ~ scale(Temper)+scale(GenerationLength_d), data = allparameters))
summary(lm(TG_ACSkew ~ log2(Temper)+log2(GenerationLength_d), data = allparameters))
summary(lm(TtoCSkew ~ scale(Temper)+scale(GenerationLength_d), data = allparameters))

summary(allparameters$Temper)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#30.70   35.90   37.00   36.80   38.08   40.10     425 


#######Lm ~ hiber

summary(lm(formula = FrT ~ scale(GenerationLength_d)+scale(Hib), data = allparameters))

#######Lm ~ Hib.unconfirmedHib

summary(lm(formula = FrT ~ scale(GenerationLength_d)+scale(Hib.unconfirmedHib), data = allparameters))

#######Lm ~ Daily

ltest = lm(formula = FrT ~ scale(GenerationLength_d)+scale(Daily), data = allparameters)
summary(ltest)

#######Lm ~ Daily
summary(lm(formula = FrT ~ scale(GenerationLength_d)+scale(Daily.unconfirmedDaily), data = allparameters))

#######Lm ~ MarsMono

allparameters$MarsMono = allparameters$Mars + allparameters$Mono
table(allparameters$MarsMono)

summary(lm(formula = FrT ~ scale(GenerationLength_d)+scale(MarsMono), data = allparameters))


#######Lm ~ Temperature

summary(lm(formula = FrT ~ scale(GenerationLength_d)+scale(Temper), data = allparameters))

#######Lm ~ cold
formediantemperature = allparameters[!is.na(allparameters$Temper),]$Temper
coldspeciesnames = allparameters[allparameters$Temper <= mean(formediantemperature) & !is.na(allparameters$Temper),]$Species
allparameters$colddummy = 0
allparameters[allparameters$Species %in% coldspeciesnames,]$colddummy = 1

summary(lm(formula = FrT ~ scale(GenerationLength_d)+scale(colddummy), data = allparameters))

#######Lm ~ allcold

allparameters$allcolddummy = allparameters$Hib.unconfirmedHib + allparameters$Daily.unconfirmedDaily + allparameters$MarsMono + allparameters$colddummy
table(allparameters$allcolddummy)
allparameters[allparameters$allcolddummy > 0,]$allcolddummy = 1

summary(lm(formula = FrT ~ log2(GenerationLength_d)+scale(allcolddummy), data = allparameters))
summary(lm(formula = FrT ~ scale(GenerationLength_d)+scale(allcolddummy), data = allparameters))

ltest = lm(formula = FrT  ~ log2(GenerationLength_d), data = allparameters)
summary(ltest)
allparameters$residuals = ltest$residuals ## residuals added

summary(lm(formula = TG_ACSkew ~ log2(GenerationLength_d)+allcolddummy, data = allparameters))
summary(lm(formula = TG_ACSkew ~ scale(GenerationLength_d)+scale(allcolddummy), data = allparameters))
summary(lm(formula = TG_ACSkew ~ GenerationLength_d+allcolddummy, data = allparameters))

nrow(allparameters[!is.na(allparameters$Temper),]) #224



#################Lm ~ Xen
allparameters$Xen = 0
allparameters[allparameters$Superorder == "Xenarthra", ]$Xen = 1
summary(lm(formula = FrT ~ scale(GenerationLength_d)+scale(Xen), data = allparameters))

summary(lm(formula = FrT ~ scale(GenerationLength_d)+scale(Xen)+scale(allcolddummy), data = allparameters))


##############################
####Fraction G
##############################

ltest = lm(formula = FrG ~ scale(GenerationLength_d), data = allparameters)
summary(ltest)

#######Lm ~ hiber

ltest = lm(formula = FrG ~ scale(GenerationLength_d)+scale(Hib.unconfirmedHib), data = allparameters)
summary(ltest)

#########Lm ~ allcold

ltest = lm(formula = FrG ~ scale(allcolddummy), data = allparameters)
summary(ltest)

ltest = lm(formula = FrG ~ scale(GenerationLength_d)*scale(allcolddummy), data = allparameters)
summary(ltest)


allparameters$TwoMammaliaGroups = as.character(allparameters$allcolddummy)
for (i in 1:nrow(allparameters)){
  if(allparameters$TwoMammaliaGroups[i] == 0){
    allparameters$TwoMammaliaGroups[i] = "Warmer mammals"
  }else{
    allparameters$TwoMammaliaGroups[i] = "Colder mammals"
  }
}


pdf("../../Body/4Figures/WholeGenomeAnalyses.EcologyAndMutSpecChordata.Mammals.KuptsovData.TemperatureHist.pdf", width = 9, height = 5.5)
gghistogram(allparameters, x = "Temper", fill = "red", ylab = " ", xlab = "Body temperature, °C", 
            add = "mean", rug = TRUE)
dev.off()


pdf("../../Body/4Figures/WholeGenomeAnalyses.EcologyAndMutSpecChordata.Mammals.KuptsovData.FIGURE2D.pdf", width = 9, height = 5.5)
ggscatter(allparameters, x = "GenerationLength_d", y = "FrT",
          color = "TwoMammaliaGroups", shape = "TwoMammaliaGroups",
          palette = c("#08519c", "#de6a85"),
          ellipse = TRUE,  xscale = "log2", xlab="Generation Length, log2", ylab="Fraction of A")
dev.off()

pdf("../../Body/4Figures/WholeGenomeAnalyses.EcologyAndMutSpecChordata.Mammals.KuptsovData.FIGURE2E.pdf", width = 9, height = 5.5)
ggscatter(allparameters, x = "GenerationLength_d", y = "TG_ACSkew",
          color = "TwoMammaliaGroups", shape = "TwoMammaliaGroups",
          palette = c("#08519c", "#de6a85"),
          ellipse = TRUE,  add = "reg.line", xscale = "log2", xlab="Generation Length, log2", ylab="Fraction of AC_TCSkew", )
dev.off()

colfunc<-colorRampPalette(c("red","yellow","springgreen","royalblue"))
plot(allparameters$Temper, allparameters$GenerationLength_d, col="#42cbf5", xlab="Temperature", ylab="GL, days", pch = 1, cex = (allparameters$TG_ACSkew*-7), ylim = c(0, 11000))


medGL = median(allparameters$GenerationLength_d)
allparameters$GLgroups = "LongGL"
allparameters[allparameters$GenerationLength_d < medGL,]$GLgroups = "ShortGL"
allparameters$DummyGL = 0
allparameters[allparameters$GenerationLength_d < medGL,]$DummyGL = 1


pdf("../../Body/4Figures/WholeGenomeAnalyses.EcologyAndMutSpecChordata.FIGURE.AllVertebrates.pdf", width = 10, height = 11)
plot(allparameters[allparameters$GLgroups == "LongGL",]$Temper, allparameters[allparameters$GLgroups == "LongGL",]$TG_ACSkew*-1, col="#121d2e", xlab="Temperature", ylab="STG-SAC", xlim = c(0, 41), ylim = c (0.1, 0.75),  pch = 16)
abline(0.362350, 0.004390, col="black", lwd = 2)
par(new=TRUE)
plot(allparameters[allparameters$GLgroups == "ShortGL",]$Temper, allparameters[allparameters$GLgroups == "ShortGL",]$TG_ACSkew*-1, col="#d7badb", xlab=" ", ylab="", xlim = c(0, 41), ylim = c (0.1, 0.75),  pch = 16)
abline((0.362350-0.043854), 0.004390, col="#d7badb", lwd = 2)
par(new=TRUE)
plot(temp[temp$Lifespan == "ShortMaturated",]$Temperature, temp[temp$Lifespan == "ShortMaturated",]$TG_ACSkew*-1, col="#4da36c", xlab=" ", ylab="  ", xlim = c(0, 41), ylim = c(0.1, 0.75),  pch = 16)
abline((0.272616+0.006709 ), 0.006184 , col="#4da36c", lwd = 2)
par(new=TRUE)
plot(temp[temp$Lifespan == "LongMaturated",]$Temperature, temp[temp$Lifespan == "LongMaturated",]$TG_ACSkew*-1, col="#42cbf5", xlab=" ", ylab= "", xlim = c(0, 41), ylim = c(0.1, 0.75),  pch = 16)
abline(0.272616, 0.006184 , col="#42cbf5", lwd = 2)
dev.off()

summary(lm(formula = -TG_ACSkew ~ Temper+DummyGL, data = allparameters))


pdf("../../Body/4Figures/WholeGenomeAnalyses.EcologyAndMutSpecChordata.Mammals.KuptsovData.FIGURE2B.withoutEllipse.pdf", width = 9, height = 5.5)
plot(log2(allparameters[allparameters$TwoMammaliaGroups == "Colder mammals",]$GenerationLength_d), allparameters[allparameters$TwoMammaliaGroups == "Colder mammals",]$TG_ACSkew, col="#08519c", xlab="Generation Length, log2", ylab="Fraction of AC_TCSkew", xlim = c(8, 14.5), ylim = c(-0.75, -0.25))
abline((-0.323919+0.043515), -0.018255, col="#08519c", lwd = 2)
par(new=TRUE)
plot(log2(allparameters[allparameters$TwoMammaliaGroups == "Warmer mammals",]$GenerationLength_d), allparameters[allparameters$TwoMammaliaGroups == "Warmer mammals",]$TG_ACSkew, col="#de6a85", xlab="Generation Length, log2", ylab="Fraction of AC_TCSkew", xlim = c(8, 14.5), ylim = c(-0.75, -0.25))
abline(-0.323919, -0.018255, col="#de6a85", lwd = 2)
dev.off()

breaks = seq(-0.75, -0.25, 0.01)
#par(mfrow = c(2,1))
hist(allparameters[allparameters$TwoMammaliaGroups == "Warmer mammals",]$TG_ACSkew, col="#de6a85", breaks = breaks, xlim = c(-0.75, -0.25), ylim = c(0, 30))
par(new=TRUE)
hist(allparameters[allparameters$TwoMammaliaGroups == "Colder mammals",]$TG_ACSkew, col="#08519c", breaks = breaks, xlim = c(-0.75, -0.25), ylim = c(0, 30))

pdf("../../Body/4Figures/WholeGenomeAnalyses.EcologyAndMutSpecChordata.Mammals.KuptsovData.FIGURE2C.pdf", width = 9, height = 5.5)
plot(density(allparameters[allparameters$TwoMammaliaGroups == "Warmer mammals",]$TG_ACSkew), col="#de6a85", xlim = c(-0.75, -0.25), ylim = c(0, 6))
par(new=TRUE)
plot(density(allparameters[allparameters$TwoMammaliaGroups == "Colder mammals",]$TG_ACSkew), col="#08519c", xlim = c(-0.75, -0.25), ylim = c(0, 6))
dev.off()
##########mutspec
MUT = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt', header = TRUE)

HG = MUT[MUT$Species == "Heterocephalus_glaber",]
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



#################################metabolic rate approximation

allparameters$TtoCSkew = (allparameters$FrT-allparameters$FrC)/(allparameters$FrT+allparameters$FrC); summary(allparameters$TtoCSkew)

allparameters$MR=(allparameters$GenerationLength_d+1)^0.75
allparameters$TemperatureK = 273.15 + allparameters$Temper
allparameters$B=allparameters$MR * exp(-1.2/((8.617*10^-5)*allparameters$TemperatureK))
cor.test(allparameters$B, allparameters$TtoCSkew, method="spearman") #rho  -0.4568971 
cor.test(allparameters$Temper, allparameters$TtoCSkew, method = "spearman")

allparameters$TemperatureK = 273.15 + allparameters$Temper
form=allparameters[!is.na(allparameters$GenerationLength_d),]
form=form[!is.na(form$TG_ACSkew),]
form=form[!is.na(form$TemperatureK),]
form$TG_ACSkew=form$TG_ACSkew * (-1)

form$x = (1/log(form$GenerationLength_d))*(log(form$TG_ACSkew)+0.2/((8.617*10^-5)*form$TemperatureK))
bodymassgenlpred=((form$GenerationLength_d/365)/(6.10*10^6))^(1/0.2)
form$x = log(form$TG_ACSkew * exp(0.62/((8.617*10^-5)*form$TemperatureK))) / log(bodymassgenlpred)
summary(form$x)
