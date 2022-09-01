rm(list=ls(all=TRUE))

if (!require(caper)) install.packages("caper")
if (!require(geiger)) install.packages("geiger")
if (!require(ggpubr)) install.packages("ggpubr")
library(caper)
library(geiger)
library("ggpubr")

### reading whole genomes database
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")) {file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")}

######reading ecology table from us + Kuptsov 
kuptsovtable = read.table("../../Body/2Derived/EcologyMammalianTable01_KuptsovA_ver2_Full.txt", sep='\t', header=TRUE)

####### obtaining neutral nucleotide fractions in whole genomes
SynNuc = SynNuc[SynNuc$Gene != 'ND6',]
SynNuc = aggregate(list(SynNuc$NeutralA,SynNuc$NeutralT,SynNuc$NeutralG,SynNuc$NeutralC), by = list(SynNuc$Species), FUN = sum)
names(SynNuc) = c('Species','NeutralA','NeutralT','NeutralG','NeutralC')
SynNuc$FrA = SynNuc$NeutralA / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC)
SynNuc$FrT = SynNuc$NeutralT / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 
SynNuc$FrG = SynNuc$NeutralG / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 
SynNuc$FrC = SynNuc$NeutralC / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 

##### merging whole genomes with ecology table
kuptsovtable$FrT= NULL
allparameters = merge(kuptsovtable, SynNuc, by="Species")
allparameters$Temper = as.numeric(gsub(",", ".", allparameters$Temperature.C._White2003.2006.other.close.species))
allparameters$GenerationLength_d = as.numeric(gsub(",", ".", allparameters$GenerationLength_d))
allparameters$TG = allparameters$FrT+allparameters$FrG
allparameters$AC = allparameters$FrA+allparameters$FrC
allparameters$TG_ACSkew = (allparameters$TG-allparameters$AC)/(allparameters$TG+allparameters$AC); summary(allparameters$TG_ACSkew)
allparameters$TtoCSkew = (allparameters$FrT-allparameters$FrC)/(allparameters$FrT+allparameters$FrC); summary(allparameters$TtoCSkew)
allparameters$AC_TGSkew = allparameters$TG_ACSkew *-1 
#allparameters$TCskew = (allparameters$FrT - allparameters$FrC)/(allparameters$FrT + allparameters$FrC)
#allparameters$CTskew = (allparameters$FrC - allparameters$FrT)/(allparameters$FrT + allparameters$FrC)
summary(allparameters$Temper)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#30.70   35.90   37.00   36.80   38.08   40.10     425 
nrow(allparameters[!is.na(allparameters$Temper) & !is.na(allparameters$GenerationLength_d),]) #224


###### Multiple models 
summary(lm(TG_ACSkew ~ log2(Temper)+log2(GenerationLength_d), data = allparameters))
summary(lm(TtoCSkew ~ scale(Temper)+scale(GenerationLength_d), data = allparameters))
summary(lm(TG_ACSkew ~ scale(Temper)+scale(GenerationLength_d), data = allparameters))
#Supl mat. 4 a
summary(lm(AC_TGSkew ~ scale(Temper)+scale(GenerationLength_d), data = allparameters)) ###PICS
#Supl mat. Fig 4.2
summary(lm(TCskew ~ scale(Temper)+scale(GenerationLength_d), data = allparameters))
summary(lm(CTskew ~ scale(Temper)+scale(GenerationLength_d), data = allparameters))

########################
###ANALYSES WITH DUMMY VARIABLES
########################

#######Lm ~ hibernators

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



#Supl mat. 4b
summary(lm(formula = AC_TGSkew ~ log2(GenerationLength_d)+scale(allcolddummy), data = allparameters))

summary(lm(formula = FrT ~ log2(GenerationLength_d)+scale(allcolddummy), data = allparameters))
summary(lm(formula = FrT ~ scale(GenerationLength_d)+scale(allcolddummy), data = allparameters))

ltest = lm(formula = FrT  ~ log2(GenerationLength_d), data = allparameters)
summary(ltest)
allparameters$residuals = ltest$residuals ## residuals added




allparameters=allparameters[allparameters$Order..or.infraorder.for.Cetacea. != "Chiroptera",]
allparameters=allparameters[allparameters$MarsMono != 1,]

##### phylogenetic inertia analysis
#allparameters = allparameters[!is.na(allparameters$Temper) & !is.na(allparameters$GenerationLength_d),]
tree = read.tree('../../Body/1Raw/mtalign.aln.treefile.rooted')


#ForG = data.frame(allparameters$Species, allparameters$Temper, allparameters$AC_TGSkew, allparameters$GenerationLength_d, allparameters$allcolddummy)
#names(ForG) = c("Species", "TemperatureC", "SAG_STC", "LongevityDAYS", "DummyCOLDER")
#write.table(ForG, file="../../Body/2Derived/ForKGMammaliaDataSet.txt", row.names = F)
row.names(allparameters) = allparameters$Species

tree_pruned = treedata(tree, allparameters, sort=T, warnings=T)$phy 

data<-as.data.frame(treedata(tree_pruned, allparameters, sort=T, warnings=T)$data)
data$Species = as.character(data$Species)

data$AC_TGSkew = as.numeric(as.character(data$AC_TGSkew))
data$Temper = as.numeric(as.character(data$Temper))
data$GenerationLength_d = as.numeric(as.character(data$GenerationLength_d))
data$allcolddummy = as.numeric(as.character(data$allcolddummy))

data_comp <- comparative.data(tree_pruned, data[, c('Species', 'AC_TGSkew',
                                                    'GenerationLength_d', "Temper")], Species, vcv=TRUE)

model = pgls(AC_TGSkew ~ scale(Temper) + scale(GenerationLength_d), data_comp, lambda="ML")
summary(model)
model = pgls(CTskew ~ scale(Temper) + scale(GenerationLength_d), data_comp, lambda="ML")
summary(model)


data_comp <- comparative.data(tree_pruned, data[, c('Species', 'AC_TGSkew',
                                                    'GenerationLength_d', "allcolddummy")], Species, vcv=TRUE)

nrow(data[!is.na(data$allcolddummy) & !is.na(data$GenerationLength_d),])

model = pgls(AC_TGSkew ~ allcolddummy + scale(GenerationLength_d), data_comp, lambda="ML")
summary(model)
model = pgls(AC_TGSkew ~ allcolddummy + GenerationLength_d, data_comp, lambda="ML")
summary(model)
model = pgls(AC_TGSkew ~ allcolddummy + log2(GenerationLength_d), data_comp, lambda="ML")
summary(model)
model = pgls(AC_TGSkew ~ allcolddummy + log2(GenerationLength_d), data_comp, lambda="ML")
summary(model)

model = pgls(AC_TGSkew ~ scale(Temper) + scale(GenerationLength_d), data_comp, lambda="ML")
summary(model)

# lambda [ ML]  : 1.000
# Coefficients:
#   Estimate Std. Error t value  Pr(>|t|)    
# (Intercept)                0.4455734  0.0508286  8.7662 4.441e-16 ***
#   scale(Temper)             -0.0048257  0.0054858 -0.8797     0.380    
# scale(GenerationLength_d)  0.0036460  0.0042781  0.8523     0.395    

model2 = pgls(AC_TGSkew ~ log2(Temper + 2) + log2(GenerationLength_d), data_comp, lambda="ML")
summary(model2)

# lambda [ ML]  : 1.000
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)               0.7338936  0.4312905  1.7016  0.09023 .
# log2(Temper + 2)         -0.0688118  0.0808772 -0.8508  0.39579  
# log2(GenerationLength_d)  0.0068201  0.0048876  1.3954  0.16430  

model3 = pgls(AC_TGSkew ~ log2(GenerationLength_d), data_comp, lambda="ML")
summary(model3)

# lambda [ ML]  : 1.000
# Coefficients:
#   Estimate Std. Error t value  Pr(>|t|)    
# (Intercept)        4.4578e-01 5.0767e-02  8.7809 4.441e-16 ***
#   GenerationLength_d 2.3930e-06 2.4945e-06  0.9593    0.3384  

# as lambda is 1.0, we can use PIC

cor.test(pic(data$AC_TGSkew, tree_pruned), pic(data$GenerationLength_d, tree_pruned))
# 0.08146519, p-value 0.03815
cor.test(pic(data$AC_TGSkew, tree_pruned), pic(data$allcolddummy, tree_pruned))

summary(lm(pic(data$AC_TGSkew, tree_pruned) ~ pic(data$Temper, tree_pruned) + pic(data$GenerationLength_d, tree_pruned)))
summary(lm(pic(data$AC_TGSkew, tree_pruned) ~ pic(data$allcolddummy, tree_pruned) + pic(data$GenerationLength_d, tree_pruned)))

model4 = pgls(AC_TGSkew ~ GenerationLength_d, data_comp, lambda="ML")
summary(model4)

# lambda [ ML]  : 1.000

# Coefficients:
#   Estimate Std. Error t value  Pr(>|t|)    
# (Intercept)              0.3692538  0.0735997  5.0171 1.078e-06 ***
#   log2(GenerationLength_d) 0.0069482  0.0048809  1.4235    0.1560    
# allcolddummy             0.0071838  0.0098121  0.7321    0.4649

##### after calculating temp dummy

# data_comp_dummy <- comparative.data(tree_pruned, data[, c('Species', 'AC_TGSkew',
#                                                     'GenerationLength_d', 'Temper',
#                                                     'allcolddummy')], Species, vcv=TRUE)
# 
# model2 = pgls(AC_TGSkew ~ allcolddummy + scale(GenerationLength_d), data_comp_dummy, lambda="ML")
# summary(model2)

# Coefficients:
#   Estimate Std. Error t value  Pr(>|t|)    
# (Intercept)               0.4451936  0.0511134  8.7099 6.661e-16 ***
#   allcolddummy              0.0074469  0.0098631  0.7550    0.4510    
# scale(GenerationLength_d) 0.0037206  0.0042780  0.8697    0.3854 









##################Figures
allparameters$TwoMammaliaGroups = as.character(allparameters$allcolddummy)
for (i in 1:nrow(allparameters)){
  if(allparameters$TwoMammaliaGroups[i] == 0){
    allparameters$TwoMammaliaGroups[i] = "Warmer mammals"
  }else{
    allparameters$TwoMammaliaGroups[i] = "Colder mammals"
  }
}




pdf("../../Body/4Figures/WholeGenomeAnalyses.EcologyAndMutSpecChordata.FIGURE.AllVertebrates.pdf", width = 7, height = 8.5)
plot(allparameters[allparameters$GLgroups == "LongGL",]$Temper, allparameters[allparameters$GLgroups == "LongGL",]$TG_ACSkew*-1, col="#121d2e", xlab="Temperature", ylab="STG-SAC", xlim = c(0, 41), ylim = c (0.1, 0.75),  pch = 16)
abline(0.362350, 0.004390, col="black", lwd = 2)
par(new=TRUE)
plot(allparameters[allparameters$GLgroups == "ShortGL",]$Temper, allparameters[allparameters$GLgroups == "ShortGL",]$TG_ACSkew*-1, col="#d7badb", xlab=" ", ylab="", xlim = c(0, 41), ylim = c (0.1, 0.75),  pch = 16)
abline((0.362350-0.043854), 0.004390, col="#d7badb", lwd = 2)
par(new=TRUE)
plot(temp[temp$Lifespan == "ShortMaturated",]$Temperature, temp[temp$Lifespan == "ShortMaturated",]$TG_ACSkew*-1, col="#4da36c", xlab=" ", ylab="  ", xlim = c(0, 41), ylim = c(0.1, 0.75),  pch = 16)
abline((0.331911-0.049196), 0.006172, col="#4da36c", lwd = 2)
par(new=TRUE)
plot(temp[temp$Lifespan == "LongMaturated",]$Temperature, temp[temp$Lifespan == "LongMaturated",]$TG_ACSkew*-1, col="#42cbf5", xlab=" ", ylab= "", xlim = c(0, 41), ylim = c(0.1, 0.75),  pch = 16)
abline(0.331911, 0.006172, col="#42cbf5", lwd = 2)
legend("bottomright", legend=c("Mammals, Long GL", "Mammals, Short GL", "Fishes, Long TM", "Fishes, Short TM"), col=c("black", "#d7badb", "#42cbf5", "#4da36c"), pch = 16)
dev.off()




pdf("../../Body/4Figures/WholeGenomeAnalyses.EcologyAndMutSpecChordata.Mammals.KuptsovData.TemperatureHist.pdf", width = 9, height = 5.5)
gghistogram(allparameters, x = "Temper", fill = "red", ylab = " ", xlab = "Body temperature, ?C", 
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
          ellipse = TRUE,  add = "reg.line", xscale = "log2", xlab="Generation Length, log2", ylab="Fraction of AC_TCSkew")
dev.off()

plot(allparameters$Temper, allparameters$GenerationLength_d, col="#42cbf5", xlab="Temperature", ylab="GL, days", pch = 1, cex = (allparameters$TG_ACSkew*-7), ylim = c(0, 11000))



pdf("../../Body/4Figures/WholeGenomeAnalyses.EcologyAndMutSpecChordata.FIGURE.AllVertebrates.pdf", width = 7, height = 8.5)
plot(allparameters[allparameters$GLgroups == "LongGL",]$Temper, allparameters[allparameters$GLgroups == "LongGL",]$TG_ACSkew*-1, col="#121d2e", xlab="Temperature", ylab="STG-SAC", xlim = c(0, 41), ylim = c (0.1, 0.75),  pch = 16)
abline(0.362350, 0.004390, col="black", lwd = 2)
par(new=TRUE)
plot(allparameters[allparameters$GLgroups == "ShortGL",]$Temper, allparameters[allparameters$GLgroups == "ShortGL",]$TG_ACSkew*-1, col="#d7badb", xlab=" ", ylab="", xlim = c(0, 41), ylim = c (0.1, 0.75),  pch = 16)
abline((0.362350-0.043854), 0.004390, col="#d7badb", lwd = 2)
par(new=TRUE)
plot(temp[temp$Lifespan == "ShortMaturated",]$Temperature, temp[temp$Lifespan == "ShortMaturated",]$TG_ACSkew*-1, col="#4da36c", xlab=" ", ylab="  ", xlim = c(0, 41), ylim = c(0.1, 0.75),  pch = 16)
abline((0.331911-0.049196), 0.006172, col="#4da36c", lwd = 2)
par(new=TRUE)
plot(temp[temp$Lifespan == "LongMaturated",]$Temperature, temp[temp$Lifespan == "LongMaturated",]$TG_ACSkew*-1, col="#42cbf5", xlab=" ", ylab= "", xlim = c(0, 41), ylim = c(0.1, 0.75),  pch = 16)
abline(0.331911, 0.006172, col="#42cbf5", lwd = 2)
legend("bottomright", legend=c("Mammals, Long GL", "Mammals, Short GL", "Fishes, Long TM", "Fishes, Short TM"), col=c("black", "#d7badb", "#42cbf5", "#4da36c"), pch = 16)
dev.off()


pdf("../../Body/4Figures/WholeGenomeAnalyses.EcologyAndMutSpecChordata.Mammals.KuptsovData.FIGURE2B.withoutEllipse.pdf", width = 7, height = 8.5)
plot(log2(allparameters[allparameters$TwoMammaliaGroups == "Colder mammals",]$GenerationLength_d), allparameters[allparameters$TwoMammaliaGroups == "Colder mammals",]$AC_TGSkew, col="#08519c", xlab="Generation Length, log2", ylab="STG-SAC skew",  pch = 16, ylim=c(0.3,0.73), xlim=c(8, 14.5))
abline((0.323919-0.043515), 0.018255, col="#08519c", lwd = 2)
par(new=TRUE)
plot(log2(allparameters[allparameters$TwoMammaliaGroups == "Warmer mammals",]$GenerationLength_d), allparameters[allparameters$TwoMammaliaGroups == "Warmer mammals",]$AC_TGSkew, col="#de6a85", xlab="", ylab="",  pch = 16, ylim=c(0.3,0.73), xlim=c(8, 14.5))
abline(0.323919, 0.018255, col="#de6a85", lwd = 2)
legend("bottomright", legend=c("Colder mammals", "Warmer mammals"), col=c("#08519c", "#de6a85"), pch = 16)
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





################### one model

onedatasetF = data.frame(temp$Species, temp$Temperature, temp$Lifespan, temp$AC_TGSkew)
onedatasetF$Longevity = 0
onedatasetF[onedatasetF$temp.Lifespan == "LongMaturated",]$Longevity=1
names(onedatasetF)=c("Species", "Temperature", "Lifespan", "AC_TGSkew", "Longevity")
onedatasetF$Ectothermy = 1

onedatasetM = data.frame(allparameters$Species, allparameters$Temper, allparameters$GLgroups, allparameters$AC_TGSkew)
onedatasetM$Longevity = 0
onedatasetM[onedatasetM$allparameters.GLgroups == "LongGL",]$Longevity = 1
names(onedatasetM)=c("Species", "Temperature", "Lifespan", "AC_TGSkew", "Longevity")
onedatasetM$Ectothermy = 0

onedata = rbind(onedatasetF, onedatasetM)
table(onedata$Ectothermy)
summary(lm(formula = AC_TGSkew ~ Temperature + Longevity + Ectothermy, data = onedata))
nrow(onedata[!is.na(onedata$Temperature) & !is.na(onedata$AC_TGSkew),])


#######Mammals with GL
GT = read.table('../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt', sep = '\t', header = TRUE)
GT$Species = gsub(' ','_',GT$Scientific_name)
ALL = merge(SynNuc, GT)
ALL$TG = ALL$FrT+ALL$FrG
ALL$AC = ALL$FrA+ALL$FrC

ALL$AC_TG = (ALL$AC-ALL$TG)/(ALL$TG+ALL$AC); summary(ALL$AC_TG)
cor.test(ALL$AC_TG, ALL$GenerationLength_d, method = "spearman")
ALL$CTSkew = (ALL$FrC - ALL$FrT)/(ALL$FrC + ALL$FrT)
summary(ALL$CTSkew)
cor.test(ALL$CTSkew, ALL$GenerationLength_d, method = "spearman")

ALL$GASkew = (ALL$FrG - ALL$FrA)/(ALL$FrG + ALL$FrA)
summary(ALL$GASkew)
cor.test(ALL$GASkew, ALL$GenerationLength_d, method = "spearman")





