rm(list=ls(all=TRUE))



MUT = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt', header = TRUE)
class(MUT$A_T)
class(MUT$Species)

######Temperature
TEMPE = read.table('../../Body/1Raw/FishBaseTemperature.txt', header = TRUE)                        
class(TEMPE$Temperature)
class(TEMPE$Species)
averageTEMPE = aggregate(Temperature ~ ., median, data = TEMPE)




######MATURITY
MATULM = read.table('../../Body/1Raw/FishBaseMaturity_Lm.txt',  header = TRUE, stringsAsFactors=FALSE)
MATUTM = read.table('../../Body/1Raw/FishBaseMaturity_Tm.txt',  header = TRUE)
class(MATULM$Lm)
MATULM$Lm = as.numeric(MATULM$Lm)
MATUTM = aggregate(Tm ~ ., median, data = MATUTM) 
MATULM = aggregate(Lm ~ ., median, data = MATULM)




##multiple reg Tm
TemperMut = merge(MUT,averageTEMPE) 
allparameters=merge(TemperMut, MATUTM)


ltest = lm(formula = T_C ~ scale(Temperature)*scale(Tm), data = allparameters)
summary(ltest)


ltest = lm(formula = T_C ~ scale(Temperature) + scale(Tm), data = allparameters)
summary(ltest)

ltest = lm(formula = scale(T_C) ~ scale(Temperature) + scale(Tm), data = allparameters)
summary(ltest)

ltest = lm(formula = scale(T_C) ~ 0 + scale(Temperature) + scale(Tm), data = allparameters)
summary(ltest)

ltest = lm(formula = scale(T_C) ~ 0 + scale(Temperature), data = allparameters)
summary(ltest)


##multiple reg Lm
TemperMut = merge(MUT,averageTEMPE) 
allparameters=merge(TemperMut, MATULM)


ltest = lm(formula = T_C ~ scale(Temperature)*scale(Lm), data = allparameters)
summary(ltest)


ltest = lm(formula = T_C ~ scale(Temperature) + scale(Lm), data = allparameters)
summary(ltest)

ltest = lm(formula = scale(T_C) ~ scale(Temperature) + scale(Lm), data = allparameters)
summary(ltest)

ltest = lm(formula = scale(T_C) ~ 0 + scale(Temperature) + scale(Lm), data = allparameters)
summary(ltest)

ltest = lm(formula = scale(T_C) ~ 0 + scale(Temperature), data = allparameters)
summary(ltest)

