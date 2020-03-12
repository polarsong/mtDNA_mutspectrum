rm(list=ls(all=TRUE))
library(ggpubr)


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


ggscatter(TemperMut, x = "Temperature", y = "T_C",
          color = "blue", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x = 3, label.sep = "\n"), xlab="Temperature, C", ylab="A>G")

#####checking A_G


ltest = lm(formula = A_G ~ scale(Temperature)*scale(Tm), data = allparameters)
summary(ltest)

ltest = lm(formula = Temperature ~ scale(T_C)*scale(A_G), data = allparameters)
summary(ltest)

ltest = lm(formula = Temperature ~ scale(T_C)+scale(A_G), data = allparameters)
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





