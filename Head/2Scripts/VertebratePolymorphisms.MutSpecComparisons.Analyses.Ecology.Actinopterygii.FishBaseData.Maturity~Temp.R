rm(list=ls(all=TRUE))




###TEMPERATURE another table
MUT = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt', header = TRUE)
class(MUT$A_T)
class(MUT$Species)

TEMPE = read.table('../../Body/1Raw/FishBaseTemperature.txt', header = TRUE)                         # tempetature table
class(TEMPE$Temperature)
class(TEMPE$Species)


### calculating average for repeating species in TEMPE
averageTEMPE = aggregate(Temperature ~ ., mean, data = TEMPE)
averageTemperMut = merge(MUT,averageTEMPE) 


######MATURITY
MATU = read.table('../../Body/1Raw/FishBaseLmTm.txt')
colnames(MATU)=c("Species", "Lm", "Tm")
MATU$Species=gsub(" ", "_", MATU$Species)
MATU=MATU[!is.na(MATU$Tm),]                           #maturity table
MATU$Lm = NULL
MATU = aggregate(Tm ~ ., mean, data = MATU) 



##multiple reg 

allparameters=merge(averageTemperMut, MATU)


ltest = lm(formula = T_C ~ scale(Temperature)*scale(Tm), data = allparameters)
summary(ltest)


ltest = lm(formula = T_C ~ scale(Temperature) + scale(Tm), data = allparameters)
summary(ltest)


