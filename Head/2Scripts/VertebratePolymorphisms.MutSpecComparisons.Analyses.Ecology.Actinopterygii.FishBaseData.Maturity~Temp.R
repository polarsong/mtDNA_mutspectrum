rm(list=ls(all=TRUE))


###TEMPERATURE another table
MUT = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt', header = TRUE)
class(MUT$A_T)
class(MUT$Species)

TEMPE = read.table('../../Body/1Raw/FishBaseTemperature2.txt')                         # tempetature table
colnames(TEMPE)=c("Species", "Temperature")
TEMPE$Species=gsub(" ", "_", TEMPE$Species)
class(TEMPE$Temperature)
class(TEMPE$Species)


### calculating average for repeating species in TEMPE
averageTEMPE = aggregate(Temperature ~ ., mean, data = TEMPE)
averageTemperMut = merge(MUT,averageTEMPE) 


######MATURITY
MATU = read.table('../../Body/1Raw/FishBaseLmTm.txt')
colnames(MATU)=c("Species", "Lm", "Tm")
MATU$Species=gsub(" ", "_", MATU$Species)
MATU=MATU[!is.na(MATU$Tm),]
MATU=MATU[MATU$Lm != "",]                                                            #maturity table


##multiple reg 

allparameters=merge(averageTemperMut, MATU)


test= lm(T_C ~ Temperature + Tm, data=allparameters)

plot(test)
