rm(list=ls(all=TRUE))


###TEMPERATURE another table
MUT = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt', header = TRUE)
class(MUT$A_T)
class(MUT$Species)

TEMPE = read.table('../../Body/1Raw/FishBaseTemperature2.txt')
colnames(TEMPE)=c("Species", "Temperature")
TEMPE$Species=gsub(" ", "_", TEMPE$Species)
class(TEMPE$Temperature)
class(TEMPE$Species)


### calculating average for repeating species in TEMPE
averageTEMPE = aggregate(Temperature ~ ., mean, data = TEMPE)
averageTemperMut = merge(MUT,averageTEMPE) 

library("ggpubr")

ggscatter(averageTemperMut, x = "Temperature", y = "T_C", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Temperature", ylab = "A>G")

######MATURITY
MATU = read.table('../../Body/1Raw/FishBaseLmTm.txt')
colnames(MATU)=c("Species", "Lm", "Tm")
MATU$Species=gsub(" ", "_", MATU$Species)
MATUTm=MATU[!is.na(MATU$Tm),]
MATULm=MATU[MATU$Lm != "",]
MATULmmut = merge(MUT, MATULm) 
MATUTmmut = merge(MUT, MATUTm) 

for (i in 1:nrow(MATULmmut)){
  MATULmmut$Lm_value[i]=as.numeric(unlist(strsplit(as.character(MATULmmut$Lm[i]), ' '))[1])
  MATULmmut$Lm_measure[i]=unlist(strsplit(as.character(MATULmmut$Lm[i]), ' '))[2]
}




###### Maturity Tm ~ T_C 
cor.test(MATUTmmut$T_C,MATUTmmut$Tm, method = 'spearman') 
ggscatter(MATUTmmut, x = "Tm", y = "T_C", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Tm", ylab = "A>G")

###### Maturity Lm TL ~ T_C 
class(MATULmmut$Lm_measure)
MATULmmutTL = MATULmmut[!is.na(MATULmmut$Lm_measure),]
MATULmmutTL = MATULmmutTL[MATULmmutTL$Lm_measure == "TL",]
ggscatter(MATULmmutTL, x = "Lm_value", y = "T_C", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Lm TL", ylab = "A>G")

###### Maturity Lm FL ~ T_C 
class(MATULmmut$Lm_measure)
MATULmmutFL = MATULmmut[!is.na(MATULmmut$Lm_measure),]
MATULmmutFL = MATULmmutFL[MATULmmutFL$Lm_measure == "FL",]
ggscatter(MATULmmutFL, x = "Lm_value", y = "T_C", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Lm FL", ylab = "A>G")

###### Maturity Lm ~ T_C 
ggscatter(MATULmmut, x = "Lm_value", y = "T_C", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Lm", ylab = "A>G")
