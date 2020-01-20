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
MATU=MATU[!is.na(MATU$Tm),]                           #maturity table
MATU$Lm = NULL
MATU = aggregate(Tm ~ ., mean, data = MATU) 

##multiple reg 

allparameters=merge(averageTemperMut, MATU)


ltest = lm(formula = T_C ~ scale(Temperature)*scale(Tm), data = allparameters)
summary(ltest)
#Call:
#lm(formula = T_C ~ scale(Temperature) * scale(Tm), data = allparameters)

#Residuals:
#  Min        1Q    Median        3Q       Max 
#-0.126180 -0.039390 -0.002721  0.029179  0.162737 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                   0.121096   0.010621  11.401 3.11e-12 ***
#  scale(Temperature)            0.025915   0.010736   2.414   0.0223 *  
#  scale(Tm)                    -0.008237   0.010698  -0.770   0.4476    
#scale(Temperature):scale(Tm)  0.011877   0.014888   0.798   0.4315    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.05871 on 29 degrees of freedom
#Multiple R-squared:  0.2016,	Adjusted R-squared:  0.119 
#F-statistic: 2.441 on 3 and 29 DF,  p-value: 0.08438


