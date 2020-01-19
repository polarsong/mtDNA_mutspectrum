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

summary(test)

#####Call:
#lm(formula = T_C ~ Temperature + Tm, data = allparameters)

#Residuals:
#  Min       1Q   Median       3Q      Max 
#-0.11116 -0.03176 -0.01421  0.02424  0.16122 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.0150943  0.0145573   1.037   0.3014    
#Temperature  0.0064884  0.0007277   8.916 1.04e-15 ***
#  Tm          -0.0036721  0.0018461  -1.989   0.0484 *  
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.04324 on 160 degrees of freedom
#Multiple R-squared:  0.3982,	Adjusted R-squared:  0.3907 
#F-statistic: 52.94 on 2 and 160 DF,  p-value: < 2.2e-16
