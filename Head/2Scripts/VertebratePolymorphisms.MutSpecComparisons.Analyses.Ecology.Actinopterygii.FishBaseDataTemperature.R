rm(list=ls(all=TRUE))


### read data
MUT = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt', header = TRUE)
class(MUT$A_T)
class(MUT$Species)

TEMPE = read.table('../../Body/1Raw/FishBaseTemperature.txt', header = TRUE)
class(TEMPE$Temperature)
class(TEMPE$Species)


### calculating average for repeating species in TEMPE
averageTEMPE = aggregate(Temperature ~ ., median, data = TEMPE)

### temperature in fishes 
TemperMut = merge(MUT,TEMPE)   # 1170 rows
averageTemperMut = merge(MUT,averageTEMPE)  # 128 rows

###### Temperature ~ T_C 
cor.test(averageTemperMut$T_C,averageTemperMut$Temperature, method = 'spearman') 

#install.packages("ggpubr")
library("ggpubr")

ggscatter(averageTemperMut, x = "Temperature", y = "T_C", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Temperature", ylab = "A>G")


cor.test(TemperMut$T_A,TemperMut$Temperature, method = 'spearman')
cor.test(TemperMut$T_G,TemperMut$Temperature, method = 'spearman') 

cor.test(TemperMut$A_T,TemperMut$Temperature, method = 'spearman')
cor.test(TemperMut$A_G,TemperMut$Temperature, method = 'spearman') 
cor.test(TemperMut$A_C,TemperMut$Temperature, method = 'spearman')

cor.test(TemperMut$G_A,TemperMut$Temperature, method = 'spearman')
cor.test(TemperMut$G_T,TemperMut$Temperature, method = 'spearman') 
cor.test(TemperMut$G_C,TemperMut$Temperature, method = 'spearman')

cor.test(TemperMut$C_A,TemperMut$Temperature, method = 'spearman')
cor.test(TemperMut$C_T,TemperMut$Temperature, method = 'spearman')
cor.test(TemperMut$C_G,TemperMut$Temperature, method = 'spearman') 

a <- lm(TemperMut$Temperature ~ TemperMut$T_C + TemperMut$T_G + TemperMut$A_G + TemperMut$A_C + TemperMut$G_T + TemperMut$C_G); summary(a)
a <- lm(TemperMut$Temperature ~ TemperMut$T_C + TemperMut$A_G + TemperMut$A_C + TemperMut$G_T + TemperMut$C_G); summary(a)
a <- lm(TemperMut$Temperature ~ TemperMut$T_C + TemperMut$A_G + TemperMut$G_T + TemperMut$C_G); summary(a)
a <- lm(TemperMut$Temperature ~ TemperMut$T_C + TemperMut$G_T + TemperMut$C_G); summary(a)
a <- lm(TemperMut$Temperature ~ scale(TemperMut$T_C) + scale(TemperMut$G_T) + scale(TemperMut$C_G)); summary(a)
# (Intercept)           15.9097     0.2905  54.759  < 2e-16 ***
# scale(TemperMut$T_C)   4.1528     0.3706  11.207  < 2e-16 ***
# scale(TemperMut$G_T)   2.0217     0.4035   5.011 9.04e-07 ***
# scale(TemperMut$C_G)  -1.3195     0.4166  -3.167  0.00169 ** 
# r2 = 0.37

