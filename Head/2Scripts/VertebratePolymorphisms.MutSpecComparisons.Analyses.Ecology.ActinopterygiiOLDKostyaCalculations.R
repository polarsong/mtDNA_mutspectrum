rm(list=ls(all=TRUE))

MUT = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt', header = TRUE)

### longevity in fishes - nothing:

AnAge = read.table('../../Body/1Raw/anage_data.txt', sep = '\t', header = TRUE)
AnAge$Species = paste(AnAge$Genus,AnAge$Species,sep = '_')
AnAge = AnAge[AnAge$Class == 'Actinopterygii',]
AnAgeMut = merge(MUT,AnAge) # 110

cor.test(AnAgeMut$A_T,AnAgeMut$Maximum.longevity..yrs., method = 'spearman')
cor.test(AnAgeMut$A_G,AnAgeMut$Maximum.longevity..yrs., method = 'spearman')
cor.test(AnAgeMut$A_C,AnAgeMut$Maximum.longevity..yrs., method = 'spearman')
cor.test(AnAgeMut$T_A,AnAgeMut$Maximum.longevity..yrs., method = 'spearman')
cor.test(AnAgeMut$T_G,AnAgeMut$Maximum.longevity..yrs., method = 'spearman')
cor.test(AnAgeMut$T_C,AnAgeMut$Maximum.longevity..yrs., method = 'spearman')
cor.test(AnAgeMut$G_A,AnAgeMut$Maximum.longevity..yrs., method = 'spearman')
cor.test(AnAgeMut$G_T,AnAgeMut$Maximum.longevity..yrs., method = 'spearman')
cor.test(AnAgeMut$G_C,AnAgeMut$Maximum.longevity..yrs., method = 'spearman')
cor.test(AnAgeMut$C_A,AnAgeMut$Maximum.longevity..yrs., method = 'spearman')
cor.test(AnAgeMut$C_T,AnAgeMut$Maximum.longevity..yrs., method = 'spearman')
cor.test(AnAgeMut$C_G,AnAgeMut$Maximum.longevity..yrs., method = 'spearman')

### temperature in fishes: 

Temper = read.table('../../Body/1Raw/BodyParametersForAllClasses.txt',header = TRUE)
table(Temper$Class)
Temper = Temper[Temper$Class == 'Actinopterygii',]
TemperMut = merge(MUT,Temper) # 320 species!!!
TemperMut$T_C.NormalOnlyByT = TemperMut$T_C / (TemperMut$T_C  + TemperMut$T_A + TemperMut$T_G )

###### Temperature, Mass, MetabolRate
###### Temperature ~ T_C 
cor.test(TemperMut$T_C,TemperMut$Temperature, method = 'spearman') # super positive!!!! rho = 0.566

install.packages("ggpubr")
library("ggpubr")
plotTCtemp = ggscatter(TemperMut, x = "Temperature", y = "T_C", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Temperature", ylab = "A>G")

cor.test(TemperMut$T_C.NormalOnlyByT,TemperMut$Temperature, method = 'spearman') # still positive!!!! rho = 0.285
cor.test(TemperMut$T_A,TemperMut$Temperature, method = 'spearman')
cor.test(TemperMut$T_G,TemperMut$Temperature, method = 'spearman') # negative

cor.test(TemperMut$A_T,TemperMut$Temperature, method = 'spearman')
cor.test(TemperMut$A_G,TemperMut$Temperature, method = 'spearman') # super negative
cor.test(TemperMut$A_C,TemperMut$Temperature, method = 'spearman') # negative

cor.test(TemperMut$G_A,TemperMut$Temperature, method = 'spearman')
cor.test(TemperMut$G_T,TemperMut$Temperature, method = 'spearman') # negative
cor.test(TemperMut$G_C,TemperMut$Temperature, method = 'spearman')

cor.test(TemperMut$C_A,TemperMut$Temperature, method = 'spearman')
cor.test(TemperMut$C_T,TemperMut$Temperature, method = 'spearman')
cor.test(TemperMut$C_G,TemperMut$Temperature, method = 'spearman') # negative

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

###### Mass:
cor.test(TemperMut$T_C,TemperMut$Mass, method = 'spearman')        # negative
cor.test(TemperMut$T_A,TemperMut$Mass, method = 'spearman')        # negative
cor.test(TemperMut$T_G,TemperMut$Mass, method = 'spearman')        # positive

cor.test(TemperMut$A_T,TemperMut$Mass, method = 'spearman')
cor.test(TemperMut$A_G,TemperMut$Mass, method = 'spearman')
cor.test(TemperMut$A_C,TemperMut$Mass, method = 'spearman')

cor.test(TemperMut$G_A,TemperMut$Mass, method = 'spearman')
cor.test(TemperMut$G_T,TemperMut$Mass, method = 'spearman')
cor.test(TemperMut$G_C,TemperMut$Mass, method = 'spearman')

cor.test(TemperMut$C_A,TemperMut$Mass, method = 'spearman')
cor.test(TemperMut$C_T,TemperMut$Mass, method = 'spearman')
cor.test(TemperMut$C_G,TemperMut$Mass, method = 'spearman')       # positive

a <- lm(log2(TemperMut$Mass) ~ TemperMut$T_C + TemperMut$T_A + TemperMut$T_G + TemperMut$C_G); summary(a)
a <- lm(log2(TemperMut$Mass) ~ TemperMut$T_C + TemperMut$T_A + TemperMut$T_G); summary(a)
a <- lm(log2(TemperMut$Mass) ~ scale(TemperMut$T_C) + scale(TemperMut$T_A) + scale(TemperMut$T_G)); summary(a)

# (Intercept)           7.80409    0.08928  87.414  < 2e-16 ***
# scale(TemperMut$T_C) -0.35655    0.09264  -3.849 0.000144 ***
# scale(TemperMut$T_A) -0.48299    0.09001  -5.366 1.56e-07 ***
# scale(TemperMut$T_G)  0.39520    0.09274   4.261 2.68e-05 ***
# r2 = 0.17

###### MetabolRate:
cor.test(TemperMut$T_C,TemperMut$MetabolRate, method = 'spearman')  # negative
cor.test(TemperMut$T_A,TemperMut$MetabolRate, method = 'spearman')  # negative
cor.test(TemperMut$T_G,TemperMut$MetabolRate, method = 'spearman')  # positive

cor.test(TemperMut$A_T,TemperMut$MetabolRate, method = 'spearman')
cor.test(TemperMut$A_G,TemperMut$MetabolRate, method = 'spearman')   # positive
cor.test(TemperMut$A_C,TemperMut$MetabolRate, method = 'spearman')

cor.test(TemperMut$G_A,TemperMut$MetabolRate, method = 'spearman')
cor.test(TemperMut$G_T,TemperMut$MetabolRate, method = 'spearman')
cor.test(TemperMut$G_C,TemperMut$MetabolRate, method = 'spearman')

cor.test(TemperMut$C_A,TemperMut$MetabolRate, method = 'spearman')
cor.test(TemperMut$C_T,TemperMut$MetabolRate, method = 'spearman')  # positive
cor.test(TemperMut$C_G,TemperMut$MetabolRate, method = 'spearman')  # positive

a <- lm(log2(TemperMut$MetabolRate) ~ TemperMut$T_C + TemperMut$T_A + TemperMut$T_G + TemperMut$A_G + TemperMut$C_T  + TemperMut$C_G); summary(a)
a <- lm(log2(TemperMut$MetabolRate) ~ TemperMut$T_C + TemperMut$T_A + TemperMut$A_G + TemperMut$C_T  + TemperMut$C_G); summary(a)
a <- lm(log2(TemperMut$MetabolRate) ~ TemperMut$T_A + TemperMut$A_G + TemperMut$C_T  + TemperMut$C_G); summary(a)
a <- lm(log2(TemperMut$MetabolRate) ~ scale(TemperMut$T_A) + scale(TemperMut$A_G) + scale(TemperMut$C_T)  + scale(TemperMut$C_G)); summary(a)
# (Intercept)           3.57036    0.07808  45.726  < 2e-16 ***
# scale(TemperMut$T_A) -0.65710    0.11770  -5.583 5.11e-08 ***
# scale(TemperMut$A_G) -0.33857    0.10539  -3.212  0.00145 ** 
# scale(TemperMut$C_T)  0.27527    0.09305   2.958  0.00333 ** 
# scale(TemperMut$C_G)  0.70015    0.11048   6.338 8.07e-10 ***
# r2 = 0.14
