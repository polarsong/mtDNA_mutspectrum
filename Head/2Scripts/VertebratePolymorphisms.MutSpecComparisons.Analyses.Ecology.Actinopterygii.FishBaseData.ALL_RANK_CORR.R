rm(list=ls(all=TRUE))


##########POLYMORHIC MUTSPEC
MUT = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt', header = TRUE)
class(MUT$A_T)
class(MUT$Species)

##########TEMPERATURE 
TEMPE = read.table('../../Body/1Raw/FishBaseTemperature.txt', header = TRUE)
class(TEMPE$Temperature)
class(TEMPE$Species)
TEMPE = aggregate(Temperature ~ ., median, data = TEMPE)
TemperMut = merge(MUT, TEMPE) 

cor.test(TemperMut$A_T,TemperMut$Temperature, method = 'spearman')   #rho  
cor.test(TemperMut$A_G,TemperMut$Temperature, method = 'spearman')   #rho     -0.3581037 p-value = 3.321e-05
cor.test(TemperMut$A_C,TemperMut$Temperature, method = 'spearman')   #rho   
cor.test(TemperMut$T_A,TemperMut$Temperature, method = 'spearman')   #rho   
cor.test(TemperMut$T_G,TemperMut$Temperature, method = 'spearman')   #rho   
cor.test(TemperMut$T_C,TemperMut$Temperature, method = 'spearman')   #rho     0.2648037 p-value = 0.002522
cor.test(TemperMut$G_A,TemperMut$Temperature, method = 'spearman')   #rho   
cor.test(TemperMut$G_T,TemperMut$Temperature, method = 'spearman')   #rho  
cor.test(TemperMut$G_C,TemperMut$Temperature, method = 'spearman')   #rho  
cor.test(TemperMut$C_A,TemperMut$Temperature, method = 'spearman')   #rho   
cor.test(TemperMut$C_T,TemperMut$Temperature, method = 'spearman')   #rho   
cor.test(TemperMut$C_G,TemperMut$Temperature, method = 'spearman')   #rho   



library("ggpubr")
ggscatter(TemperMut, x = "Temperature", y = "A_G", 
                       add = "reg.line", conf.int = TRUE, 
                       cor.coef = TRUE, cor.method = "spearman",
                       xlab = "Temperature", ylab = "A_G")


##########MATURITY Lm (mean length at first maturity in  )  and Tm (Mean or median age at first maturity)
MATULM = read.table('../../Body/1Raw/FishBaseMaturity_Lm.txt',  header = TRUE, stringsAsFactors=FALSE)
MATUTM = read.table('../../Body/1Raw/FishBaseMaturity_Tm.txt',  header = TRUE)
class(MATULM$Lm)
MATULM$Lm = as.numeric(MATULM$Lm)
MATUTM = aggregate(Tm ~ ., median, data = MATUTM) 
MATULM = aggregate(Lm ~ ., median, data = MATULM)
MATULmmut = merge(MUT,MATULM) 
MATUTmmut = merge(MUT,MATUTM)

cor.test(MATULmmut$A_T,MATULmmut$Lm, method = 'spearman')   #rho  
cor.test(MATULmmut$A_G,MATULmmut$Lm, method = 'spearman')   #rho            
cor.test(MATULmmut$A_C,MATULmmut$Lm, method = 'spearman')   #rho   
cor.test(MATULmmut$T_A,MATULmmut$Lm, method = 'spearman')   #rho   
cor.test(MATULmmut$T_G,MATULmmut$Lm, method = 'spearman')   #rho   
cor.test(MATULmmut$T_C,MATULmmut$Lm, method = 'spearman')   #rho      
cor.test(MATULmmut$G_A,MATULmmut$Lm, method = 'spearman')   #rho   
cor.test(MATULmmut$G_T,MATULmmut$Lm, method = 'spearman')   #rho   
cor.test(MATULmmut$G_C,MATULmmut$Lm, method = 'spearman')   #rho  -0.266859         p-value = 0.008235
cor.test(MATULmmut$C_A,MATULmmut$Lm, method = 'spearman')   #rho     
cor.test(MATULmmut$C_T,MATULmmut$Lm, method = 'spearman')   #rho    
cor.test(MATULmmut$C_G,MATULmmut$Lm, method = 'spearman')   #rho   

cor.test(MATUTmmut$A_T,MATUTmmut$Tm, method = 'spearman')   #rho  
cor.test(MATUTmmut$A_G,MATUTmmut$Tm, method = 'spearman')   #rho            
cor.test(MATUTmmut$A_C,MATUTmmut$Tm, method = 'spearman')   #rho   
cor.test(MATUTmmut$T_A,MATUTmmut$Tm, method = 'spearman')   #rho   
cor.test(MATUTmmut$T_G,MATUTmmut$Tm, method = 'spearman')   #rho   
cor.test(MATUTmmut$T_C,MATUTmmut$Tm, method = 'spearman')   #rho       
cor.test(MATUTmmut$G_A,MATUTmmut$Tm, method = 'spearman')   #rho   
cor.test(MATUTmmut$G_T,MATUTmmut$Tm, method = 'spearman')   #rho   
cor.test(MATUTmmut$G_C,MATUTmmut$Tm, method = 'spearman')   #rho  -0.2219755     p-value = 0.0222
cor.test(MATUTmmut$C_A,MATUTmmut$Tm, method = 'spearman')   #rho     
cor.test(MATUTmmut$C_T,MATUTmmut$Tm, method = 'spearman')   #rho    
cor.test(MATUTmmut$C_G,MATUTmmut$Tm, method = 'spearman')   #rho   




#########ANAGE Maximum.longevity
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



MATULMTM = merge(MATULM, MATUTM)
cor.test(MATULMTM$Lm,MATULMTM$Tm, method = 'spearman')
ggscatter(MATULMTM, x = "Lm", y = "Tm", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Lm", ylab = "Tm")



TEMPMATULM = merge(MATULM, TEMPE)
cor.test(TEMPMATULM$Lm,TEMPMATULM$Temperature, method = 'spearman')
ggscatter(TEMPMATULM, x = "Lm", y = "Temperature", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Lm", ylab = "Temperature")



TEMPMATUTM = merge(MATUTM, TEMPE)
cor.test(TEMPMATUTM$Tm,TEMPMATUTM$Temperature, method = 'spearman')
ggscatter(TEMPMATUTM, x = "Tm", y = "Temperature", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Tm", ylab = "Temperature")

ANAGEMATULM = merge(MATULM, AnAge)
cor.test(ANAGEMATULM$Lm,ANAGEMATULM$Maximum.longevity..yrs., method = 'spearman')
ggscatter(ANAGEMATULM, x = "Lm", y = "Maximum.longevity..yrs.", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Lm", ylab = "Maximum.longevity..yrs.")

ANAGETEMP = merge(TEMPE, AnAge)
cor.test(ANAGETEMP$Temperature,ANAGETEMP$Maximum.longevity..yrs., method = 'spearman')
ggscatter(ANAGETEMP, x = "Temperature", y = "Maximum.longevity..yrs.", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Temperature", ylab = "Maximum.longevity..yrs.")



