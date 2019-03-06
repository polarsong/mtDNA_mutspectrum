rm(list=ls(all=TRUE))

MUT = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt', header = TRUE)

GT = read.table('../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt', sep = '\t', header = TRUE)
GT$Species = gsub(' ','_',GT$Scientific_name)

#########
#### MAMMALS: Generation length ~ T>C + A>G
#########

MamGt = merge(GT,MUT, by = 'Species')

cor.test(MamGt$GenerationLength_d,MamGt$A_T, method = 'spearman') # negative
cor.test(MamGt$GenerationLength_d,MamGt$A_G, method = 'spearman') # positive
cor.test(MamGt$GenerationLength_d,MamGt$A_C, method = 'spearman') # negative

cor.test(MamGt$GenerationLength_d,MamGt$T_C, method = 'spearman') # positive
cor.test(MamGt$GenerationLength_d,MamGt$T_A, method = 'spearman')
cor.test(MamGt$GenerationLength_d,MamGt$T_G, method = 'spearman')

cor.test(MamGt$GenerationLength_d,MamGt$G_A, method = 'spearman')
cor.test(MamGt$GenerationLength_d,MamGt$G_T, method = 'spearman') # negative
cor.test(MamGt$GenerationLength_d,MamGt$G_C, method = 'spearman')

cor.test(MamGt$GenerationLength_d,MamGt$C_A, method = 'spearman') # negative
cor.test(MamGt$GenerationLength_d,MamGt$C_T, method = 'spearman')
cor.test(MamGt$GenerationLength_d,MamGt$C_G, method = 'spearman')

a<-lm(MamGt$GenerationLength_d ~ MamGt$A_T +  MamGt$A_G + MamGt$A_C + MamGt$T_C + MamGt$G_T + MamGt$C_A); summary(a)
a<-lm(MamGt$GenerationLength_d ~ MamGt$A_T +  MamGt$A_G + MamGt$A_C + MamGt$T_C + MamGt$G_T); summary(a)
a<-lm(scale(MamGt$GenerationLength_d) ~ scale(MamGt$A_T) +  scale(MamGt$A_G) + scale(MamGt$A_C) + scale(MamGt$T_C) + scale(MamGt$G_T)); summary(a)
a<-lm(MamGt$GenerationLength_d ~ scale(MamGt$A_T) +  scale(MamGt$A_G) + scale(MamGt$A_C) + scale(MamGt$T_C) + scale(MamGt$G_T)); summary(a)
#  (Intercept)       2002.92      75.94  26.376  < 2e-16 ***
#  scale(MamGt$A_T)  -254.80      78.19  -3.259  0.00120 ** 
#  scale(MamGt$A_G)   345.73      78.19   4.422 1.22e-05 *** !!!
#  scale(MamGt$A_C)  -239.33      79.14  -3.024  0.00263 ** 
#  scale(MamGt$T_C)   326.73      80.25   4.071 5.48e-05 *** !!!
#  scale(MamGt$G_T)   256.02      78.44   3.264  0.00118 ** 

### remove effect of ancestral nucleotide frequency (WORKS!!!!!)

MamGt$T_C.NoEffectOfTFreq = MamGt$T_C / (MamGt$T_C + MamGt$T_A + MamGt$T_G); summary(MamGt$T_C.NoEffectOfTFreq)
cor.test(MamGt$GenerationLength_d,MamGt$T_C.NoEffectOfTFreq, method = 'spearman') # positive and significant!!!!

MamGt$A_G.NoEffectOfAFreq = MamGt$A_G / (MamGt$A_G + MamGt$A_T + MamGt$A_C); summary(MamGt$A_G.NoEffectOfAFreq)
cor.test(MamGt$GenerationLength_d,MamGt$A_G.NoEffectOfAFreq, method = 'spearman') # super positive and super significant!!!!

par(mfrow=c(2,2))
plot(log2(MamGt$GenerationLength_d),MamGt$T_C)
plot(log2(MamGt$GenerationLength_d),MamGt$A_G)

#########
#### MAMMALS: Body Temperature  ~ NOTHING (N = 84)
#########

Temper = read.table('../../Body/1Raw/MammalianTemperature.ods.txt', sep = '\t', header = TRUE)
Temper$Species = paste(Temper$FirstName,Temper$SecondName,sep = '_')

MamTemper = merge(Temper,MUT, by = 'Species') # 84 just

cor.test(MamTemper$Temp,MamTemper$A_T, method = 'spearman')
cor.test(MamTemper$Temp,MamTemper$A_G, method = 'spearman')
cor.test(MamTemper$Temp,MamTemper$A_C, method = 'spearman') #  a bit negative
cor.test(MamTemper$Temp,MamTemper$T_A, method = 'spearman')
cor.test(MamTemper$Temp,MamTemper$T_G, method = 'spearman')
cor.test(MamTemper$Temp,MamTemper$T_C, method = 'spearman')
cor.test(MamTemper$Temp,MamTemper$G_A, method = 'spearman'); plot(MamTemper$Temp,MamTemper$G_A)
cor.test(MamTemper$Temp,MamTemper$G_T, method = 'spearman')
cor.test(MamTemper$Temp,MamTemper$G_C, method = 'spearman')
cor.test(MamTemper$Temp,MamTemper$C_A, method = 'spearman')
cor.test(MamTemper$Temp,MamTemper$C_T, method = 'spearman')
cor.test(MamTemper$Temp,MamTemper$C_G, method = 'spearman')

#########
#### MAMMALS: MR  ~ 
#########

Temper = read.table('../../Body/1Raw/MammalianTemperature.ods.txt', sep = '\t', header = TRUE)
Temper$Species = paste(Temper$FirstName,Temper$SecondName,sep = '_')

MamTemper = merge(Temper,MUT, by = 'Species') # 84 just

cor.test(MamTemper$G_A,MamTemper$MR, method = 'spearman'); # a bit negative (non significant)
cor.test(MamTemper$G_T,MamTemper$MR, method = 'spearman');
cor.test(MamTemper$G_C,MamTemper$MR, method = 'spearman');

cor.test(MamTemper$A_T,MamTemper$MR, method = 'spearman');
cor.test(MamTemper$A_G,MamTemper$MR, method = 'spearman');
cor.test(MamTemper$A_C,MamTemper$MR, method = 'spearman');

cor.test(MamTemper$T_A,MamTemper$MR, method = 'spearman');
cor.test(MamTemper$T_G,MamTemper$MR, method = 'spearman');
cor.test(MamTemper$T_C,MamTemper$MR, method = 'spearman');

cor.test(MamTemper$C_A,MamTemper$MR, method = 'spearman'); # a bit negative
cor.test(MamTemper$C_T,MamTemper$MR, method = 'spearman'); # a bit positive
cor.test(MamTemper$C_G,MamTemper$MR, method = 'spearman');

plot(MamTemper$Temp,MamTemper$G_A)



#########
#### MAMMALS: HIBERNATION  ~ NOTHING
#########

Hib = read.table('../../Body/1Raw/HibernatingMammals.txt', sep = '\t')
names(Hib) = c('Species')
Hib$Species = gsub(' ','_',Hib$Species)

MamGt$G_A.NoEffectOfGFreq = MamGt$G_A / (MamGt$G_A + MamGt$G_T + MamGt$G_C); summary(MamGt$G_A.NoEffectOfGFreq)

nrow(MamGt[MamGt$Species %in% Hib$Species,]) # N = 35

boxplot(MamGt[MamGt$Species %in% Hib$Species,]$G_A,MamGt[!MamGt$Species %in% Hib$Species,]$G_A, outline = FALSE, notch = TRUE)
boxplot(MamGt[MamGt$Species %in% Hib$Species,]$G_A.NoEffectOfGFreq,MamGt[!MamGt$Species %in% Hib$Species,]$G_A.NoEffectOfGFreq, outline = FALSE, notch = TRUE)
MamGt = MamGt[MamGt$G_A.NoEffectOfGFreq < 1,]
boxplot(MamGt[MamGt$Species %in% Hib$Species,]$G_A.NoEffectOfGFreq,MamGt[!MamGt$Species %in% Hib$Species,]$G_A.NoEffectOfGFreq, outline = FALSE, notch = TRUE)
wilcox.test(MamGt[MamGt$Species %in% Hib$Species,]$G_A.NoEffectOfGFreq,MamGt[!MamGt$Species %in% Hib$Species,]$G_A.NoEffectOfGFreq, alternative = 'less')

MamGtTemper = merge(MamGt,Temper) # 78
a<-lm(MamGtTemper$G_A ~ MamGtTemper$GenerationLength_d*MamGtTemper$Temp); summary(a)

#########
#### MAMMALS: BMR from AN AGE
#########

AnAge = read.table('../../Body/1Raw/anage_data.txt', sep = '\t', header = TRUE)
AnAge$Species = paste(AnAge$Genus,AnAge$Species,sep = '_')
AnAgeMammals = AnAge[AnAge$Class == 'Mammalia',]
AnAgeMammals = merge(MUT,AnAgeMammals)

cor.test(AnAgeMammals$G_A,AnAgeMammals$Temperature..K., method = 'spearman')
cor.test(AnAgeMammals$G_A,AnAgeMammals$Metabolic.rate..W., method = 'spearman') # a bit negative
cor.test(AnAgeMammals$A_G,AnAgeMammals$Metabolic.rate..W., method = 'spearman') #
cor.test(AnAgeMammals$T_C,AnAgeMammals$Metabolic.rate..W., method = 'spearman') # a bit positive

cor.test(AnAgeMammals$G_A,AnAgeMammals$Body.mass..g., method = 'spearman') # a bit negative!

#########
#### MAMMALS: Body Mass 
#########

GT = read.table('../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt', sep = '\t', header = TRUE)
GT$Species = gsub(' ','_',GT$Scientific_name)

MamGt = merge(GT,MUT, by = 'Species')

cor.test(MamGt$G_A,MamGt$AdultBodyMass_g, method = 'spearman') # nothing (negative sign)
cor.test(MamGt$T_C,MamGt$AdultBodyMass_g, method = 'spearman') # a bit positive
