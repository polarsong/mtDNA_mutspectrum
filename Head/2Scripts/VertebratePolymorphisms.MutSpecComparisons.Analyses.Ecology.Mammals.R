rm(list=ls(all=TRUE))

pdf('../../Body/4Figures/VertebratePolymorphisms.MutSpecComparisons.Analyses.Ecology.Mammals.R.01.pdf')

# MUT = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt', header = TRUE)
MUT = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegCytb.txt', header = TRUE) # OnlyFourFoldDegCytb
MUT$TsTv = (MUT$T_C + MUT$C_T + MUT$G_A + MUT$A_G) / (MUT$T_A + MUT$A_T + MUT$G_C + MUT$C_G + MUT$G_T + MUT$T_G + MUT$C_A + MUT$A_C)
summary(MUT$TsTv)
MUT = MUT[MUT$TsTv < Inf,]

GT = read.table('../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt', sep = '\t', header = TRUE)
GT$Species = gsub(' ','_',GT$Scientific_name)

Temper = read.table('../../Body/1Raw/BodyParametersForAllClasses.txt', header = TRUE)

MamGt = merge(GT,MUT, by = 'Species') # 426
Qual = data.frame(table(MamGt$Species))
UniqueSpecies = Qual[Qual$Freq == 1,]$Var1; length(UniqueSpecies);
MamGt = MamGt[MamGt$Species %in% UniqueSpecies,]

######## TsTv and Generation length - significant

nrow(MamGt)  # 424
cor.test(MamGt$TsTv,MamGt$GenerationLength_d, method = 'spearman') # rho = 0.228, p = 2.021e-06
plot(log2(MamGt$GenerationLength_d),log2(MamGt$TsTv))

#### Generation length ~ T>C + A>G

cor.test(MamGt$GenerationLength_d,MamGt$A_T, method = 'spearman') # negative, - 0.22
cor.test(MamGt$GenerationLength_d,MamGt$A_G, method = 'spearman') 
cor.test(MamGt$GenerationLength_d,MamGt$A_C, method = 'spearman') # negative, - 0.17

cor.test(MamGt$GenerationLength_d,MamGt$T_C, method = 'spearman') # positive, rhp = 0.25
cor.test(MamGt$GenerationLength_d,MamGt$T_A, method = 'spearman')
cor.test(MamGt$GenerationLength_d,MamGt$T_G, method = 'spearman') 

cor.test(MamGt$GenerationLength_d,MamGt$G_A, method = 'spearman')
cor.test(MamGt$GenerationLength_d,MamGt$G_T, method = 'spearman') # negative, -0.27
cor.test(MamGt$GenerationLength_d,MamGt$G_C, method = 'spearman') 

cor.test(MamGt$GenerationLength_d,MamGt$C_A, method = 'spearman') # negative, - 0.23
cor.test(MamGt$GenerationLength_d,MamGt$C_T, method = 'spearman')
cor.test(MamGt$GenerationLength_d,MamGt$C_G, method = 'spearman')

a<-lm(MamGt$GenerationLength_d ~ MamGt$A_T +  MamGt$A_G + MamGt$A_C + MamGt$T_C + MamGt$T_G + MamGt$G_T + MamGt$C_A); summary(a) # T>C is the most significant
a<-lm(MamGt$GenerationLength_d ~ MamGt$A_T +  MamGt$A_G + MamGt$A_C + MamGt$T_C + MamGt$G_T + MamGt$C_A); summary(a) # T>C is the most significant
a<-lm(MamGt$GenerationLength_d ~ MamGt$A_T +  MamGt$A_G + MamGt$T_C + MamGt$G_T + MamGt$C_A); summary(a) # T>C is the most significant
a<-lm(MamGt$GenerationLength_d ~ MamGt$A_T +  MamGt$T_C + MamGt$G_T + MamGt$C_A); summary(a) # T>C is the most significant
a<-lm(MamGt$GenerationLength_d ~ MamGt$A_T +  MamGt$T_C + MamGt$G_T); summary(a) # T>C is the most significant
a<-lm(log2(MamGt$GenerationLength_d) ~ scale(MamGt$A_T) +  scale(MamGt$T_C) + scale(MamGt$G_T)); summary(a) # T>C is the most significant#
# (Intercept)      10.39025    0.05584 186.064  < 2e-16 ***
#  scale(MamGt$A_T) -0.18037    0.05603  -3.219  0.00138 ** 
#  scale(MamGt$T_C)  0.26205    0.05734   4.570 6.42e-06 ***
#  scale(MamGt$G_T) -0.18386    0.05723  -3.213  0.00142 ** 
  ---

### remove effect of ancestral nucleotide frequency (WORKS!!!!!)

MamGt$T_C.NoEffectOfTFreq = MamGt$T_C / (MamGt$T_C + MamGt$T_A + MamGt$T_G); summary(MamGt$T_C.NoEffectOfTFreq)
cor.test(MamGt$GenerationLength_d,MamGt$T_C.NoEffectOfTFreq, method = 'spearman') # positive and significant!!!! rho = 0.164, p = 0.0007114

plot(log2(MamGt$GenerationLength_d),MamGt$T_C)

##### Body Mass (N = 426)

cor.test(MamGt$AdultBodyMass_g,MamGt$A_T, method = 'spearman')
cor.test(MamGt$AdultBodyMass_g,MamGt$A_C, method = 'spearman') # negative
cor.test(MamGt$AdultBodyMass_g,MamGt$A_G, method = 'spearman')

cor.test(MamGt$AdultBodyMass_g,MamGt$T_A, method = 'spearman')
cor.test(MamGt$AdultBodyMass_g,MamGt$T_G, method = 'spearman')
cor.test(MamGt$AdultBodyMass_g,MamGt$T_C, method = 'spearman') # positive

cor.test(MamGt$AdultBodyMass_g,MamGt$G_A, method = 'spearman')
cor.test(MamGt$AdultBodyMass_g,MamGt$G_T, method = 'spearman') # negative
cor.test(MamGt$AdultBodyMass_g,MamGt$G_C, method = 'spearman')

cor.test(MamGt$AdultBodyMass_g,MamGt$C_A, method = 'spearman') # negative
cor.test(MamGt$AdultBodyMass_g,MamGt$C_T, method = 'spearman')
cor.test(MamGt$AdultBodyMass_g,MamGt$C_G, method = 'spearman')

a<-lm(log2(MamGt$AdultBodyMass_g) ~ MamGt$A_C + MamGt$T_C + MamGt$G_T + MamGt$C_A); summary(a)   # T>C is the most significant
a<-lm(log2(MamGt$AdultBodyMass_g) ~ MamGt$T_C + MamGt$G_T + MamGt$C_A); summary(a)   # T>C is the most significant
a<-lm(log2(MamGt$AdultBodyMass_g) ~ MamGt$T_C + MamGt$C_A); summary(a)   # T>C is the most significant
a<-lm(log2(MamGt$AdultBodyMass_g) ~ scale(MamGt$T_C) + scale(MamGt$C_A)); summary(a)   # T>C is the most significant


#########
#### MAMMALS: Body Temperature  (N = 84) - nothing
#########

table(Temper$Class)
Temper = Temper[Temper$Class == 'Mammalia',]
MamTemper = merge(Temper,MUT, by = 'Species') # 84 just

cor.test(MamTemper$Temp,MamTemper$A_T, method = 'spearman')
cor.test(MamTemper$Temp,MamTemper$A_G, method = 'spearman')
cor.test(MamTemper$Temp,MamTemper$A_C, method = 'spearman') #  a bit negative
cor.test(MamTemper$Temp,MamTemper$T_A, method = 'spearman')
cor.test(MamTemper$Temp,MamTemper$T_G, method = 'spearman')
cor.test(MamTemper$Temp,MamTemper$T_C, method = 'spearman')
cor.test(MamTemper$Temp,MamTemper$G_A, method = 'spearman')
cor.test(MamTemper$Temp,MamTemper$G_T, method = 'spearman')
cor.test(MamTemper$Temp,MamTemper$G_C, method = 'spearman')
cor.test(MamTemper$Temp,MamTemper$C_A, method = 'spearman')
cor.test(MamTemper$Temp,MamTemper$C_T, method = 'spearman')
cor.test(MamTemper$Temp,MamTemper$C_G, method = 'spearman')


#########
#### MAMMALS: HIBERNATION  ~ NOTHING
#########

Hib = read.table('../../Body/1Raw/HibernatingMammals.txt', sep = '\t')
names(Hib) = c('Species')
Hib$Species = gsub(' ','_',Hib$Species)

MamGt$G_A.NoEffectOfGFreq = MamGt$G_A / (MamGt$G_A + MamGt$G_T + MamGt$G_C); summary(MamGt$G_A.NoEffectOfGFreq)
nrow(MamGt[MamGt$Species %in% Hib$Species,]) # N = 35
MamGt = MamGt[MamGt$G_A.NoEffectOfGFreq < 1,]
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
#### MAMMALS: PRCOMP
#########

MUT = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegCytb.txt', header = TRUE) # OnlyFourFoldDegCytb
MUT$TsTv = (MUT$T_C + MUT$C_T + MUT$G_A + MUT$A_G) / (MUT$T_A + MUT$A_T + MUT$G_C + MUT$C_G + MUT$G_T + MUT$T_G + MUT$C_A + MUT$A_C)
MUT = MUT[MUT$TsTv < Inf,]
GT = read.table('../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt', sep = '\t', header = TRUE)
GT$Species = gsub(' ','_',GT$Scientific_name)
MamGt = merge(GT,MUT, by = 'Species') # 426

Qual = data.frame(table(MamGt$Species))
UniqueSpecies = Qual[Qual$Freq == 1,]$Var1; length(UniqueSpecies);
MamGt = MamGt[MamGt$Species %in% UniqueSpecies,]

###### PCA 

names(MamGt)
MATRIX = MamGt[,c(14:25)]
row.names(MATRIX)=MamGt$Species
matrix = MATRIX

PCA = prcomp(matrix, center = TRUE, scale = TRUE) #FALSE) # I don't scale because we analyze the same units (fraction from MutSpec) 
print(PCA)  
summary(PCA)
MATRIX$Pca1 = PCA$x[,1]
MATRIX$Pca2 = PCA$x[,2]
MATRIX$Pca3 = PCA$x[,3]
MATRIX$Pca4 = PCA$x[,4]

MATRIX = cbind(MATRIX,MamGt[,c(1:14,26)])

###### ANALYSIS OF COMPONENTS
plot(log2(MATRIX$GenerationLength_d),MATRIX$Pca2)
cor.test(MATRIX$Pca1,MATRIX$GenerationLength_d, method = 'spearman') # rho = 0.004288973, p = 0.9298
cor.test(MATRIX$Pca2,MATRIX$GenerationLength_d, method = 'spearman') # rho = -0.3897472, p =  < 2.2e-16
cor.test(MATRIX$Pca3,MATRIX$GenerationLength_d, method = 'spearman') # rho = -0.3897472, p =  < 2.2e-16

biplot(PCA, choices=c(1,2), col = c('white','black'), cex = 0.8) #  biplot(princomp(USArrests),choices=c(1,3))

MATRIX$SIZE = (log10(MATRIX$GenerationLength_d) - min(log10(MATRIX$GenerationLength_d))) / (max(log10(MATRIX$GenerationLength_d))   - min(log10(MATRIX$GenerationLength_d)))  #  normilaze data to 0-1 range
summary(MATRIX$SIZE)

plot(MATRIX$Pca1,MATRIX$Pca2, pch = 16, cex = 1.5, col = rgb(MATRIX$SIZE,0,0,1))
plot(MATRIX$Pca1,MATRIX$Pca2, pch = 16, cex = 1.5, col = rgb(0,MATRIX$SIZE,0,1), xlim = c(-11,3), ylim = c(-3,4.5)); 
biplot(PCA, choices=c(1,2), scale = 0.5, col = c('grey','black'), cex = 0.5) # , xlim = c(-11,3), ylim = c(-3,4.5))

summary(MATRIX$Pca1)
summary(PCA$x[,1])

# arrows(0,0,0.50722224,0.078411901, col = 'red', lwd = 2)

biplot(PCA, choices=c(1,2), scale = 1, col = c('grey','black'), cex = 0.5)
plot(MATRIX$Pca1,MATRIX$Pca2, pch = 16, cex = 1.5, col = rgb(0,0,MATRIX$SIZE,1))

plot(PCA$x[,1],PCA$x[,2])

PCA$x # PC's
PCA$sdev # the eigenvalues (res$sdev) giving information on the magnitude of each PC, 
PCA$rotation # and the loadings (res$rotation).

dev.off()


