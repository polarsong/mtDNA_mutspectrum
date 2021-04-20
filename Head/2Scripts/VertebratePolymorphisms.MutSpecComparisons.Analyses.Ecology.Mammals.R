rm(list=ls(all=TRUE))

pdf('../../Body/4Figures/VertebratePolymorphisms.MutSpecComparisons.Analyses.Ecology.Mammals.R.01.pdf', width = 70, height = 50)
par(mfrow=c(2,2))

MUT = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegCytb.txt', header = TRUE) # OnlyFourFoldDegCytb
MUT$TsTv = (MUT$T_C + MUT$C_T + MUT$G_A + MUT$A_G) / (MUT$T_A + MUT$A_T + MUT$G_C + MUT$C_G + MUT$G_T + MUT$T_G + MUT$C_A + MUT$A_C)
MUT$TC_TATGTC = (MUT$T_C) / (MUT$T_A + MUT$T_G + MUT$T_C)
summary(MUT$TC_TATGTC)
summary(MUT$TsTv)
MUT = MUT[MUT$TsTv < Inf,]
MUT = MUT[!is.na(MUT$TC_TATGTC),]

GT = read.table('../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt', sep = '\t', header = TRUE)
GT$Species = gsub(' ','_',GT$Scientific_name)

MamGt = merge(GT,MUT, by = 'Species') # 426
Qual = data.frame(table(MamGt$Species))
UniqueSpecies = Qual[Qual$Freq == 1,]$Var1; length(UniqueSpecies);
MamGt = MamGt[MamGt$Species %in% UniqueSpecies,]
MamGt$TC_TCGA = MamGt$T_C / (MamGt$T_C + MamGt$G_A) #  424

######## add the number of mutations per species
Number = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.txt')
Number = data.frame(Number$Species,Number$NumOfFourFoldMutInCytB); names(Number)=c('Species','NumOfFourFoldMutInCytB')
Number = unique(Number)
nrow(MamGt)
MamGt = merge(MamGt,Number, by ='Species')
nrow(MamGt)

summary(MamGt$NumOfFourFoldMutInCytB); # dev.off()
hist(MamGt$NumOfFourFoldMutInCytB,breaks=100)
MamGt = MamGt[MamGt$NumOfFourFoldMutInCytB >= 15,]
hist(MamGt$NumOfFourFoldMutInCytB,breaks=100)

######## TC_TCGA and Generation length
summary(MamGt$TC_TCGA) 
cor.test(MamGt$TC_TCGA,MamGt$GenerationLength_d, method = 'spearman')

######## TsTv and Generation length -
nrow(MamGt) 
cor.test(MamGt$TsTv,MamGt$GenerationLength_d, method = 'spearman')
a <- lm(log2(MamGt$TsTv) ~ log2(MamGt$GenerationLength_d)); 
summary(a)
plot(log2(MamGt$GenerationLength_d),log2(MamGt$TsTv))   ### PAPER
abline(a, col = 'red', lwd = 2)

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

cor.test(MamGt$GenerationLength_d,MamGt$T_C, method = 'spearman')
a<-lm(MamGt$T_C ~ log2(MamGt$GenerationLength_d))
summary(a)
plot(log2(MamGt$GenerationLength_d),MamGt$T_C)
abline(a, col = 'red', lwd = 2)

#### ADD NUMBER OF MUTATIONS USED TO ESTIMATE MUT SPEC
summary(MamGt$NumOfFourFoldMutInCytB)
a<-lm(log2(MamGt$GenerationLength_d) ~ scale(MamGt$A_T) +  scale(MamGt$T_C) + scale(MamGt$G_T) + scale(MamGt$NumOfFourFoldMutInCytB)); summary(a)

### remove effect of ancestral nucleotide frequency (WORKS!!!!!)
MamGt$T_C.NoEffectOfTFreq = MamGt$T_C / (MamGt$T_C + MamGt$T_A + MamGt$T_G); summary(MamGt$T_C.NoEffectOfTFreq)
cor.test(MamGt$GenerationLength_d,MamGt$T_C.NoEffectOfTFreq, method = 'spearman') # positive and significant!!!! rho = 0.164, p = 0.0007114

########### BOXPLOTS BY FAMILIES
par(mar = c(10, 10, 10, 10))
par(cex.lab=10) # is for y-axis
par(cex.axis=10) # is for x-axis
##### boxplots by quartiles
boxplot(MamGt[MamGt$GenerationLength_d<=quantile(MamGt$GenerationLength_d,0.25),]$TsTv,
        MamGt[MamGt$GenerationLength_d>quantile(MamGt$GenerationLength_d,0.25) & MamGt$GenerationLength_d<=quantile(MamGt$GenerationLength_d,0.5),]$TsTv,
        MamGt[MamGt$GenerationLength_d>quantile(MamGt$GenerationLength_d,0.5) & MamGt$GenerationLength_d<=quantile(MamGt$GenerationLength_d,0.75),]$TsTv,
        MamGt[MamGt$GenerationLength_d>quantile(MamGt$GenerationLength_d,0.75),]$TsTv,
        names=c('1','2','3','4'), outline = FALSE, notch = TRUE)

wilcox.test(
        MamGt[MamGt$GenerationLength_d<=quantile(MamGt$GenerationLength_d,0.25),]$TsTv,
        MamGt[MamGt$GenerationLength_d>quantile(MamGt$GenerationLength_d,0.25) & MamGt$GenerationLength_d<=quantile(MamGt$GenerationLength_d,0.5),]$TsTv)
wilcox.test(
        MamGt[MamGt$GenerationLength_d<=quantile(MamGt$GenerationLength_d,0.25),]$TsTv,
        MamGt[MamGt$GenerationLength_d>quantile(MamGt$GenerationLength_d,0.5) & MamGt$GenerationLength_d<=quantile(MamGt$GenerationLength_d,0.75),]$TsTv)
wilcox.test(
        MamGt[MamGt$GenerationLength_d<=quantile(MamGt$GenerationLength_d,0.25),]$TsTv,
        MamGt[MamGt$GenerationLength_d>quantile(MamGt$GenerationLength_d,0.75),]$TsTv)
wilcox.test(
        MamGt[MamGt$GenerationLength_d>quantile(MamGt$GenerationLength_d,0.25) & MamGt$GenerationLength_d<=quantile(MamGt$GenerationLength_d,0.5),]$TsTv,
        MamGt[MamGt$GenerationLength_d>quantile(MamGt$GenerationLength_d,0.5) & MamGt$GenerationLength_d<=quantile(MamGt$GenerationLength_d,0.75),]$TsTv)
wilcox.test(
        MamGt[MamGt$GenerationLength_d>quantile(MamGt$GenerationLength_d,0.25) & MamGt$GenerationLength_d<=quantile(MamGt$GenerationLength_d,0.5),]$TsTv,
        MamGt[MamGt$GenerationLength_d>quantile(MamGt$GenerationLength_d,0.75),]$TsTv)


##### boxplots by median
boxplot(MamGt[MamGt$GenerationLength_d<=quantile(MamGt$GenerationLength_d,0.5),]$TsTv,
        MamGt[MamGt$GenerationLength_d>quantile(MamGt$GenerationLength_d,0.5),]$TsTv,
        names=c('short-lived','long-lived'), outline = FALSE, notch = TRUE, cex = 3)
quantile(MamGt$GenerationLength_d,0.5) # 1497 days
wilcox.test(MamGt[MamGt$GenerationLength_d<=quantile(MamGt$GenerationLength_d,0.5),]$TsTv,MamGt[MamGt$GenerationLength_d>quantile(MamGt$GenerationLength_d,0.5),]$TsTv) # 1.687e-06

##### boxplots by families
par(mfrow=c(1,1))
### derive families 
taxa = read.table("../../Body/1Raw/TaxaFromKostya.Names.stat", sep = "\t") 
taxa = as.data.frame(taxa[grepl('Mammalia',taxa$V1),])
names(taxa) = 'taxa'
taxa$Species = gsub(";(.*)",'',taxa$taxa);
taxa$Species = gsub(" ",'_',taxa$Species);
taxa$Family = gsub(";Mammalia(.*)",'',taxa$taxa)
taxa$Family = gsub("(.*);",'',taxa$Family)
table(taxa$Family)
### merge (perfect merge! no species without taxa)
nrow(MamGt)
MamGt = merge(MamGt,taxa,by='Species') #, all.x=TRUE)
nrow(MamGt)
for (i in 1:nrow(MamGt))   {  MamGt$FamilyShort[i] = paste(unlist(strsplit(MamGt$Family[i],''))[c(1:3)],collapse = '')  }

FamFreq = data.frame(table(MamGt$Family));
FrequentFamilies = FamFreq[FamFreq$Freq >= 3,]$Var1; length(FrequentFamilies) # 3!!!
FamFreq = FamFreq[FamFreq$Var1 %in% FrequentFamilies,]
names(FamFreq)=c('Family','NumberOfSpecies')

Mammalia = MamGt[MamGt$Family %in% FrequentFamilies,]
agg = aggregate(list(Mammalia$TsTv,Mammalia$GenerationLength_d), by = list(Mammalia$Family), FUN = median)
names(agg) = c('Family','TsTv','GenerationLength_d')
cor.test(agg$TsTv,agg$GenerationLength_d,method = 'spearman') ### 0.001649, Rho = 0.8545455 PAPER!!!
plot(agg$GenerationLength_d,agg$TsTv)
agg = merge(agg,FamFreq)
agg = agg[order(agg$GenerationLength_d),]
agg$TsTv = round(agg$TsTv,2)
agg$GenerationLength_d = round(agg$GenerationLength_d,0)


library(boxplotdbl) # install.packages('boxplotdbl')
X = data.frame(as.factor(Mammalia$FamilyShort),Mammalia$TsTv)
Y = data.frame(as.factor(Mammalia$FamilyShort),Mammalia$GenerationLength_d)
par(mar = c(10, 10, 10, 10))
boxplotdou(Y,X,ylim = c(0,50), xlim = c(0,8000), name.on.axis = FALSE, cex = 6, pch = 0, cex.lab = 5, cex.axis = 5, col = rainbow(11)) # name.on.axis = FALSE factor.labels = FALSE,  draw.legend = TRUE

X = data.frame(as.factor(Mammalia$FamilyShort),Mammalia$T_C)
Y = data.frame(as.factor(Mammalia$FamilyShort),Mammalia$GenerationLength_d)
par(mar = c(10, 10, 10, 10))
boxplotdou(Y,X, xlim = c(0,8000), name.on.axis = FALSE, cex = 6, pch = 0, cex.lab = 5, cex.axis = 5, col = rainbow(11)) # name.on.axis = FALSE factor.labels = FALSE,  draw.legend = TRUE

X = data.frame(as.factor(Mammalia$FamilyShort),Mammalia$TC_TCGA)
Y = data.frame(as.factor(Mammalia$FamilyShort),Mammalia$GenerationLength_d)
par(mar = c(10, 10, 10, 10))
boxplotdou(Y,X, xlim = c(0,8000), name.on.axis = FALSE, cex = 6, pch = 0, cex.lab = 5, cex.axis = 5, col = rainbow(11)) # name.on.axis = FALSE factor.labels = FALSE,  draw.legend = TRUE

X = data.frame(as.factor(Mammalia$FamilyShort),Mammalia$TC_TATGTC)
Y = data.frame(as.factor(Mammalia$FamilyShort),Mammalia$GenerationLength_d)
par(mar = c(10, 10, 10, 10))
boxplotdou(Y,X, xlim = c(0,8000), name.on.axis = FALSE, cex = 6, pch = 0, cex.lab = 5, cex.axis = 5, col = rainbow(11)) # name.on.axis = FALSE factor.labels = FALSE,  draw.legend = TRUE

# we can see that the strong increasing effect is starting from the animals with > 1000 days of generation length
dev.off()

##########################################################################################
#################### PICs, Alina+Kostya

library(ape)
library(geiger)
library(caper)

names(MamGt)
summary(MamGt$GenerationLength_d)
tree <- read.tree("../../Body/1Raw/mtalign.aln.treefile.rooted")
row.names(MamGt) = MamGt$Species
tree_w = treedata(tree, MamGt[, c('Species', 'TsTv', 'T_C', 'TC_TCGA', 'G_A','TC_TATGTC','GenerationLength_d','NumOfFourFoldMutInCytB')], 
                sort=T, warnings=T)$phy
data<-as.data.frame(treedata(tree_w, MamGt[, c('Species', 'TsTv', 'T_C', 'TC_TCGA','G_A','TC_TATGTC','GenerationLength_d','NumOfFourFoldMutInCytB')], 
                 sort=T, warnings=T)$data)
nrow(data)
data$Species = as.character(data$Species)
data$TsTv = as.numeric(as.character(data$TsTv))
data$T_C = as.numeric(as.character(data$T_C))
data$G_A = as.numeric(as.character(data$G_A))
data$TC_TCGA = as.numeric(as.character(data$TC_TCGA))
data$TC_TATGTC = as.numeric(as.character(data$TC_TATGTC))
data$GenerationLength_d = as.numeric(as.character(data$GenerationLength_d))
data$NumOfFourFoldMutInCytB = as.numeric(as.character(data$NumOfFourFoldMutInCytB))

MutComp = comparative.data(tree_w, data, Species, vcv=TRUE)
summary(pgls(TsTv ~ log2(GenerationLength_d), MutComp, lambda="ML")) # 1A
summary(pgls(TsTv ~ log2(GenerationLength_d) + log2(NumOfFourFoldMutInCytB), MutComp, lambda="ML")) # 1B
summary(pgls(TsTv ~ 0 + log2(GenerationLength_d) + log2(NumOfFourFoldMutInCytB), MutComp, lambda="ML")) # 1C
summary(pgls(T_C ~ log2(GenerationLength_d), MutComp, lambda="ML")) # 2A
summary(pgls(T_C ~ 0 + log2(GenerationLength_d), MutComp, lambda="ML")) # 2B
summary(pgls(T_C ~ 0 + log2(GenerationLength_d)  + log2(NumOfFourFoldMutInCytB), MutComp, lambda="ML")) # 2C
summary(pgls(TC_TCGA ~ log2(GenerationLength_d), MutComp, lambda="ML")) # 3A
summary(pgls(TC_TCGA ~ 0 + log2(GenerationLength_d), MutComp, lambda="ML")) # 3B
summary(pgls(TC_TCGA ~ 0 + log2(GenerationLength_d) + log2(NumOfFourFoldMutInCytB), MutComp, lambda="ML")) # 3C
summary(pgls(TC_TATGTC ~ log2(GenerationLength_d), MutComp, lambda="ML")) # 4A
summary(pgls(TC_TATGTC ~ log2(GenerationLength_d) + log2(NumOfFourFoldMutInCytB), MutComp, lambda="ML")) # 4B

####################################################################################
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

######## add the number of mutations per species
Number = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.txt')
Number = data.frame(Number$Species,Number$NumOfFourFoldMutInCytB); names(Number)=c('Species','NumOfFourFoldMutInCytB')
Number = unique(Number) # 2118
nrow(MamGt) # 424
MamGt = merge(MamGt,Number, by ='Species')
nrow(MamGt) # 424

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

MATRIX = cbind(MATRIX,MamGt[,c(1:14,26,27)])

###### ANALYSIS OF COMPONENTS
plot(log2(MATRIX$GenerationLength_d),MATRIX$Pca2)
cor.test(MATRIX$Pca1,MATRIX$GenerationLength_d, method = 'spearman') # rho = 0.004288973, p = 0.9298
cor.test(MATRIX$Pca2,MATRIX$GenerationLength_d, method = 'spearman') # rho = -0.3897472, p =  < 2.2e-16
cor.test(MATRIX$Pca3,MATRIX$GenerationLength_d, method = 'spearman') # rho = -0.1037697, p =  0.03266

cor.test(MATRIX$Pca1,MATRIX$NumOfFourFoldMutInCytB, method = 'spearman') # rho = -0.02800981, p = 0.5652
cor.test(MATRIX$Pca2,MATRIX$NumOfFourFoldMutInCytB, method = 'spearman') # rho = 0.2247093, p = 2.964e-06 PAPER

a<-lm(MATRIX$Pca2 ~ scale(MATRIX$GenerationLength_d) + scale(MATRIX$NumOfFourFoldMutInCytB)); summary(a)  # PAPER
a<-lm(MATRIX$Pca2 ~ scale(MATRIX$GenerationLength_d)*scale(MATRIX$NumOfFourFoldMutInCytB)); summary(a)    # PAPER

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
# 1.603726e+00 1.227148e+00 1.204381e+00 1.057175e+00 1.036204e+00 9.831485e-01 9.157673e-01 8.374509e-01 8.318695e-01 7.727840e-01 6.961052e-01 2.487852e-15
PCA$rotation # and the loadings (res$rotation).

dev.off()

#### IF WE RERUN PRCOMP WITH MANY MUTATIONS (N = 211)

summary(MamGt$NumOfFourFoldMutInCytB) # why minumum is 5!!???
MamGt = MamGt[MamGt$NumOfFourFoldMutInCytB > 60,]
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

MATRIX = cbind(MATRIX,MamGt[,c(1:14,26,27)])

PCA$sdev # the eigenvalues (res$sdev) giving information on the magnitude of each PC, 
PCA$rotation # still the first component 
cor.test(MATRIX$Pca2,MATRIX$GenerationLength_d, method = 'spearman') # rho = -0.3897472, p =  < 2.2e-16

#########################################################################################
######################### PIC wit PCA2

library(ape)

tree <- read.tree("../../Body/1Raw/mtalign.aln.treefile.rooted")

data = MATRIX[which(as.character(MATRIX$Species) %in% tree$tip.label),]

df_vec <- as.character(MATRIX$Species)
tree_vec <- tree$tip.label

a <- setdiff(df_vec, tree_vec)
b <- setdiff(tree_vec, df_vec)

tree2 <- drop.tip(tree, b)

TempData = data[, c('Pca2', 'GenerationLength_d')]
contrasts <- as.data.frame(apply(TempData, 2, pic, tree2))
names(contrasts) = names(TempData)

cor.test(contrasts$Pca2, log(contrasts$GenerationLength_d), method = 'spearman')
# rho = -0.1240444, pvalue = 0.1487
