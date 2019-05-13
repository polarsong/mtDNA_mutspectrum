rm(list=ls(all=TRUE))

Som = read.table("../../Body/2Derived/HealthySomaticHumanMutations.GTEx.PatientSpecificDerive.txt", header = TRUE, sep = '\t')

################################################################################################
######################## general statistice (to comapre with TCGA)
################################################################################################

par(mfrow = c(2,2))

Subs = data.frame(table(Som$Substitution))
Subs = Subs[order(Subs$Freq),]; Summ = sum(Subs$Freq); Subs$Freq = Subs$Freq/Summ
Subs
#Var1        Freq
#  C_G 0.002339181
#  T_G 0.010136452
#  T_A 0.016374269
#  A_C 0.016764133
#  A_T 0.019103314
#  C_A 0.019103314
#  G_T 0.022222222
#  G_C 0.024561404
#  C_T 0.130994152
#  A_G 0.184405458
#  T_C 0.270175439
#  G_A 0.283820663

Subs = data.frame(table(Som[Som$TurnOverRate <= 64,]$Substitution)); Subs = Subs[order(Subs$Freq),]; Summ = sum(Subs$Freq); Subs$Freq = Subs$Freq/Summ
Subs
#   C_G 0.003353454
#   G_T 0.008048290
#   T_G 0.008718981
#   A_C 0.012072435
#   T_A 0.012072435
#   A_T 0.012743125
#   C_A 0.017437961
#   G_C 0.026827632
#   C_T 0.120724346
#   A_G 0.191817572
#   T_C 0.273641851
#   G_A 0.312541918

Subs = data.frame(table(Som[Som$TurnOverRate > 64,]$Substitution)); Subs = Subs[order(Subs$Freq),]; Summ = sum(Subs$Freq); Subs$Freq = Subs$Freq/Summ
Subs
#  Var1       Freq
#  C_G 0.0009310987
#  T_G 0.0121042831
#  C_A 0.0214152700
#  G_C 0.0214152700
#  T_A 0.0223463687
#  A_C 0.0232774674
#  A_T 0.0279329609
#  G_T 0.0418994413
#  C_T 0.1452513966
#  A_G 0.1741154562
#  G_A 0.2439478585
#  T_C 0.2653631285

# slow cells have: lower T_C, lower A_G (they are together in pca), higher A_T and G_T (they are also together in pca). 

# run prcomp for GTEX (26 tissues)


boxplot(Som[Som$Substitution == 'T_C',]$Position, Som[Som$Substitution == 'G_A',]$Position, Som[Som$Substitution == 'C_T',]$Position, Som[Som$Substitution == 'A_G',]$Position, notch = TRUE, names = c('T>C','G>A','C>T','A>G'))

################################################################################################
######################## count MutSpek and correlate it with Turnover
################################################################################################

Som$G_A = 0; Som[Som$Substitution == 'G_A',]$G_A = 1;
Som$T_C = 0; Som[Som$Substitution == 'T_C',]$T_C = 1;
Som$A_G = 0; Som[Som$Substitution == 'A_G',]$A_G = 1;
Som$C_T = 0; Som[Som$Substitution == 'C_T',]$C_T = 1;
Som$G_T = 0; Som[Som$Substitution == 'G_T',]$G_T = 1;
Som$G_C = 0; Som[Som$Substitution == 'G_C',]$G_C = 1;
Som$T_G = 0; Som[Som$Substitution == 'T_G',]$T_G = 1;
Som$T_A = 0; Som[Som$Substitution == 'T_A',]$T_A = 1;
Som$A_T = 0; Som[Som$Substitution == 'A_T',]$A_T = 1;
Som$A_C = 0; Som[Som$Substitution == 'A_C',]$A_C = 1;
Som$C_A = 0; Som[Som$Substitution == 'C_A',]$C_A = 1;
Som$C_G = 0; Som[Som$Substitution == 'C_G',]$C_G = 1;

Agg = aggregate(list(Som$G_A,Som$G_T,Som$G_C,Som$C_A,Som$C_T,Som$C_G,Som$T_A,Som$T_G,Som$T_C,Som$A_T,Som$A_G,Som$A_C), by = list(Som$TissueShortName,Som$TurnOverRate), FUN = mean)
names(Agg) = c('TissueShortName','TurnOverRate','G_A','G_T','G_C','C_A','C_T','C_G','T_A','T_G','T_C','A_T','A_G','A_C')

cor.test(Agg$TurnOverRate,Agg$G_A, method = 'spearman') # negative:rho -0.4989313; p = 0.008068
ObservedPval = as.numeric(cor.test(Agg$TurnOverRate,Agg$G_A, method = 'spearman')[3]) 
plot(log2(Agg$TurnOverRate),Agg$G_A)

cor.test(Agg$TurnOverRate,Agg$C_T, method = 'spearman') # positive
cor.test(Agg$TurnOverRate,Agg$T_C, method = 'spearman') # zero
cor.test(Agg$TurnOverRate,Agg$A_G, method = 'spearman') # zero

a<-lm(log2(Agg$TurnOverRate) ~ Agg$G_A + Agg$C_T + Agg$T_C + Agg$A_G); summary(a)
a<-lm(log2(Agg$TurnOverRate) ~ Agg$G_A + Agg$C_T + Agg$T_C); summary(a)
a<-lm(log2(Agg$TurnOverRate) ~ Agg$G_A + Agg$C_T); summary(a)
a<-lm(log2(Agg$TurnOverRate) ~ Agg$G_A); summary(a) # R

################################################################################################
######################## make equal (and small) the number of observed mutations in each tissue:
################################################################################################

ShortRes = c()
table(Som$TissueShortName) # the minimum is 8 - uterus
#Adipose Tissue   Adrenal Gland           Blood    Blood Vessel           Brain          Breast           Colon       Esophagus           Heart          Kidney 
#88              83             614             132              54              35              52             212              21              36 
#Liver            Lung          Muscle           Nerve           Ovary        Pancreas       Pituitary        Prostate  Salivary Gland            Skin 
#94              70              81             113              38              36              43              40              13             489 
#Small Intestine          Spleen         Stomach          Testis         Thyroid          Uterus          Vagina 
#15              15              19              40             110               8              14 

### take 10 random mutations (with replacement!) from each tissue and repeat the analysis 1000 times
VecOfTissues = unique(Som$TissueShortName); length(VecOfTissues)

for (permut in 1:1000)
{
  for (tissue in 1:length(VecOfTissues))
    { # tissue = 26
    temp = Som[Som$TissueShortName == VecOfTissues[tissue],]
    temp = temp[sample(nrow(temp),10, replace = TRUE),]
    if (tissue == 1) {Short = temp;}
    if (tissue >  1) {Short = rbind(Short,temp);}
    }
    
  Agg = aggregate(list(Short$G_A,Short$G_T,Short$G_C,Short$C_A,Short$C_T,Short$C_G,Short$T_A,Short$T_G,Short$T_C,Short$A_T,Short$A_G,Short$A_C), by = list(Short$TissueShortName,Short$TurnOverRate), FUN = mean)
  names(Agg) = c('TissueShortName','TurnOverRate','G_A','G_T','G_C','C_A','C_T','C_G','T_A','T_G','T_C','A_T','A_G','A_C')
  
  cor.test(Agg$TurnOverRate,Agg$G_A, method = 'spearman') # negative:rho -0.4989313; p = 0.008068
  ObservedPval = as.numeric(cor.test(Agg$TurnOverRate,Agg$G_A, method = 'spearman')[3]) 
  ObservedRho = as.numeric(cor.test(Agg$TurnOverRate,Agg$G_A, method = 'spearman')[4]) 
  
  ShortRes = rbind(ShortRes,c(ObservedPval,ObservedRho))
}

ShortRes = data.frame(ShortRes); names(ShortRes)=c('Pvalue','Rho')
plot(ShortRes$Rho,-log10(ShortRes$Pvalue))
hist(ShortRes$Rho, breaks = 100)

################################################################################################
######################## what is happening within the same tissue of the same patient (VAF(T_C) ~ VAF(G_A)...)
################################################################################################

#summary(SomSrr$TurnOverRate) # median = 64
# SRR is unique for sample (tissue X patient), so, I need to filter a subset of SRRs with at least one G_A and T_C and compare VAFs
SomSrr = Som[Som$Substitution == 'G_A' | Som$Substitution == 'T_C',] # better correlate the ratio with tissue-specific properties
SomSrr$Number = 1
AggSrr = aggregate(SomSrr$Number, by=list(SomSrr$SRR), FUN = sum); summary(AggSrr$x)
SrrVec = unique(AggSrr[AggSrr$x > 1,]$Group.1); length(SrrVec) # 262 tissues with > 1 substitutions (T_C,G_A)
SomSrr = SomSrr[SomSrr$SRR %in% SrrVec,]
AggVaf = aggregate(SomSrr$AF, by = list(SomSrr$Substitution,SomSrr$SRR,SomSrr$TurnOverRate), FUN = mean)
AggVafGa = AggVaf[AggVaf$Group.1 =='G_A',]; AggVafGa = AggVafGa[,2:4]; names(AggVafGa)=c('SRR','TurnOverRate','GaVaf')
AggVafTc = AggVaf[AggVaf$Group.1 =='T_C',]; AggVafTc = AggVafTc[,2:4]; names(AggVafTc)=c('SRR','TurnOverRate','TcVaf')
Vaf = merge(AggVafGa,AggVafTc, by = c('SRR','TurnOverRate')) # 166 tissues with at least one T_C and G_A
summary(log2(Vaf$TcVaf/Vaf$GaVaf))
hist(log2(Vaf$TcVaf/Vaf$GaVaf), breaks = 20) # take tissues with low cell division - they decrease cell division with time and this effect might be more pronounced
wilcox.test(log2(Vaf$TcVaf/Vaf$GaVaf), mu = 0) # 0.037
Vaf$TcToGaVaf = Vaf$TcVaf/Vaf$GaVaf; summary(Vaf$TcToGaVaf)
cor.test(Vaf$TcToGaVaf,Vaf$TurnOverRate, method = 'spearman')
plot(Vaf$TurnOverRate,Vaf$TcToGaVaf); abline(h = 1, col = 'red', lwd = 2)
cor.test(Vaf$TcVaf,Vaf$TurnOverRate, method = 'spearman')
cor.test(Vaf$GaVaf,Vaf$TurnOverRate, method = 'spearman')

plot(Vaf$TcVaf,Vaf$GaVaf)
cor.test(Vaf$TcVaf,Vaf$GaVaf, method ='spearman') # positive - good
abline(0,1, col = 'red')
plot(Vaf$TcVaf,Vaf$GaVaf, method ='spearman', xlim = c(0,0.2), ylim = c(0,0.2)) # positive - good
abline(0,1, col = 'red')

# the higher the VAF the shorter the turnover rate
cor.test(Vaf$TcVaf,Vaf$TurnOverRate, method = 'spearman')
cor.test(Vaf$GaVaf,Vaf$TurnOverRate, method = 'spearman')
wilcox.test(Vaf$TcVaf,Vaf$GaVaf, paired = TRUE)

cor.test(Vaf$TcToGaVaf,Vaf$TurnOverRate, method = 'spearman')
summary(Vaf$TurnOverRate) # 4.25    30.00    64.00  1958.10   270.00 30000.00 
boxplot(Vaf[Vaf$TurnOverRate <= 30,]$TcToGaVaf,Vaf[Vaf$TurnOverRate > 30 & Vaf$TurnOverRate <= 64,]$TcToGaVaf,Vaf[Vaf$TurnOverRate > 64 & Vaf$TurnOverRate <= 270,]$TcToGaVaf,Vaf[Vaf$TurnOverRate > 270,]$TcToGaVaf, ylim = c(0,5))
abline(h = 1, col = 'red', lwd = 2)
cor.test(Vaf$TcToGaVaf,Vaf$TurnOverRate, method = 'spearman')

################################################################################################
############ permute 1000 times and estimate p value of spearman rank correlation
################################################################################################

PvalueVec = c()
for (i in 1:1000) 
  {
  Som$Substitution = sample(Som$Substitution);
  Som$G_A = 0; Som[Som$Substitution == 'G_A',]$G_A = 1;
  Som$T_C = 0; Som[Som$Substitution == 'T_C',]$T_C = 1;
  Som$A_G = 0; Som[Som$Substitution == 'A_G',]$A_G = 1;
  Som$C_T = 0; Som[Som$Substitution == 'C_T',]$C_T = 1;
  Som$G_T = 0; Som[Som$Substitution == 'G_T',]$G_T = 1;
  Som$G_C = 0; Som[Som$Substitution == 'G_C',]$G_C = 1;
  Som$T_G = 0; Som[Som$Substitution == 'T_G',]$T_G = 1;
  Som$T_A = 0; Som[Som$Substitution == 'T_A',]$T_A = 1;
  Som$A_T = 0; Som[Som$Substitution == 'A_T',]$A_T = 1;
  Som$A_C = 0; Som[Som$Substitution == 'A_C',]$A_C = 1;
  Som$C_A = 0; Som[Som$Substitution == 'C_A',]$C_A = 1;
  Som$C_G = 0; Som[Som$Substitution == 'C_G',]$C_G = 1;
  Agg = aggregate(list(Som$G_A,Som$G_T,Som$G_C,Som$C_A,Som$C_T,Som$C_G,Som$T_A,Som$T_G,Som$T_C,Som$A_T,Som$A_G,Som$A_C), by = list(Som$TissueShortName,Som$TurnOverRate), FUN = mean)
  names(Agg) = c('TissueShortName','TurnOverRate','G_A','G_T','G_C','C_A','C_T','C_G','T_A','T_G','T_C','A_T','A_G','A_C')
  Pvalue = as.numeric(cor.test(Agg$TurnOverRate,Agg$G_A, method = 'spearman')[3])
  PvalueVec = c(PvalueVec,Pvalue)
  }
  
summary(PvalueVec)
hist(PvalueVec, breaks = 50)
length(PvalueVec[PvalueVec<=ObservedPval]) / length(PvalueVec)
