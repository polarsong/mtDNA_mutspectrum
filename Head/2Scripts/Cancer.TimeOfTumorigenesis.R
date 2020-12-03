###################################
###### 26.03.2018: Cancer mut Spectrum - before and after tumorogenesis (05.06.2018)
###################################

#### nucleotide content of the human genome: A T G C: 4993	3871	2159	5357

rm(list=ls(all=TRUE))

ALL = read.table("../../Body/1Raw/mtDNA_snv_Oct2016.txt", head = TRUE, sep = '\t')  # 7611

### DERIVE NECESSARY TRAITS:
ALL$TumorVarFreq = ALL$tumor_reads2/(ALL$tumor_reads1 + ALL$tumor_reads2); summary(ALL$TumorVarFreq)  # 0.01000 0.01738 0.04540 0.20268 0.26278 0.99864
ALL$NormalVarFreq = ALL$normal_reads2/(ALL$normal_reads1 + ALL$normal_reads2); summary(ALL$NormalVarFreq) 

#### pdf out
pdf("../../Body/4Figures/Cancer.TimeOfTumorigenesis01.pdf" , height = 20, width = 15)
par(mfrow=c(2,2), cex = 1.5)

#### DO VAF IN CANCER APPROXIMATE THE TIME OF ORIGIN?
#### if there is a correlation between frequency in normal tissue and in cancer? yes

nrow(ALL)                           # 7611
nrow(ALL[ALL$NormalVarFreq == 0,])  # 2265
nrow(ALL[ALL$NormalVarFreq >  0,])  # 5436
wilcox.test(ALL[ALL$NormalVarFreq >  0,]$TumorVarFreq,ALL[ALL$NormalVarFreq == 0,]$TumorVarFreq) # 0.005325 PAPER
boxplot(ALL[ALL$NormalVarFreq ==  0,]$TumorVarFreq,ALL[ALL$NormalVarFreq > 0,]$TumorVarFreq, outline = FALSE, notch = TRUE, names = c("ZeroInNormal\nN=2265","NonZeroInNormal\nN=5436"), ylab = 'VafInCancer')

cor.test(ALL$NormalVarFreq,ALL$TumorVarFreq,method = 'spearman') # Rho = 0.06958285, p =  1.226e-09 PAPER
cor.test(ALL[ALL$NormalVarFreq >  0,]$NormalVarFreq,ALL[ALL$NormalVarFreq >  0,]$TumorVarFreq, method = 'spearman') # Rho = 0.0883715, p =  9.638e-11 PAPER

ALL$Subs = paste(ALL$ref,ALL$var, sep = '_'); table(ALL$Subs)
VecOfTransitionSubstitutions = c('A_G','G_A','C_T','T_C') # all tr
VecOfTransversionSubstitutions = c('C_A','A_C','C_G','G_C','G_T','T_G','T_A','A_T') # ALL Tv

# CancerType = unique(ALL$Tier2); length(CancerType)
CancerType = unique(ALL$tissue); length(CancerType)

####### TEST 1: Ts/Tv in early and late variants - global analysis
str(ALL$Subs)
str(VecOfTransitionSubstitutions)

summary(ALL$TumorVarFreq) # 0.01000 0.01738 0.04540 0.20268 0.26278 0.99864
A = nrow(ALL[ALL$TumorVarFreq <= 0.04540 & ALL$Subs %in% VecOfTransitionSubstitutions,])
B = nrow(ALL[ALL$TumorVarFreq <= 0.04540 & ALL$Subs %in% VecOfTransversionSubstitutions,])
C = nrow(ALL[ALL$TumorVarFreq > 0.04540 & ALL$Subs %in% VecOfTransitionSubstitutions,])
D = nrow(ALL[ALL$TumorVarFreq > 0.04540 & ALL$Subs %in% VecOfTransversionSubstitutions,])
X = rbind(c(A,B),c(C,D))
fisher.test(X) # 0.689, p = 2.524e-05
mosaicplot(X)
X
# 3469  337
# 3566  239

## many transversions are low quality "noise", so it is important to show, that 
## results are the same if I work with extremely low P-values.
summary(-log10(ALL$pval))  # this is general p-value - that this variant is not mistake?
summary(ALL$somatic_p_value)  # this is the somatic p-value? 
cor.test(ALL$pval,ALL$somatic_p_value, method = 'spearman') # rho = 0.83, p-value < 2.2e-16

# pval
summary(-log10(ALL$pval))
NewAll = ALL[-log10(ALL$pval) > 59,]
nrow(NewAll)                # 5709
median(NewAll$TumorVarFreq) # 0.094

A = nrow(NewAll[NewAll$TumorVarFreq <= 0.094 & NewAll$Subs %in% VecOfTransitionSubstitutions,])
B = nrow(NewAll[NewAll$TumorVarFreq <= 0.094 & NewAll$Subs %in% VecOfTransversionSubstitutions,])
C = nrow(NewAll[NewAll$TumorVarFreq > 0.094 & NewAll$Subs %in% VecOfTransitionSubstitutions,])
D = nrow(NewAll[NewAll$TumorVarFreq > 0.094 & NewAll$Subs %in% VecOfTransversionSubstitutions,])
X = rbind(c(A,B),c(C,D))
fisher.test(X) # 0.667, p = 0.0002002

## instead of Ts/Tv we analyse T>C/G>A - both are high quality transitions - distribution of frequencies is similar, T>C has a bit higher VAF
summary(ALL[ALL$Subs == 'T_C',]$TumorVarFreq)
summary(ALL[ALL$Subs == 'G_A',]$TumorVarFreq)
wilcox.test(ALL[ALL$Subs == 'T_C',]$TumorVarFreq,ALL[ALL$Subs == 'G_A',]$TumorVarFreq) # 0.1577
hist(ALL[ALL$Subs == 'T_C',]$TumorVarFreq, breaks = 50)
hist(ALL[ALL$Subs == 'G_A',]$TumorVarFreq, breaks = 50)

## instead of Ts/Tv we analyse T>C/Tv 
A = nrow(ALL[ALL$TumorVarFreq <= 0.04540 & ALL$Subs == 'T_C',])
B = nrow(ALL[ALL$TumorVarFreq <= 0.04540 & ALL$Subs %in% VecOfTransversionSubstitutions,])
C = nrow(ALL[ALL$TumorVarFreq > 0.04540 & ALL$Subs == 'T_C',])
D = nrow(ALL[ALL$TumorVarFreq > 0.04540 & ALL$Subs %in% VecOfTransversionSubstitutions,])
X = rbind(c(A,B),c(C,D))
fisher.test(X) # 0.681, p = 5.26e-05
X

## instead of Ts/Tv we analyse G>A/Tv 
A = nrow(ALL[ALL$TumorVarFreq <= 0.04540 & ALL$Subs == 'G_A',])
B = nrow(ALL[ALL$TumorVarFreq <= 0.04540 & ALL$Subs %in% VecOfTransversionSubstitutions,])
C = nrow(ALL[ALL$TumorVarFreq > 0.04540 & ALL$Subs == 'G_A',])
D = nrow(ALL[ALL$TumorVarFreq > 0.04540 & ALL$Subs %in% VecOfTransversionSubstitutions,])
X = rbind(c(A,B),c(C,D))
fisher.test(X) # 0.715092, p = 0.00023185
X
# both T>C and G>A similarly contribute to the trend

####### TEST 2: Ts/Tv in early and late - separately for each cancer type

Final = c();
for (Tissue in CancerType)
{  # Tissue = 'Bladder'
    # Temp = ALL[ALL$Tier2 == Tissue,]
    Temp = ALL[ALL$tissue == Tissue,]
    MedVaf = median(Temp$TumorVarFreq)
    TsMedianVaf = median(Temp[Temp$Subs %in% VecOfTransitionSubstitutions,]$TumorVarFreq)
    TvMedianVaf = median(Temp[Temp$Subs %in% VecOfTransversionSubstitutions,]$TumorVarFreq)
    TempLate = Temp[Temp$TumorVarFreq <= median(Temp$TumorVarFreq),]
    TempEarly = Temp[Temp$TumorVarFreq  > median(Temp$TumorVarFreq),]
    LateMed = median(TempLate$TumorVarFreq)
    EarlyMed = median(TempEarly$TumorVarFreq)
    EarlyTs = nrow(TempEarly[TempEarly$Subs %in% VecOfTransitionSubstitutions,])
    # EarlyTs = nrow(TempEarly[TempEarly$Subs == 'G_A',])
    EarlyTv = nrow(TempEarly[TempEarly$Subs %in% VecOfTransversionSubstitutions,])
    LateTs = nrow(TempLate[TempLate$Subs %in% VecOfTransitionSubstitutions,])
    # LateTs = nrow(TempLate[TempLate$Subs == 'G_A',])
    LateTv = nrow(TempLate[TempLate$Subs %in% VecOfTransversionSubstitutions,])
    
    X = rbind(c(LateTs,LateTv),c(EarlyTs,EarlyTv))
    Pvalue = as.numeric(fisher.test(X)[1])
    Odds = as.numeric(fisher.test(X)[3])
    
    Line = c(Tissue,EarlyTs,EarlyTv,LateTs,LateTv,LateMed,EarlyMed,Pvalue,Odds,TsMedianVaf,TvMedianVaf,MedVaf); 
    Final = rbind(Final,Line)
}
Final = data.frame(Final); names(Final)=c('Tissue','EarlyTs','EarlyTv','LateTs','LateTv','LateMed','EarlyMed','Pvalue','Odds','TsMedianVaf','TvMedianVaf','MedVaf')

for (i in 2:12) {Final[,i] = as.numeric(as.character(Final[,i]))}
plot(Final$Odds,Final$Pvalue)
Final = Final[order(Final$Pvalue),]
Final$TsMinusTv = Final$TsMedianVaf - Final$TvMedianVaf
summary(Final$TsMinusTv)
nrow(Final[!is.na(Final$TsMinusTv) & Final$TsMinusTv > 0,]) # 25 
nrow(Final[!is.na(Final$TsMinusTv) & Final$TsMinusTv <= 0,]) # 13  
Final[!is.na(Final$Pvalue) & Final$Pvalue < 0.01,] # only one cancer is nominally significant: PBCA 1.139457e-06 0.06051041; Pediatric Brain Cancer

Final$EarlyTsTv = Final$EarlyTs/Final$EarlyTv
Final$LateTsTv = Final$LateTs/Final$LateTv
Final$AllTsTv = (Final$LateTs+Final$EarlyTs)/(Final$LateTv + Final$EarlyTv)
Final$OddsRatio = Final$LateTsTv/Final$EarlyTsTv
FinalShort = Final[!is.na(Final$OddsRatio) & Final$OddsRatio != Inf,]  # 36 instead of 40
wilcox.test(FinalShort$LateTsTv,FinalShort$EarlyTsTv, paired = TRUE) # 0.0006272,

wilcox.test(FinalShort[FinalShort$Tissue != 'PBCA',]$LateTsTv,FinalShort[FinalShort$Tissue != 'PBCA',]$EarlyTsTv, paired = TRUE) # 0.001051

# plot: one segment - one cancer
plot(NA, xlim=c(0,1), ylim=c(0,72), xlab='VAF', ylab="Ts/Tv")
for (Tissue in CancerType)
{ 
  Temp = FinalShort[FinalShort$Tissue == Tissue,]
  if (nrow(Temp) == 1)
  {
    if (Temp$Pvalue > 0.00001) { segments(Temp$LateMed, Temp$LateTsTv, Temp$EarlyMed, Temp$EarlyTsTv, col = rgb(0.1,0.1,0.1,0.2), lwd = 4)}   # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)}
    if (Temp$Pvalue < 0.00001) { segments(Temp$LateMed, Temp$LateTsTv, Temp$EarlyMed, Temp$EarlyTsTv, col = rgb(1,0.1,0.1,1), lwd = 4)}   # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)}
  }
}
CancerTypeSpecificData = Final

## only T>C

Final = c();
for (Tissue in CancerType)
{  # Tissue = 'Bladder'
  # Temp = ALL[ALL$Tier2 == Tissue,]
  Temp = ALL[ALL$tissue == Tissue,]
  TempLate = Temp[Temp$TumorVarFreq <= median(Temp$TumorVarFreq),]
  TempEarly = Temp[Temp$TumorVarFreq  > median(Temp$TumorVarFreq),]
  LateMed = median(TempLate$TumorVarFreq)
  EarlyMed = median(TempEarly$TumorVarFreq)
  # EarlyTs = nrow(TempEarly[TempEarly$Subs %in% VecOfTransitionSubstitutions,])
  EarlyTs = nrow(TempEarly[TempEarly$Subs == 'T_C',])
  EarlyTv = nrow(TempEarly[TempEarly$Subs %in% VecOfTransversionSubstitutions,])
  # LateTs = nrow(TempLate[TempLate$Subs %in% VecOfTransitionSubstitutions,])
  LateTs = nrow(TempLate[TempLate$Subs == 'T_C',])
  LateTv = nrow(TempLate[TempLate$Subs %in% VecOfTransversionSubstitutions,])
  Line = c(Tissue,EarlyTs,EarlyTv,LateTs,LateTv,LateMed,EarlyMed); 
  Final = rbind(Final,Line)
}
Final = data.frame(Final); names(Final)=c('Tissue','EarlyTs','EarlyTv','LateTs','LateTv','LateMed','EarlyMed')

for (i in 2:7) {Final[,i] = as.numeric(as.character(Final[,i]))}
Final$EarlyTsTv = Final$EarlyTs/Final$EarlyTv
Final$LateTsTv = Final$LateTs/Final$LateTv
Final$AllTsTv = (Final$LateTs+Final$EarlyTs)/(Final$LateTv + Final$EarlyTv)
Final$OddsRatio = Final$LateTsTv/Final$EarlyTsTv
FinalShort = Final[!is.na(Final$OddsRatio) & Final$OddsRatio != Inf,]
nrow(FinalShort)  # 36 instead of 40
wilcox.test(FinalShort$LateTsTv,FinalShort$EarlyTsTv, paired = TRUE) # 0.0004558

## only G>A

Final = c();
for (Tissue in CancerType)
{  # Tissue = 'Bladder'
  # Temp = ALL[ALL$Tier2 == Tissue,]
  Temp = ALL[ALL$tissue == Tissue,]
  TempLate = Temp[Temp$TumorVarFreq <= median(Temp$TumorVarFreq),]
  TempEarly = Temp[Temp$TumorVarFreq  > median(Temp$TumorVarFreq),]
  LateMed = median(TempLate$TumorVarFreq)
  EarlyMed = median(TempEarly$TumorVarFreq)
  # EarlyTs = nrow(TempEarly[TempEarly$Subs %in% VecOfTransitionSubstitutions,])
  EarlyTs = nrow(TempEarly[TempEarly$Subs == 'G_A',])
  EarlyTv = nrow(TempEarly[TempEarly$Subs %in% VecOfTransversionSubstitutions,])
  # LateTs = nrow(TempLate[TempLate$Subs %in% VecOfTransitionSubstitutions,])
  LateTs = nrow(TempLate[TempLate$Subs == 'G_A',])
  LateTv = nrow(TempLate[TempLate$Subs %in% VecOfTransversionSubstitutions,])
  Line = c(Tissue,EarlyTs,EarlyTv,LateTs,LateTv,LateMed,EarlyMed); 
  Final = rbind(Final,Line)
}
Final = data.frame(Final); names(Final)=c('Tissue','EarlyTs','EarlyTv','LateTs','LateTv','LateMed','EarlyMed')

for (i in 2:7) {Final[,i] = as.numeric(as.character(Final[,i]))}
Final$EarlyTsTv = Final$EarlyTs/Final$EarlyTv
Final$LateTsTv = Final$LateTs/Final$LateTv
Final$AllTsTv = (Final$LateTs+Final$EarlyTs)/(Final$LateTv + Final$EarlyTv)
Final$OddsRatio = Final$LateTsTv/Final$EarlyTsTv
FinalShort = Final[!is.na(Final$OddsRatio) & Final$OddsRatio != Inf,]
nrow(FinalShort)  # 36 instead of 40
wilcox.test(FinalShort$LateTsTv,FinalShort$EarlyTsTv, paired = TRUE) # 0.0006133

##### TEST 3: per patient: rank all mutations according to VAF and check probability to have transition as a function of the order

length(unique(ALL$index))
length(unique(ALL$sample))
ALL$Number = 1
AGG = aggregate(ALL$Number, by = list(ALL$sample), FUN = sum); names(AGG)=c('sample','NumberOfMut'); summary(AGG$NumberOfMut) # median = 3
VecOfPatientsWithManyMut = AGG[AGG$NumberOfMut>=2,]$sample; length(VecOfPatientsWithManyMut) # 1253 1715
Final = c()
for (i in 1:length(VecOfPatientsWithManyMut))
{  # i = 3
  Temp = ALL[ALL$sample == VecOfPatientsWithManyMut[i],]
  if (nrow(Temp[Temp$Subs %in% VecOfTransitionSubstitutions,]) > 0 & nrow(Temp[Temp$Subs %in% VecOfTransversionSubstitutions,]) > 0)
  # if (nrow(Temp[Temp$Subs == 'T_C',]) > 0 & nrow(Temp[Temp$Subs == 'G_A',]) > 0)
  {
    VafTs = mean(Temp[Temp$Subs %in% VecOfTransitionSubstitutions,]$TumorVarFreq); 
    # VafTs = mean(Temp[Temp$Subs == 'T_C',]$TumorVarFreq); 
    VafTv = mean(Temp[Temp$Subs %in% VecOfTransversionSubstitutions,]$TumorVarFreq);
    Final=rbind(Final,c(VecOfPatientsWithManyMut[i],VafTs,VafTv))
  }
}

Final = data.frame(Final); names(Final)=c('sample','VafTs','VafTv')
wilcox.test(Final$VafTs,Final$VafTv, paired = TRUE) # p = p-value = 5.109e-15 (419 observations), T>C (p-value = 1.501e-06), G>A (p-value = 2.651e-09), T>C vs G>A (0.009667)

Final$RatioVafTsVafTv = log2(Final$VafTs/Final$VafTv); summary(Final$RatioVafTsVafTv) # higher than zero!!!
wilcox.test(Final$RatioVafTsVafTv, mu = 0)
hist(Final$RatioVafTsVafTv, breaks = 50, main = '', xlab = 'log2(VAF(Ts)/VAF(Tv))') # I can color or plot separately different types of cancers
abline(v = 0, col = 'red', lwd = 3)
Final$DiffVafTsVafTv = Final$VafTs - Final$VafTv; summary(Final$DiffVafTsVafTv) # 0.054


### T>C only:
AGG = aggregate(ALL$Number, by = list(ALL$sample), FUN = sum); names(AGG)=c('sample','NumberOfMut'); summary(AGG$NumberOfMut) # median = 3
VecOfPatientsWithManyMut = AGG[AGG$NumberOfMut>=2,]$sample; length(VecOfPatientsWithManyMut) # 1253 1715
Final = c()
for (i in 1:length(VecOfPatientsWithManyMut))
{  # i = 3
  Temp = ALL[ALL$sample == VecOfPatientsWithManyMut[i],]
  if (nrow(Temp[Temp$Subs == 'T_C',]) > 0 & nrow(Temp[Temp$Subs %in% VecOfTransversionSubstitutions,]) > 0)
  {
    VafTs = mean(Temp[Temp$Subs == 'T_C',]$TumorVarFreq); 
    VafTv = mean(Temp[Temp$Subs %in% VecOfTransversionSubstitutions,]$TumorVarFreq);
    Final=rbind(Final,c(VecOfPatientsWithManyMut[i],VafTs,VafTv))
  }
}
Final = data.frame(Final); names(Final)=c('sample','VafTs','VafTv')
wilcox.test(Final$VafTs,Final$VafTv, paired = TRUE) # p = p-value = 1.501e-06 (286 observations),

### G>A only:
AGG = aggregate(ALL$Number, by = list(ALL$sample), FUN = sum); names(AGG)=c('sample','NumberOfMut'); summary(AGG$NumberOfMut) # median = 3
VecOfPatientsWithManyMut = AGG[AGG$NumberOfMut>=2,]$sample; length(VecOfPatientsWithManyMut) # 1253 1715
Final = c()
for (i in 1:length(VecOfPatientsWithManyMut))
{  # i = 3
  Temp = ALL[ALL$sample == VecOfPatientsWithManyMut[i],]
  if (nrow(Temp[Temp$Subs == 'G_A',]) > 0 & nrow(Temp[Temp$Subs %in% VecOfTransversionSubstitutions,]) > 0)
  {
    VafTs = mean(Temp[Temp$Subs == 'G_A',]$TumorVarFreq); 
    VafTv = mean(Temp[Temp$Subs %in% VecOfTransversionSubstitutions,]$TumorVarFreq);
    Final=rbind(Final,c(VecOfPatientsWithManyMut[i],VafTs,VafTv))
  }
}
Final = data.frame(Final); names(Final)=c('sample','VafTs','VafTv')
wilcox.test(Final$VafTs,Final$VafTv, paired = TRUE) # p = p-value = 2.651e-09 (346 observations),

### T>C/G>A:
AGG = aggregate(ALL$Number, by = list(ALL$sample), FUN = sum); names(AGG)=c('sample','NumberOfMut'); summary(AGG$NumberOfMut) # median = 3
VecOfPatientsWithManyMut = AGG[AGG$NumberOfMut>=2,]$sample; length(VecOfPatientsWithManyMut) # 1253 1715
Final = c()
for (i in 1:length(VecOfPatientsWithManyMut))
{  # i = 3
  Temp = ALL[ALL$sample == VecOfPatientsWithManyMut[i],]
  if (nrow(Temp[Temp$Subs == 'T_C',]) > 0 & nrow(Temp[Temp$Subs == 'G_A',]) > 0)
  {
    VafTs = mean(Temp[Temp$Subs == 'T_C',]$TumorVarFreq); 
    VafTv = mean(Temp[Temp$Subs == 'G_A',]$TumorVarFreq);
    Final=rbind(Final,c(VecOfPatientsWithManyMut[i],VafTs,VafTv))
  }
}
Final = data.frame(Final); names(Final)=c('sample','VafTs','VafTv')
wilcox.test(Final$VafTs,Final$VafTv, paired = TRUE) # p-value = 0.009667 (983 observations),
nrow(Final) # 983

Final$DiffVafTsVafTv = Final$VafTs - Final$VafTv; summary(Final$DiffVafTsVafTv) # 0.045 - higher than zero, 
wilcox.test(Final$DiffVafTsVafTv, mu = 0) # but not significant
Final$RatioVafTsVafTv = log2(Final$VafTs/Final$VafTv); summary(Final$RatioVafTsVafTv) # 0.045 - higher than zero, 
wilcox.test(Final$RatioVafTsVafTv, mu = 0) # not significant

dev.off()
