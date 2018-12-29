###### 3.04.2018: Cancer cases with extremely high number of mutations (breast cancer - cataesis)
###### 27.03.2018: Cancer mut Spectrum - asymmetry (strand bias) PAPER
###### 26.03.2018: Cancer mut Spectrum - before and after tumorogenesis 
###### 18.03.2018: Cancer number of divisions (from Table S1 one before last column Tomasetti and Vogelstein 2015 science) and T>C fraction (from bioarchive Campbell 2018): PAPER

###################################
###### 3.04.2018: Cancer cases with extremely high number of mutations (breast cancer - cataesis)
###################################

rm(list=ls(all=TRUE))
C = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/CancerMtDnaSomaticMutations/mtDNA_snv_Oct2016.txt', head = TRUE, sep = '\t')
head(C)
C$Subs = paste(C$ref,C$var, sep = '_')
C$TotalNumberOfMut = 1
AGG = aggregate(C$TotalNumberOfMut, by = list(C$sample), FUN = sum)
AGG = AGG[order(AGG$x, decreasing = TRUE),]
VecOfHighlyMutatedSamples = c(as.character(AGG[AGG$x > 17,]$Group.1))
for (i in 1:length(VecOfHighlyMutatedSamples))
{ # i = 1
  temp = C[C$sample == VecOfHighlyMutatedSamples[i],]
  FractionT_C = nrow(temp[temp$Subs == 'T_C' | temp$Subs == 'A_G',])/nrow(temp)
  OneLine = data.frame(VecOfHighlyMutatedSamples[i], nrow(temp),FractionT_C)
  if (i == 1) {FINAL = OneLine}
  if (i >  1) {FINAL = rbind(FINAL,OneLine)}
}
FINAL
# VecOfHighlyMutatedSamples.i. nrow.temp. FractionT_C
# 1 2779fa01-ac93-4e80-a997-3385f72172c3         33   0.8181818
# 2                            ICGC_MB34         24   0.3750000
# 3                            ICGC_0441         18   0.7222222

C = C[!C$sample %in% FINAL$VecOfHighlyMutatedSamples.i.,]
FractionT_C = nrow(C[C$Subs == 'T_C' | C$Subs == 'A_G',])/nrow(C)
## fisher test?

###################################
###### 27.03.2018: Cancer mut Spectrum - asymmetry (strand bias) PAPER
###################################

rm(list=ls(all=TRUE))
C = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/CancerMtDnaSomaticMutations/mtDNA_snv_Oct2016.txt', head = TRUE, sep = '\t')
head(C)
table(C$Levin2012) 

Ref = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/HumanMtDna.txt', head= TRUE)
HumanMtDna = unlist(strsplit(as.character(Ref$HumanGenome[1]),''))
length(HumanMtDna)
table(HumanMtDna)

par(mfrow=c(1,3))
for (i in 1:length(HumanMtDna))
{ # i = 1
  temp = HumanMtDna[1:i]
  AtoT = length(temp[temp == 'T'])/(length(temp[temp == 'A']) + length(temp[temp == 'T']))
  OneLine = data.frame(i,AtoT)
  if (i == 1) {Final = OneLine}
  if (i >  1) {Final = rbind(Final,OneLine)}
}
plot(Final$i,Final$AtoT, main = 'T/A', ylim = c(0,1))
abline(h=0.5, col = 'red')
TtoAPosition = Final$i
TtoARatio = Final$AtoT

for (i in 1:length(HumanMtDna))
{ # i = 1
  temp = HumanMtDna[1:i]
  AtoT = length(temp[temp == 'A'])/(length(temp[temp == 'A']) + length(temp[temp == 'T']))
  OneLine = data.frame(i,AtoT)
  if (i == 1) {Final = OneLine}
  if (i >  1) {Final = rbind(Final,OneLine)}
}
plot(Final$i,Final$AtoT, main = 'A/T', ylim = c(0,1))
abline(h=0.5, col = 'red')
AtoTPosition = Final$i
AtoTRatio = Final$AtoT


for (i in 1:length(HumanMtDna))
{ # i = 1
  temp = HumanMtDna[1:i]
  AtoT = length(temp[temp == 'G'])/(length(temp[temp == 'G']) + length(temp[temp == 'C']))
  OneLine = data.frame(i,AtoT)
  if (i == 1) {Final = OneLine}
  if (i >  1) {Final = rbind(Final,OneLine)}
}
plot(Final$i,Final$AtoT, main = 'G/C', ylim = c(0,1))
abline(h=0.5, col = 'red')
GtoCPosition = Final$i
GtoCRatio = Final$AtoT

FrA = 5124;
FrT = 4094;
FrG = 2169;
FrC = 5181;

C$Subs = paste(C$ref,C$var, sep = '_')
table(C$Subs)
C$NumberOfMut = 1

AncA = C[C$ref == 'A',]; AncA$NumberOfMut = 1/FrA;
AncT = C[C$ref == 'T',]; AncT$NumberOfMut = 1/FrT;
AncG = C[C$ref == 'G',]; AncG$NumberOfMut = 1/FrG;
AncC = C[C$ref == 'C',]; AncC$NumberOfMut = 1/FrC;

FIN = rbind(AncA,AncT,AncG,AncC)
# FIN = FIN[FIN$position > 6000,]

pdf('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/4_FIGURES/CancerAsymmetryAlongGenome.pdf', width = 14, height = 14)
par(mfrow=c(3,2))
##### asymmetry along genome - TcToAg: 
TEMP = FIN[FIN$Subs == 'T_C' | FIN$Subs == 'A_G',]
TEMP = TEMP[order(TEMP$position),]
TcToAg = 0.5 # 
for (region in 1:nrow(TEMP))
{ # region = 1
  ShortTemp = TEMP[c(1:region),]
  TcToAg = nrow(ShortTemp[ShortTemp$Subs == 'T_C',])/nrow(ShortTemp[ShortTemp$Subs == 'T_C' | ShortTemp$Subs == 'A_G',])
  OneLine = data.frame(ShortTemp$position[region],TcToAg)
  if (region == 1) {Final = OneLine}
  if (region >  1) {Final = rbind(Final,OneLine)}
}
str(Final)
plot(Final$ShortTemp.position.region.,Final$TcToAg, main = 'level of asymmetry: TcToAg', ylim = c(0,1), xlim = c(0,16000))
par(new=TRUE)
plot(TtoAPosition,TtoARatio, ylab='', xlab='', ylim = c(0,1), xlim = c(0,16000), col = 'red', pch='.')
abline(h=0.5, col = 'red')

##### asymmetry along genome - GaToCt: 
TEMP = FIN[FIN$Subs == 'G_A' | FIN$Subs == 'C_T',]
TEMP = TEMP[order(TEMP$position),]
TcToAg = 0.5 # 
for (region in 1:nrow(TEMP))
{ # region = 1
  ShortTemp = TEMP[c(1:region),]
  TcToAg = nrow(ShortTemp[ShortTemp$Subs == 'G_A',])/nrow(ShortTemp[ShortTemp$Subs == 'C_T' | ShortTemp$Subs == 'G_A',])
  OneLine = data.frame(ShortTemp$position[region],TcToAg)
  if (region == 1) {Final = OneLine}
  if (region >  1) {Final = rbind(Final,OneLine)}
}
str(Final)
plot(Final$ShortTemp.position.region.,Final$TcToAg, main = 'level of asymmetry: GaToCt', ylim = c(0,1), xlim = c(0,16000))
par(new = TRUE)
plot(GtoCPosition,GtoCRatio, main = 'level of asymmetry: GaToCt',ylab='', xlab='', ylim = c(0,1), xlim = c(0,16000), col = 'red', pch = '.')
abline(h=0.5, col = 'red')

##### asymmetry along genome - GcToCg: 
TEMP = FIN[FIN$Subs == 'G_C' | FIN$Subs == 'C_G',]
TEMP = TEMP[order(TEMP$position),]
TcToAg = 0.5 # 
for (region in 1:nrow(TEMP))
{ # region = 1
  ShortTemp = TEMP[c(1:region),]
  TcToAg = nrow(ShortTemp[ShortTemp$Subs == 'G_C',])/nrow(ShortTemp[ShortTemp$Subs == 'C_G' | ShortTemp$Subs == 'G_C',])
  OneLine = data.frame(ShortTemp$position[region],TcToAg)
  if (region == 1) {Final = OneLine}
  if (region >  1) {Final = rbind(Final,OneLine)}
}
str(Final)
plot(Final$ShortTemp.position.region.,Final$TcToAg, main = 'level of asymmetry: GcToCg', ylim = c(0,1), xlim = c(0,16000))
par(new = TRUE)
plot(GtoCPosition,GtoCRatio, ylab='', xlab='', ylim = c(0,1), xlim = c(0,16000), col = 'red', pch = '.')
abline(h=0.5, col = 'red')

##### asymmetry along genome - GtToCa: 
TEMP = FIN[FIN$Subs == 'G_T' | FIN$Subs == 'C_A',]
TEMP = TEMP[order(TEMP$position),]
TcToAg = 0.5 # 
for (region in 1:nrow(TEMP))
{ # region = 1
  ShortTemp = TEMP[c(1:region),]
  TcToAg = nrow(ShortTemp[ShortTemp$Subs == 'G_T',])/nrow(ShortTemp[ShortTemp$Subs == 'C_A' | ShortTemp$Subs == 'G_T',])
  OneLine = data.frame(ShortTemp$position[region],TcToAg)
  if (region == 1) {Final = OneLine}
  if (region >  1) {Final = rbind(Final,OneLine)}
}
str(Final)
plot(Final$ShortTemp.position.region.,Final$TcToAg, main = 'level of asymmetry: GtToCa', ylim = c(0,1), xlim = c(0,16000))
par(new = TRUE)
plot(GtoCPosition,GtoCRatio, ylab='', xlab='', ylim = c(0,1), xlim = c(0,16000), col = 'red', pch = '.')
abline(h=0.5, col = 'red')


##### asymmetry along genome - TaToAt: 
TEMP = FIN[FIN$Subs == 'T_A' | FIN$Subs == 'A_T',]
TEMP = TEMP[order(TEMP$position),]
TcToAg = 0.5 # 
for (region in 1:nrow(TEMP))
{ # region = 1
  ShortTemp = TEMP[c(1:region),]
  TcToAg = nrow(ShortTemp[ShortTemp$Subs == 'T_A',])/nrow(ShortTemp[ShortTemp$Subs == 'T_A' | ShortTemp$Subs == 'A_T',])
  OneLine = data.frame(ShortTemp$position[region],TcToAg)
  if (region == 1) {Final = OneLine}
  if (region >  1) {Final = rbind(Final,OneLine)}
}
str(Final)
plot(Final$ShortTemp.position.region.,Final$TcToAg, main = 'level of asymmetry: TaToAt', ylim = c(0,1), xlim = c(0,16000))
par(new=TRUE)
plot(TtoAPosition,TtoARatio, ylab='', xlab='', ylim = c(0,1), xlim = c(0,16000), col = 'red', pch='.')
abline(h=0.5, col = 'red')


##### asymmetry along genome - AcToTg: 
TEMP = FIN[FIN$Subs == 'A_C' | FIN$Subs == 'T_G',]
TEMP = TEMP[order(TEMP$position),]
TcToAg = 0.5 # 
for (region in 1:nrow(TEMP))
{ # region = 1
  ShortTemp = TEMP[c(1:region),]
  TcToAg = nrow(ShortTemp[ShortTemp$Subs == 'A_C',])/nrow(ShortTemp[ShortTemp$Subs == 'A_C' | ShortTemp$Subs == 'T_G',])
  OneLine = data.frame(ShortTemp$position[region],TcToAg)
  if (region == 1) {Final = OneLine}
  if (region >  1) {Final = rbind(Final,OneLine)}
}
str(Final)
plot(Final$ShortTemp.position.region.,Final$TcToAg, main = 'level of asymmetry: AcToTg', ylim = c(0,1), xlim = c(0,16000))
par(new=TRUE)
plot(AtoTPosition,AtoTRatio, ylab='', xlab='', ylim = c(0,1), xlim = c(0,16000), col = 'red', pch='.')
abline(h=0.5, col = 'red')
dev.off()

# log2(10000) = 13.29
# log2(100000) = 16.6
# log2(1000000) = 19.3
# exp(16)

##### aggregate without taking into account different types of cancer:
AGG = aggregate(FIN$NumberOfMut, by = list(FIN$Subs), FUN = sum)

TcToAg = AGG[AGG$Group.1 == 'T_C',]$x/AGG[AGG$Group.1 == 'A_G',]$x; TcToAg # 5
GaToCt = AGG[AGG$Group.1 == 'G_A',]$x/AGG[AGG$Group.1 == 'C_T',]$x; GaToCt # 14
GcToCg = AGG[AGG$Group.1 == 'G_C',]$x/AGG[AGG$Group.1 == 'C_G',]$x; GcToCg # 12.92047
TaToAt = AGG[AGG$Group.1 == 'T_A',]$x/AGG[AGG$Group.1 == 'A_T',]$x; TaToAt # 0.5041117
GtToCa = AGG[AGG$Group.1 == 'G_T',]$x/AGG[AGG$Group.1 == 'C_A',]$x; GtToCa # 0.7008638
AcToTg = AGG[AGG$Group.1 == 'A_C',]$x/AGG[AGG$Group.1 == 'T_G',]$x; AcToTg # 2.692876

##### aggregate with taking into account different types of cancer:
AGG = aggregate(FIN$NumberOfMut, by = list(FIN$Subs,FIN$Tier2), FUN = sum)
names(AGG) = c('Subs','Tissue','Freq')

TC = AGG[AGG$Subs == 'T_C',]; TC=TC[,c(2,3)]; names(TC)=c('Tissue','FreqTC')
AG = AGG[AGG$Subs == 'A_G',]; AG=AG[,c(2,3)]; names(AG)=c('Tissue','FreqAG')
GA = AGG[AGG$Subs == 'G_A',]; GA=GA[,c(2,3)]; names(GA)=c('Tissue','FreqGA')
CT = AGG[AGG$Subs == 'C_T',]; CT=CT[,c(2,3)]; names(CT)=c('Tissue','FreqCT')
GC = AGG[AGG$Subs == 'G_C',]; GC=GC[,c(2,3)]; names(GC)=c('Tissue','FreqGC')
CG = AGG[AGG$Subs == 'C_G',]; CG=CG[,c(2,3)]; names(CG)=c('Tissue','FreqCG')
TA = AGG[AGG$Subs == 'T_A',]; TA=TA[,c(2,3)]; names(TA)=c('Tissue','FreqTA')
AT = AGG[AGG$Subs == 'A_T',]; AT=AT[,c(2,3)]; names(AT)=c('Tissue','FreqAT')
GT = AGG[AGG$Subs == 'G_T',]; GT=GT[,c(2,3)]; names(GT)=c('Tissue','FreqGT')
CA = AGG[AGG$Subs == 'C_A',]; CA=CA[,c(2,3)]; names(CA)=c('Tissue','FreqCA')
AC = AGG[AGG$Subs == 'A_C',]; AC=AC[,c(2,3)]; names(AC)=c('Tissue','FreqAC')
TG = AGG[AGG$Subs == 'T_G',]; TG=TG[,c(2,3)]; names(TG)=c('Tissue','FreqTG')

ALL = merge(TC,AG, all = TRUE); ALL = merge(ALL,GA,all = TRUE); ALL = merge(ALL,CT,all = TRUE); ALL = merge(ALL,GC,all = TRUE); ALL = merge(ALL,CG,all = TRUE);
ALL = merge(ALL,TA,all = TRUE); ALL = merge(ALL,AT,all = TRUE); ALL = merge(ALL,GT,all = TRUE); ALL = merge(ALL,CA,all = TRUE); ALL = merge(ALL,AC,all = TRUE);  
ALL = merge(ALL,TG,all = TRUE); 

ALL$TcToAg = ALL$FreqTC / ALL$FreqAG; summary(ALL$TcToAg); wilcox.test(ALL$TcToAg, mu = 1) # 6.403e-05
ALL$GaToCt = ALL$FreqGA / ALL$FreqCT; summary(ALL$GaToCt); wilcox.test(ALL$GaToCt, mu = 1) # 6.403e-05
ALL$GcToCg = ALL$FreqGC / ALL$FreqCG; summary(ALL$GcToCg); wilcox.test(ALL$GcToCg, mu = 1) # 0.001592 * 6 = 0.0095
ALL$TaToAt = ALL$FreqTA / ALL$FreqAT; summary(ALL$TaToAt); wilcox.test(ALL$TaToAt, mu = 1) # 0.1219
ALL$GtToCa = ALL$FreqGT / ALL$FreqCA; summary(ALL$GtToCa); wilcox.test(ALL$GtToCa, mu = 1) # 0.3931
ALL$AcToTg = ALL$FreqAC / ALL$FreqTG; summary(ALL$AcToTg); wilcox.test(ALL$AcToTg, mu = 1) # 0.008653 * 6 = 0.0519

pdf('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/4_FIGURES/StrandBias6RatiosCancers.pdf')
boxplot(ALL$TcToAg,ALL$GaToCt,ALL$GcToCg,ALL$TaToAt,ALL$GtToCa,ALL$AcToTg, names = c('TC/AG','GA/CT','GC/CG','TA/AT','GT/CA','AC/TG'))
abline(h = 1, col = 'red')
dev.off()

###################################
###### 26.03.2018: Cancer mut Spectrum - before and after tumorogenesis (05.06.2018)
###################################

#### nucleotide content of the human genome: A T G C: 4993	3871	2159	5357

rm(list=ls(all=TRUE))
ALL = read.table('/media/konstantinpopadin/ac45df81-e084-4d30-9653-5c57cc9b58fd/konstantinpopadin/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/CancerMtDnaSomaticMutations/mtDNA_snv_Oct2016.txt', head = TRUE, sep = '\t')

### DERIVE NECESSARY TRAITS:
ALL$TumorVarFreq = ALL$tumor_reads2/(ALL$tumor_reads1 + ALL$tumor_reads2); summary(ALL$TumorVarFreq)  # 0.01000 0.01738 0.04540 0.20268 0.26278 0.99864
ALL$NormalVarFreq = ALL$normal_reads2/(ALL$normal_reads1 + ALL$normal_reads2); summary(ALL$NormalVarFreq) 
table(ALL$Levin2012) 
#exist_inLevin      notinLevin present_inLevin 
#2542            5012              57 
ALL$Subs = paste(ALL$ref,ALL$var, sep = '_'); table(ALL$Subs)
VecOfTransitionSubstitutions = c('A_G','G_A','C_T','T_C') # all tr
VecOfTransversionSubstitutions = c('C_A','A_C','C_G','G_C','G_T','T_G','T_A','A_T') # ALL Tv
Ts = ALL[ALL$Subs %in% VecOfTransitionSubstitutions,]; nrow(Ts) 
Tv = ALL[ALL$Subs %in% VecOfTransversionSubstitutions,]; nrow(Tv)
Ts$Ts = 1; Tv$Ts = 0;  
ALL = rbind(Ts,Tv)

CancerType = unique(ALL$Tier2); length(CancerType)

####### printing pdf
pdf('/media/konstantinpopadin/ac45df81-e084-4d30-9653-5c57cc9b58fd/konstantinpopadin/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/4_FIGURES/CancerTsTvNewOld01.pdf', height = 15, width = 15)
par(mfrow=c(3,2))

####### APPROACH 1: TsTv in TWO GROUPS of mutations: ZERO IN NORMAL; NON ZERO IN NORMAL
ALL$GradationFromNewToOld = 0
for (i in 1:nrow(ALL))
{
  if (ALL$NormalVarFreq[i] == 0)  {ALL$GradationFromNewToOld[i] = 0}
  if (ALL$NormalVarFreq[i] > 0)   {ALL$GradationFromNewToOld[i] = 1}
}
table(ALL$GradationFromNewToOld)
# 0    1 
# 2265 5346

Final = c();
for (i in 0:max(ALL$GradationFromNewToOld))
{ # i = 0
  for (JackNife in 0:100)
  {
    Temp = ALL[ALL$GradationFromNewToOld == i,]
    if (JackNife > 0) 
    {
      Temp = ALL[ALL$GradationFromNewToOld == i,]
      Temp = Temp[sample(nrow(Temp),500),]
    }
    TempTs = Temp[Temp$Subs %in% VecOfTransitionSubstitutions,]
    TempTv = Temp[Temp$Subs %in% VecOfTransversionSubstitutions,]
    TsTv = nrow(TempTs)/nrow(TempTv)
    Line = c(i,TsTv,JackNife); names(Line)=c('Group','TsTv','JackNife')
    Final = rbind(Final,Line)
  }
}
Final = data.frame(Final)  #
Final[Final$JackNife == 0,]
# 1       0 11.37705        0
# 102     1 12.60305        0

FinalJackNife = Final[Final$JackNife > 0,]
boxplot(FinalJackNife[FinalJackNife$Group == 0,]$TsTv,FinalJackNife[FinalJackNife$Group == 1,]$TsTv,notch = TRUE, outline = FALSE, names = c('zero in normal (N = 2265)','more-than-zero in normal (N=5346)'), ylab = 'TsTv', main = 'TsTv in new versus old classes of mutations')

####### APPROACH 2: TsTv in four GROUPS of mutations: from rare in cancer (young) to common in cancer (old)
summary(ALL$TumorVarFreq) # 0.01000 0.01738 0.04540 0.20268 0.26278 0.99864
ALL$GradationFromNewToOld = 0
for (i in 1:nrow(ALL))
{
  if (ALL$TumorVarFreq[i] <= 0.01738)  {ALL$GradationFromNewToOld[i] = 0}
  if (ALL$TumorVarFreq[i] > 0.01738 & ALL$TumorVarFreq[i] <= 0.04540)  {ALL$GradationFromNewToOld[i] = 1}
  if (ALL$TumorVarFreq[i] > 0.04540 & ALL$TumorVarFreq[i] <= 0.26278)  {ALL$GradationFromNewToOld[i] = 2}
  if (ALL$TumorVarFreq[i] > 0.26278)  {ALL$GradationFromNewToOld[i] = 3}
}
table(ALL$GradationFromNewToOld)
# 0    1    2    3 
# 1903 1903 1902 1903 

Final = c();
for (i in 0:max(ALL$GradationFromNewToOld))
{ # i = 0
  for (JackNife in 0:100)
  {
    Temp = ALL[ALL$GradationFromNewToOld == i,]
    if (JackNife > 0) 
    {
      Temp = ALL[ALL$GradationFromNewToOld == i,]
      Temp = Temp[sample(nrow(Temp),500),]
    }
    TempTs = Temp[Temp$Subs %in% VecOfTransitionSubstitutions,]
    TempTv = Temp[Temp$Subs %in% VecOfTransversionSubstitutions,]
    TsTv = nrow(TempTs)/nrow(TempTv)
    Line = c(i,TsTv,JackNife); names(Line)=c('Group','TsTv','JackNife')
    Final = rbind(Final,Line)
  }
}
Final = data.frame(Final)  #
Final[Final$JackNife == 0,]
# Group      TsTv JackNife
# 1       0  9.015789        0
# 102     1 11.945578        0
# 203     2 12.208333        0
# 304     3 19.031579        0
FinalJackNife = Final[Final$JackNife > 0,]
boxplot(FinalJackNife[FinalJackNife$Group == 0,]$TsTv,FinalJackNife[FinalJackNife$Group == 1,]$TsTv,FinalJackNife[FinalJackNife$Group == 2,]$TsTv,FinalJackNife[FinalJackNife$Group == 3,]$TsTv,notch = TRUE, outline = FALSE, names = c('first 25% (N=1903)','second 25% (N=1903)','third 25% (N=1902)','fourth 25% (N=1903)'), ylab = 'TsTv', main = 'TsTv in rare (young) versus common (old) cancer mutations')

####### APPROACH 3: reduce 4 groups to 2 (in relation to median) and run it cancer-specifically
summary(ALL$TumorVarFreq) # 0.01000 0.01738 0.04540 0.20268 0.26278 0.99864
ALL$GradationFromNewToOld = 0
for (i in 1:nrow(ALL))
{
  if (ALL$TumorVarFreq[i] <= 0.04540)  {ALL$GradationFromNewToOld[i] = 0}
  if (ALL$TumorVarFreq[i] > 0.04540)  {ALL$GradationFromNewToOld[i] = 1}
}
table(ALL$GradationFromNewToOld)
# 0    1 
# 3806 3805 

Final = c();
for (i in 0:max(ALL$GradationFromNewToOld))
{ # i = 0
  for (JackNife in 0:100)
  {
    Temp = ALL[ALL$GradationFromNewToOld == i,]
    if (JackNife > 0) 
    {
      Temp = ALL[ALL$GradationFromNewToOld == i,]
      Temp = Temp[sample(nrow(Temp),500),]
    }
    TempTs = Temp[Temp$Subs %in% VecOfTransitionSubstitutions,]
    TempTv = Temp[Temp$Subs %in% VecOfTransversionSubstitutions,]
    TsTv = nrow(TempTs)/nrow(TempTv)
    Line = c(i,TsTv,JackNife); names(Line)=c('Group','TsTv','JackNife')
    Final = rbind(Final,Line)
  }
}
Final = data.frame(Final)  #
Final[Final$JackNife == 0,]
#Group     TsTv JackNife
#1       0 10.29377        0
#102     1 14.92050        0
#FinalJackNife = Final[Final$JackNife > 0,]
boxplot(FinalJackNife[FinalJackNife$Group == 1,]$TsTv,FinalJackNife[FinalJackNife$Group == 0,]$TsTv,notch = TRUE, outline = FALSE, names = c('high VAF (N=3806)','low VAF (N=3805)'), ylab = 'TsTv', main = 'TsTv in mutations with high and low VAFs')
wilcox.test(FinalJackNife[FinalJackNife$Group == 1,]$TsTv,FinalJackNife[FinalJackNife$Group == 0,]$TsTv); # 4.636e-07

####### APPROACH 4: do it separately for each cancer type: median in each case should be cancer- specific!!!!

Final = c();
for (Tissue in CancerType)
{  # Tissue = 'Bladder'
    Temp = ALL[ALL$Tier2 == Tissue,]
    TempLate = Temp[Temp$TumorVarFreq <= median(Temp$TumorVarFreq),]
    TempEarly = Temp[Temp$TumorVarFreq  > median(Temp$TumorVarFreq),]
    EarlyTs = nrow(TempEarly[TempEarly$Subs %in% VecOfTransitionSubstitutions,])
    EarlyTv = nrow(TempEarly[TempEarly$Subs %in% VecOfTransversionSubstitutions,])
    LateTs = nrow(TempLate[TempLate$Subs %in% VecOfTransitionSubstitutions,])
    LateTv = nrow(TempLate[TempLate$Subs %in% VecOfTransversionSubstitutions,])
    Line = c(Tissue,EarlyTs,EarlyTv,LateTs,LateTv); 
    Final = rbind(Final,Line)
}
Final = data.frame(Final); names(Final)=c('Tissue','EarlyTs','EarlyTv','LateTs','LateTv')
for (i in 2:5) {Final[,i] = as.numeric(as.character(Final[,i]))}
Final$EarlyTsTv = Final$EarlyTs/Final$EarlyTv
Final$LateTsTv = Final$LateTs/Final$LateTv
### there is one Infinity (19/0) - paired test will be able to deal with it? Inf is very big number.. should be ok.
wilcox.test(Final$LateTsTv,Final$EarlyTsTv, paired = TRUE) # 0.0106

## plot segments: dev.off() one cancer is missing here
plot(NA, xlim=c(-0.2,1.2), ylim=c(0,70), xlab='', ylab="TsTv")
for (Tissue in CancerType)
{ # Tissue = 'Bladder'
  Temp = Final[Final$Tissue == Tissue,]
  segments(0, Temp$LateTsTv, 1, Temp$EarlyTsTv, col = rgb(0.1,0.1,0.1,0.1), lwd = 3) # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
}
segments(0, mean(Final$LateTsTv), 1, mean(Final[!is.infinite(Final$EarlyTsTv),]$EarlyTsTv), col = rgb(1,0.1,0.1,1), lwd = 3) # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
Final
#Tissue EarlyTs EarlyTv LateTs LateTv  LateTsTv EarlyTsTv
#1          Bladder      47       2     41      9  4.555556 23.500000
#2  Bone/SoftTissue      72       5     74      4 18.500000 14.400000
#3           Breast     321      23    316     29 10.896552 13.956522
#4          Biliary      43       2     41      4 10.250000 21.500000
#5           Cervix      19       0     16      3  5.333333       Inf
#6         Lymphoid     145      13    145     13 11.153846 11.153846
#7          Myeloid      32       7     37      2 18.500000  4.571429
#8     Colon/Rectum     105       7     98     14  7.000000 15.000000
#9         Prostate     355      14    359     11 32.636364 25.357143
#10       Esophagus     187      17    180     25  7.200000 11.000000
#11         Stomach      93      12     95     10  9.500000  7.750000
#12             CNS     128       3    105     26  4.038462 42.666667
#13       Head/Neck      67       6     67      6 11.166667 11.166667
#14          Kidney     372      39    349     62  5.629032  9.538462
#15           Liver     596      40    589     48 12.270833 14.900000
#16            Lung     130       4    126      8 15.750000 32.500000
#17           Ovary     194       7    188     13 14.461538 27.714286
#18        Pancreas     404      26    397     34 11.676471 15.538462
#19            Skin     106       2    101      8 12.625000 53.000000
#20         Thyroid      79       5     78      6 13.000000 15.800000
#21          Uterus      71       1     67      6 11.166667 71.000000

###### APPROACH 5: RUN IT PATIENT-WISE: rank all mutations according to VAF and check probability to have transition as a function of the order.

length(unique(ALL$index))
length(unique(ALL$sample))
ALL$Number = 1
AGG = aggregate(ALL$Number, by = list(ALL$sample), FUN = sum); names(AGG)=c('sample','NumberOfMut'); summary(AGG$NumberOfMut) # median = 3
VecOfPatientsWithManyMut = AGG[AGG$NumberOfMut>=2,]$sample; length(VecOfPatientsWithManyMut) # 1253 1715
Final = c()
for (i in 1:length(VecOfPatientsWithManyMut))
{ 
  Temp = ALL[ALL$sample == VecOfPatientsWithManyMut[i],]
  if (nrow(Temp[Temp$Ts == 1,]) > 0 & nrow(Temp[Temp$Ts == 0,]) > 0)
  {
    VafTs = mean(Temp[Temp$Ts == '1',]$TumorVarFreq); 
    VafTv = mean(Temp[Temp$Ts == '0',]$TumorVarFreq);
    Final=rbind(Final,c(VecOfPatientsWithManyMut[i],VafTs,VafTv))
  }
}

Final = data.frame(Final); names(Final)=c('sample','VafTs','VafTv')
Final$RatioVafTsVafTv = log2(Final$VafTs/Final$VafTv); summary(Final$RatioVafTsVafTv) # higher than zero!!!
hist(Final$RatioVafTsVafTv, breaks = 50, main = '', xlab = 'log2(VAF(Ts)/VAF(Tv))') # I can color or plot separately different types of cancers
abline(v = 0, col = 'red', lwd = 3)
wilcox.test(Final$VafTs,Final$VafTv, paired = TRUE) # p = p-value = 5.109e-15 (419 observations)
dev.off()  


FreqA = 4993;
FreqT = 3871;
FreqG = 2159;
FreqC = 5357;
total = FreqA + FreqT + FreqG + FreqC
FreqA = FreqA/total; FreqT = FreqT/total; FreqG = FreqG/total; FreqC = FreqC/total;
FreqA  # 0.304823 
FreqT  # 0.2363248
FreqG  # 0.1318071
FreqC  # 0.3270452
FreqA+FreqT+FreqG+FreqC # 1

plot(ALL$NormalVarFreq, ALL$TumorVarFreq)

summary(ALL$NormalVarFreq)                         #  0.0000000 0.0000000 0.0003605 0.0006286 0.0007522 0.0285714
summary(ALL[ALL$NormalVarFreq > 0,]$NormalVarFreq) #  3.116e-05 3.315e-04 5.513e-04 8.950e-04 1.011e-03 2.857e-02
nrow(ALL[ALL$NormalVarFreq == 0,]) # 2265
nrow(ALL[ALL$NormalVarFreq > 0 & ALL$NormalVarFreq < 5.513e-04,]) # 2673
nrow(ALL[ALL$NormalVarFreq > 5.513e-04,]) # 2673


ALL$NumberOfMut = 1
ALL$NumberOfMutNormalized = 1
#### normilaze by the number of the first nucleotide (the higher the frequency of the first nucleotide, the less should be the probability of a mutation)
#### so, here I just divide by the number of ancestral nucleotide
for (i in 1:nrow(ALL))
{
  if (ALL$ref[i] == 'A') {ALL$NumberOfMutNormalized[i] = 1/4993;}
  if (ALL$ref[i] == 'T') {ALL$NumberOfMutNormalized[i] = 1/3871;}
  if (ALL$ref[i] == 'G') {ALL$NumberOfMutNormalized[i] = 1/2159;}
  if (ALL$ref[i] == 'C') {ALL$NumberOfMutNormalized[i] = 1/5357;}
}

#### CAN WE PREDICT TIME OF ORIGIN OF THE mtDNA SOMATIC MUTATIONS USING 3 TYPES OF INFORMATION: mtDB, Levin and the frequency in normal
#### If yes, overlap should be significant:
ALL$MtInDbNumeric = 0
ALL$MtInLivenNumeric = 0
ALL$MtInBothDbNumeric = 0
ALL$FreqMoreThanZeroInNormal = 0
for (i in 1:nrow(ALL))
{ # i = 1
  if (ALL$mtDB[i] != 'notinDB')    {ALL$MtInDbNumeric[i] = 1; ALL$MtInBothDbNumeric[i] = 1;}
  if (ALL$Levin2012[i] == 'exist_inLevin')    {ALL$MtInLivenNumeric[i] = 1;  ALL$MtInBothDbNumeric[i] = 1;}
  if (ALL$NormalVarFreq[i] > 0)    {ALL$FreqMoreThanZeroInNormal[i] = 1;}
}
table(ALL$MtInLivenNumeric) # 0 - 5069; 1 - 2542
table(ALL$MtInDbNumeric) # 0 - 5609; 1 - 2002
table(ALL$MtInBothDbNumeric) # 0 - 4957; 1 - 2654
summary(ALL$NormalVarFreq)
summary(ALL[ALL$MtInDbNumeric == 0,]$NormalVarFreq) # 0.0003490
summary(ALL[ALL$MtInDbNumeric == 1,]$NormalVarFreq) # 0.0004102


boxplot(ALL[ALL$MtInDbNumeric == 0,]$NormalVarFreq,ALL[ALL$MtInDbNumeric == 1,]$NormalVarFreq, outline = FALSE, notch = TRUE, names = c('NotInDatabases (N=5609)','InDatabase (N=2002)'), ylab = 'Variant Frequency in Normal Tissue', ylim = c(0,0.003), main = 'All Variants')
wilcox.test(ALL[ALL$MtInDbNumeric == 0,]$NormalVarFreq,ALL[ALL$MtInDbNumeric == 1,]$NormalVarFreq) # p = 4.044e-05
nrow(ALL[ALL$MtInDbNumeric == 0 & ALL$NormalVarFreq > 0,]) # 3968
nrow(ALL[ALL$MtInDbNumeric == 1 & ALL$NormalVarFreq > 0,]) # 1378
boxplot(ALL[ALL$MtInDbNumeric == 0 & ALL$NormalVarFreq > 0,]$NormalVarFreq,ALL[ALL$MtInDbNumeric == 1 & ALL$NormalVarFreq > 0,]$NormalVarFreq, outline = FALSE, notch = TRUE, names = c('NotInDatabases (N=3968)','InDatabase (N=1378)'), ylab = 'Variant Frequency in Normal Tissue', ylim = c(0,0.003), main = 'Variants With Freq In Normal Tissue > 0')
wilcox.test(ALL[ALL$MtInDbNumeric == 0 & ALL$NormalVarFreq > 0,]$NormalVarFreq,ALL[ALL$MtInDbNumeric == 1 & ALL$NormalVarFreq > 0,]$NormalVarFreq) # p-value < 2.2e-16

boxplot(ALL[ALL$MtInLivenNumeric == 0,]$NormalVarFreq,ALL[ALL$MtInLivenNumeric == 1,]$NormalVarFreq, outline = FALSE, notch = TRUE, names = c('NotInLevin (N=5069)','InLevin (N=2542)'), ylab = 'Variant Frequency in Normal Tissue', ylim = c(0,0.003), main = 'All Variants')
nrow(ALL[ALL$MtInLivenNumeric == 0 & ALL$NormalVarFreq > 0,]) # 3600
nrow(ALL[ALL$MtInLivenNumeric == 1 & ALL$NormalVarFreq > 0,]) # 1746
wilcox.test(ALL[ALL$MtInLivenNumeric == 0,]$NormalVarFreq,ALL[ALL$MtInLivenNumeric == 1,]$NormalVarFreq) # 0.002115
boxplot(ALL[ALL$MtInLivenNumeric == 0  & ALL$NormalVarFreq > 0,]$NormalVarFreq,ALL[ALL$MtInLivenNumeric == 1  & ALL$NormalVarFreq > 0,]$NormalVarFreq, outline = FALSE, notch = TRUE, names = c('NotInLevin (N=3600)','InLevin (N=1746)'), ylab = 'Variant Frequency in Normal Tissue', ylim = c(0,0.003), main = 'Variants With Freq In Normal Tissue > 0')
wilcox.test(ALL[ALL$MtInLivenNumeric == 0  & ALL$NormalVarFreq > 0,]$NormalVarFreq,ALL[ALL$MtInLivenNumeric == 1  & ALL$NormalVarFreq > 0,]$NormalVarFreq) # 1.038e-15

boxplot(ALL[ALL$MtInBothDbNumeric == 0,]$NormalVarFreq,ALL[ALL$MtInBothDbNumeric == 1,]$NormalVarFreq, outline = FALSE, notch = TRUE, names = c('NotInBothDb (N=4957)','InAtLeastOneDB (N=2654)'), ylab = 'Variant Frequency in Normal Tissue', ylim = c(0,0.003), main = 'All Variants')
boxplot(ALL[ALL$MtInBothDbNumeric == 0  & ALL$NormalVarFreq > 0,]$NormalVarFreq,ALL[ALL$MtInBothDbNumeric == 1  & ALL$NormalVarFreq > 0,]$NormalVarFreq, outline = FALSE, notch = TRUE, names = c('NotInBothDb (N=3529)','InAtLeastOneDB (N=1817)'), ylab = 'Variant Frequency in Normal Tissue', ylim = c(0,0.003), main = 'Variants With Freq In Normal Tissue > 0')

# ALL$GradationFromNewToOld = ALL$MtInDbNumeric  + ALL$MtInLivenNumeric + ALL$FreqMoreThanZeroInNormal
table(ALL$GradationFromNewToOld) # 1428 4366 1817
#  0    1    2    3 
# 1428 3783 1093 1307 


###### by cancer type
Final = c();
CancerType = unique(ALL$Tier2); length(CancerType)

for (i in 0:1) # 2, 3
{ # i = 0
  for (Tissue in CancerType)
    {
    # Temp = ALL[ALL$GradationFromNewToOld == i,] # 
    Temp = ALL[ALL$NormalFreqGroups == i,]
    Temp = Temp[Temp$Tier2 == Tissue,]
    
    TempTs = Temp[Temp$Subs %in% VecOfTransitionSubstitutions,]
    TempTv = Temp[Temp$Subs %in% VecOfTransversionSubstitutions,]
    TsTv = nrow(TempTs)/nrow(TempTv)
    Line = c(i,nrow(TempTs),nrow(TempTv),Tissue); names(Line)=c('Group','Ts','Tv','Tissue')
    Final = rbind(Final,Line)
    }
}
Final = data.frame(Final)  # can we jacknife it to derive variations? Can we get p-values of the deviation from the neutral (randomly shuffle 0,1,2,3 and for each rank get boxplot of expectations?)!!!!
Final$Ts = as.numeric(as.character(Final$Ts))
Final$Tv = as.numeric(as.character(Final$Tv))
Final$Group = as.numeric(as.character(Final$Group))
Final$TsTv = Final$Ts/Final$Tv
Final = Final[order(Final$Tissue,Final$Group),]
FinalNoInf = Final[Final$TsTv != Inf,]

summary(FinalNoInf$TsTv)
FinalNoInf$TotalMuts = FinalNoInf$Ts + FinalNoInf$Tv; 
summary(FinalNoInf$TotalMuts)  # 
quantile(FinalNoInf$TotalMuts,0.1)  # 28
FinalNoInf = FinalNoInf[FinalNoInf$TotalMuts >= 30,]
nrow(FinalNoInf) # 36

## plot segments:

AGG = aggregate(FinalNoInf$TsTv, by = list(FinalNoInf$Group), FUN = mean); names(AGG)=c('Group','TsTv')
plot(NA, xlim=c(0,1), ylim=c(0,100), xlab='', ylab="TsTv")
for (Tissue in CancerType)
{ 
  Temp = FinalNoInf[FinalNoInf$Tissue == Tissue,]
  for (count in 1:(nrow(Temp)-1))
   {#  count = 1 
   segments(Temp$Group[count], Temp$TsTv[count], Temp$Group[count+1], Temp$TsTv[count+1], col = rgb(0.1,0.1,0.1,0.1), lwd = 3) # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
   }
}
# mean
for (count in 1:(nrow(AGG)-1))
{#  count = 1 
  segments(AGG$Group[count], AGG$TsTv[count], AGG$Group[count+1], AGG$TsTv[count+1], col = rgb(1,0.1,0.1,1), lwd = 3) # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
}
dev.off()

### paired test (TsTv in the left part compare with TsTv on the right part of each segment)
Final = c()
for (Tissue in CancerType)
{
  Temp = FinalNoInf[FinalNoInf$Tissue == Tissue,]
  Temp = Temp[order(Temp$Group),]
  Vec = Temp$TsTv
  if (length(Vec) == 2) {Line = c(Vec[1],Vec[2]); Final = rbind(Final,Line)}
  if (length(Vec) == 3) {Line = c(Vec[1],Vec[2]); Final = rbind(Final,Line); Line = c(Vec[2],Vec[3]); Final = rbind(Final,Line);}
  if (length(Vec) == 4) {Line = c(Vec[1],Vec[2]); Final = rbind(Final,Line); Line = c(Vec[2],Vec[3]); Final = rbind(Final,Line); Line = c(Vec[3],Vec[4]); Final = rbind(Final,Line);}
}
Final = data.frame(Final); names(Final)=c('LeftPartOfSegment','RightPartOfSegment'); nrow(Final) # 40
wilcox.test(Final$LeftPartOfSegment,Final$RightPartOfSegment, paired = TRUE) # P = 0.0007414 # if TC only

###################################### barplots (GC & TA are zero in 2-nd and 3-rd classes!!!! it is not good. Why? these classes are expected to have spectrum similar to normal cells... some bias...)
###################################### do I see changes in spectrum within each patient?
table(ALL$GradationFromNewToOld)
summary(ALL$NumberOfMutNormalized)
table(ALL$Subs)

GoodOrder = c('T_C','G_A','A_G','C_T','A_C','C_G','G_C','T_A','T_G','A_T','C_A','G_T');  # as in mammals
VecOfColors = c('green4','green','green3','greenyellow','steelblue1','steelblue3','royalblue1','royalblue4','orchid1','orchid4','purple','purple4') # 'mediumblue','navy') # as in mammals

Group0 = ALL[ALL$GradationFromNewToOld == 0,]; Group0Agg = aggregate(Group0$NumberOfMutNormalized, by = list(Group0$Subs), FUN = sum); names(Group0Agg)=c('Subs','Class0'); Total = sum(Group0Agg$Class0); Group0Agg$Class0 = Group0Agg$Class0/Total; sum(Group0Agg$Class0); 
Group1 = ALL[ALL$GradationFromNewToOld == 1,]; Group1Agg = aggregate(Group1$NumberOfMutNormalized, by = list(Group1$Subs), FUN = sum); names(Group1Agg)=c('Subs','Class1'); Total = sum(Group1Agg$Class1); Group1Agg$Class1 = Group1Agg$Class1/Total; sum(Group1Agg$Class1); 
Group2 = ALL[ALL$GradationFromNewToOld == 2,]; Group2Agg = aggregate(Group2$NumberOfMutNormalized, by = list(Group2$Subs), FUN = sum); names(Group2Agg)=c('Subs','Class2'); Total = sum(Group2Agg$Class2); Group2Agg$Class2 = Group2Agg$Class2/Total; sum(Group2Agg$Class2); 
# Group3 = ALL[ALL$GradationFromNewToOld == 3,]; Group3Agg = aggregate(Group3$NumberOfMutNormalized, by = list(Group3$Subs), FUN = sum); names(Group3Agg)=c('Subs','Class3'); Total = sum(Group3Agg$Class3); Group3Agg$Class3 = Group3Agg$Class3/Total; sum(Group3Agg$Class3); 
Groups = merge(Group0Agg,Group1Agg);Groups = merge(Groups,Group2Agg, all.x = TRUE);  
row.names(Groups) = Groups$Subs;
Groups[is.na(Groups)]<-0
Groups = Groups[match(GoodOrder, Groups$Subs),] # order Substitutions as I din for mammals
Groups = Groups[,-c(1)]

par(las=2)
par(mai=c(2.4,0.82,0.82,0.42)) # default par(mai=c(1.02,0.82,0.82,0.42))
FinalForBarPLot = as.matrix(Groups)
barplot(FinalForBarPLot,col=VecOfColors, legend.text=TRUE)
dev.off()


######################################
#### PLOTS:
pdf('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/4_FIGURES/CancerTsTvNewOld.pdf', height = 15, width = 15)
par(mfrow=c(2,2))

for (method in 1:3)
{# method = 1
  if (method == 1) 
    { 
    #### mtDB
    C1 = ALL[ALL$mtDB == 'notinDB',]; C1$Groups = 'new'; 
    C2 = ALL[ALL$mtDB != 'notinDB',]; C2$Groups = 'old';
    title = paste('Not in mtDB: ',nrow(C1),' ; In mtDB: ', nrow(C2),sep = '')
    }
  if (method == 2) 
    {
    #### Levin2012
    C1 = ALL[ALL$Levin2012 == 'notinLevin',]; C1$Groups = 'new';
    C2 = ALL[ALL$Levin2012 != 'notinLevin',]; C2$Groups = 'old';
    title = paste('Not in Levin: ',nrow(C1),' ; In Levin: ', nrow(C2),sep = '')
    }
  if (method == 3)   
    {
    #### C$NormalVarFreq
    C1 = ALL[ALL$NormalVarFreq == 0,]; #  'new';
    C2 = ALL[ALL$NormalVarFreq >  0,]; #  'old';
    title = paste('Normal Var Freq == 0: ',nrow(C1),' ; Normal Var Freq > 0: ', nrow(C2),sep = '')
    }
  
#### somatic_p_value - read about it???
#C1 = C[C$somatic_p_value == 0,] # 1793;
#C2 = C[C$somatic_p_value >  0,] # 5818;
#cor.test(C$somatic_p_value,C$NormalVarFreq, method = 'spearman')

  C1$Groups = C1$Tier2;
  C2$Groups = C2$Tier2;

  VecOfTransitionSubstitutions = c('A_G','G_A','C_T','T_C') # all tr
  VecOfTransversionSubstitutions = c('C_A','A_C','C_G','G_C','G_T','T_G','T_A','A_T') # ALL Tv 

  for (groups in 1:2)
{ # groups = 1
  if (groups == 1) {C = C1}
  if (groups == 2) {C = C2}

  MUT_TR = C[C$Subs %in% VecOfTransitionSubstitutions,]
  MUT_TV = C[C$Subs %in% VecOfTransversionSubstitutions,]

  AGG1 = aggregate(MUT_TR$NumberOfMut, by = list(MUT_TR$Groups), FUN = sum); names(AGG1) = c('Groups','Tr')
  AGG2 = aggregate(MUT_TV$NumberOfMut, by = list(MUT_TV$Groups), FUN = sum); names(AGG2) = c('Groups','Tv')
  TrTv = merge(AGG1,AGG2, by = ('Groups'))
  TrTv$TrTv = TrTv$Tr/TrTv$Tv
  if (groups == 1) {TrTv1 = TrTv}
  if (groups == 2) {TrTv2 = TrTv}
}
names(TrTv1) =c('Groups','TrNew','TsNew','TsTvNew')
names(TrTv2) =c('Groups','TrOld','TsOld','TsTvOld')
TrTv = merge(TrTv1,TrTv2)
MINIMUM =  min(min(TrTv$TsTvNew),min(TrTv$TsTvOld))
MAXIMUM =  max(max(TrTv$TsTvNew),max(TrTv$TsTvOld))
plot(TrTv$TsTvNew,TrTv$TsTvOld, xlim=c(MINIMUM,MAXIMUM),ylim=c(MINIMUM,MAXIMUM), pch = '', main = title)
text(TrTv$TsTvNew,TrTv$TsTvOld, TrTv$Groups, xlim=c(MINIMUM,MAXIMUM),ylim=c(MINIMUM,MAXIMUM))
abline(a = 0, b = 1, col = 'red')
cor.test(TrTv$TsTvNew,TrTv$TsTvOld,method = 'spearman') # p = 0.055, rho = 0.49
}
dev.off();

###################################
###### 18.03.2018: Cancer number of divisions (from Table S1 one before last column Tomasetti and Vogelstein 2015 science) and T>C fraction (from bioarchive Campbell 2018)
###################################

rm(list=ls(all=TRUE))
ALL = read.table('/media/konstantinpopadin/ac45df81-e084-4d30-9653-5c57cc9b58fd/konstantinpopadin/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/CancerMtDnaSomaticMutations/mtDNA_snv_Oct2016.txt', head = TRUE, sep = '\t')

head(ALL)
ALL$TissueTier2 = paste(ALL$tissue,ALL$Tier2, sep = '_')
table(ALL$tissue)
table(ALL$Tier2)
table(ALL$TissueTier2)
ALL$Subs = paste(ALL$ref,ALL$var, sep = '_')

pdf('/media/konstantinpopadin/ac45df81-e084-4d30-9653-5c57cc9b58fd/konstantinpopadin/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/4_FIGURES/CancerTsTvCellDivisions.pdf', height = 15, width = 15)
par(mfcol = c(2,3))
for (methods in 1:3)
{ # methods = 1
  if (methods == 1)  {C = ALL; title = paste('ALL: ',nrow(C),sep = '') }
  if (methods == 2)  {C = ALL[ALL$mtDB == 'notinDB',]; title = paste('NEW: ',nrow(C),sep = '') } # new only (5609)
  if (methods == 3)  {C = ALL[ALL$mtDB != 'notinDB',];  title = paste('OLD: ',nrow(C),sep = '')} # old only (2002)
  
  table(C$Subs)
  C$NumberOfMut = 1

VecOfTransitionSubstitutions = c('T_C')  # c('A_G','G_A','C_T','T_C') # all tr
VecOfTransversionSubstitutions = c('C_A','A_C','C_G','G_C','G_T','T_G','T_A','A_T') # ALL Tv 

MUT_TR = C[C$Subs %in% VecOfTransitionSubstitutions,]
MUT_TV = C[C$Subs %in% VecOfTransversionSubstitutions,]

AGG1 = aggregate(MUT_TR$NumberOfMut, by = list(MUT_TR$Tier2), FUN = sum); names(AGG1) = c('CancerTypes','Tr')
AGG2 = aggregate(MUT_TV$NumberOfMut, by = list(MUT_TV$Tier2), FUN = sum); names(AGG2) = c('CancerTypes','Tv')
TrTv = merge(AGG1,AGG2, by = 'CancerTypes')
TrTv$TrTv = TrTv$Tr/TrTv$Tv
TrTv = TrTv[order(TrTv$TrTv),]
TrTv
VecOfTissues = TrTv$CancerTypes 
# Kidney          Myeloid         Bladder         CNS             Stomach         Esophagus       Colon/Rectum    Lymphoid        Head/Neck       Cervix          Breast         
# Pancreas        Liver           Biliary         Thyroid         Bone/SoftTissue Ovary           Uterus          Skin            Lung            Prostate

Cells = data.frame(matrix( c("Thyroid","Head/Neck","Ovary","Esophagus","Pancreas","Skin","Lung","Liver","Colon/Rectum","Lymphoid","CNS","Bone/SoftTissue","Myeloid",
                     7,1720,0,1390,80,199,5.6,88,5840,960,0,5,960), ncol = 2))
names(Cells) = c('CancerTypes','NumOfCellDivPerLife')
Cells$NumOfCellDivPerLife = as.numeric(as.character(Cells$NumOfCellDivPerLife))
str(Cells)
TrTv = merge(TrTv,Cells)
TrTv = TrTv[order(TrTv$TrTv),]
str(TrTv)

plot(TrTv$NumOfCellDivPerLife,TrTv$TrTv, pch = '', ylab = 'Ts/Tv', xlab = 'Number Of Divisoins Of Each Stemm Cell Per Lifetime', main = title, xlim = c(-500, 6500)); 
cor.test(TrTv$TrTv,TrTv$NumOfCellDivPerLife, method = 'spearman') # -0.6804434, p-value = 0.01048
text(TrTv$NumOfCellDivPerLife,TrTv$TrTv,TrTv$CancerTypes)
boxplot(TrTv[TrTv$TrTv < median(TrTv$TrTv),]$NumOfCellDivPerLife,TrTv[TrTv$TrTv > median(TrTv$TrTv),]$NumOfCellDivPerLife,names =c('LowTsTv','HighTsTv'), ylab = 'Number Of Divisoins Of Each Stemm Cell Per Lifetime')
}
dev.off()

###################################
###### 16.03.2018: Parsimony data
###################################

rm(list=ls(all=TRUE))
GT = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/GenerationLenghtforAllMammals/GenerationLenghtforMammals.xlsx.txt', header = TRUE, sep = '\t')
head(GT)
str(GT)
GT$Species = gsub(' ','_',GT$Scientific_name)
length(unique(GT$Species))
summary(GT$AdultBodyMass_g)     # 2        21        71    136058       614 154321304 
summary(GT$GenerationLength_d)  # 129.0   624.4  1101.1  1578.8  2064.3 18980.0 # max = 18980 days => 52 years. ok
GT = GT[,c(11,13)]
summary(GT$GenerationLength_d)

############ MUTATION SPECTRUM 
MUT = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/2_DERIVED/Table_fixed_parsymony.txt', header = TRUE)
nrow(MUT)
VecOfNormalSubstitutions = c('A_C','C_A','A_G','G_A','C_G','G_C','C_T','T_C','G_T','T_G','T_A','A_T')

##### FILTER: normal Substitutions
MUT = MUT[MUT$Subs %in% VecOfNormalSubstitutions,]
nrow(MUT)

##### FILTER: Synonymous Substitutions
MUT = MUT[MUT$AncestralAA == MUT$DescendantAA,]
nrow(MUT)

##### FILTER: Gene
table(MUT$Gene)
MUT = MUT[MUT$Gene == 'CytB',]

##### NORMALIZATION: the third position of four-fold synonymous substitutions:
NUC = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/2_DERIVED/ATGC_counts_in_SYN_codons_wit_full_gene.txt', header = TRUE)
NUC$Gene = gsub("(.*)\\.",'',NUC$Species)
NUC$Species = gsub("\\.(.*)",'',NUC$Species)
MUT = merge(MUT,NUC, by = c("Species","Gene"))  # compare CountA.x and CountA.y  - they should be identical.
nrow(MUT) # A bit less !!!! WHY????

XTRACT = function(x) {first = unlist(strsplit(as.character(x),'_'))[1]; return(first);}; MUT$AncestralNuc = apply(as.matrix(MUT$Subs), 1, EXTRACT)
MUT$NumberOfSynMutPerSpecies = 1
MUT_A = MUT[MUT$AncestralNuc == 'A',]; MUT_T = MUT[MUT$AncestralNuc == 'T',]; MUT_G = MUT[MUT$AncestralNuc == 'G',]; MUT_C = MUT[MUT$AncestralNuc == 'C',] # 64145+123587+97657+128195=413584 
MUT_A$NumberOfSynMutPerSpecies = MUT_A$NumberOfSynMutPerSpecies/MUT_A$CountA_Syn;
MUT_T$NumberOfSynMutPerSpecies = MUT_T$NumberOfSynMutPerSpecies/MUT_T$CountT_Syn;
MUT_G$NumberOfSynMutPerSpecies = MUT_G$NumberOfSynMutPerSpecies/MUT_G$CountG_Syn;
MUT_C$NumberOfSynMutPerSpecies = MUT_C$NumberOfSynMutPerSpecies/MUT_C$CountC_Syn;
MUT = rbind(MUT_A,MUT_T,MUT_G,MUT_C)

##### FILTER: fourfold degenerate sites:
VecOfSynFourFoldDegenerateSites = c('CTT','CTC','CTA','CTG','GTT','GTC','GTA','GTG','TCT','TCC','TCA','TCG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','GCT','GCC','GCA','GCG','CGT','CGC','CGA','CGG','GGT','GGC','GGA','GGG'); length(unique(VecOfSynFourFoldDegenerateSites)) # 32
MUT = MUT[MUT$AncestorCodon %in% VecOfSynFourFoldDegenerateSites & MUT$DescendantCodon %in% VecOfSynFourFoldDegenerateSites,]
nrow(MUT)

######## CALCULATE MATRIX
### CALCULATE MATRIX 1: TAKE ONLY SPECIES WITH MANY SUBSTITUTIONS:
MUT$NumberOfAnyMutPerSpecies = 1
AGG1 = aggregate(MUT$NumberOfAnyMutPerSpecies, by = list(MUT$Species), FUN = sum)
summary(AGG1$x) # 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.0    30.0   100.0   239.5   285.0  2629.0 
ListOfSpeciesWithManySubst = AGG1[AGG1$x >= 20,]$Group.1; length(ListOfSpeciesWithManySubst) # 232
MUT = MUT[MUT$Species %in% ListOfSpeciesWithManySubst,]
nrow(MUT)

###### Tr and Tv
VecOfTransitionSubstitutions = c('A_G','G_A','C_T','T_C') # all tr
VecOfTransversionSubstitutions = c('C_A','A_C','C_G','G_C','G_T','T_G','T_A','A_T') # ALL Tv 

MUT_TR = MUT[MUT$Subs %in% VecOfTransitionSubstitutions,]
MUT_TV = MUT[MUT$Subs %in% VecOfTransversionSubstitutions,]

AGG1 = aggregate(MUT_TR$NumberOfSynMutPerSpecies, by = list(MUT_TR$Species), FUN = sum); names(AGG1) = c('Species','Tr')
AGG2 = aggregate(MUT_TV$NumberOfSynMutPerSpecies, by = list(MUT_TV$Species), FUN = sum); names(AGG2) = c('Species','Tv')
TrTv = merge(AGG1,AGG2, by = 'Species')
TrTv = merge(TrTv,GT,by = 'Species'); summary(TrTv$GenerationLength_d)
TrTv$TrTv = TrTv$Tr/TrTv$Tv; summary(TrTv$TrTv); TrTv = TrTv[!is.na(TrTv$TrTv),]

hist(TrTv$TrTv, breaks = 100)
cor.test(TrTv$GenerationLength_d,TrTv$TrTv, method = 'spearman')
a<-lm(scale(TrTv$TrTv) ~ scale(TrTv$GenerationLength_d)); summary(a) # all genes: 1.668e-01; CYTB: 2.666e-01; COX1: 2.785e-01;
plot(log2(TrTv$GenerationLength_d),log2(TrTv$TrTv))

###################################
###### 08.03.2018: From GC and to GC for mammals
###################################

############# GENERATION LENGTH FOR ALL MAMMALS
#  https://natureconservation.pensoft.net/article/1343/  data for generation length of ALL mammals!!!
rm(list=ls(all=TRUE))
GT = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/GenerationLenghtforAllMammals/GenerationLenghtforMammals.xlsx.txt', header = TRUE, sep = '\t')
head(GT)
str(GT)
GT$Species = gsub(' ','_',GT$Scientific_name)
length(unique(GT$Species))
summary(GT$AdultBodyMass_g)     # 2        21        71    136058       614 154321304 
summary(GT$GenerationLength_d)  # 129.0   624.4  1101.1  1578.8  2064.3 18980.0 # max = 18980 days => 52 years. ok
GT = GT[,c(11,13)]
summary(GT$GenerationLength_d)

############ MUTATION SPECTRUM 
MUT = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/2_DERIVED/Table_fixed.txt', header = TRUE)
nrow(MUT)
VecOfNormalSubstitutions = c('A_C','C_A','A_G','G_A','C_G','G_C','C_T','T_C','G_T','T_G','T_A','A_T')

##### FILTER: normal Substitutions
MUT = MUT[MUT$Subs %in% VecOfNormalSubstitutions,]
nrow(MUT)

##### FILTER: Synonymous Substitutions
MUT = MUT[MUT$AncestralAA == MUT$DescendantAA,]
nrow(MUT)

##### FILTER: Gene
table(MUT$Gene)
MUT = MUT[MUT$Gene == 'CytB',]

##### NORMALIZATION: the third position of four-fold synonymous substitutions:
NUC = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/2_DERIVED/ATGC_counts_in_SYN_codons_wit_full_gene.txt', header = TRUE)
NUC$Gene = gsub("(.*)\\.",'',NUC$Species)
NUC$Species = gsub("\\.(.*)",'',NUC$Species)
MUT = merge(MUT,NUC, by = c("Species","Gene"))  # compare CountA.x and CountA.y  - they should be identical.
nrow(MUT) # A bit less !!!! WHY????

EXTRACT = function(x) {first = unlist(strsplit(as.character(x),'_'))[1]; return(first);}; MUT$AncestralNuc = apply(as.matrix(MUT$Subs), 1, EXTRACT)
MUT$NumberOfSynMutPerSpecies = 1
MUT_A = MUT[MUT$AncestralNuc == 'A',]; MUT_T = MUT[MUT$AncestralNuc == 'T',]; MUT_G = MUT[MUT$AncestralNuc == 'G',]; MUT_C = MUT[MUT$AncestralNuc == 'C',] # 64145+123587+97657+128195=413584 
MUT_A$NumberOfSynMutPerSpecies = MUT_A$NumberOfSynMutPerSpecies/MUT_A$CountA_Syn;
MUT_T$NumberOfSynMutPerSpecies = MUT_T$NumberOfSynMutPerSpecies/MUT_T$CountT_Syn;
MUT_G$NumberOfSynMutPerSpecies = MUT_G$NumberOfSynMutPerSpecies/MUT_G$CountG_Syn;
MUT_C$NumberOfSynMutPerSpecies = MUT_C$NumberOfSynMutPerSpecies/MUT_C$CountC_Syn;
MUT = rbind(MUT_A,MUT_T,MUT_G,MUT_C)

##### FILTER: fourfold degenerate sites:
VecOfSynFourFoldDegenerateSites = c('CTT','CTC','CTA','CTG','GTT','GTC','GTA','GTG','TCT','TCC','TCA','TCG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','GCT','GCC','GCA','GCG','CGT','CGC','CGA','CGG','GGT','GGC','GGA','GGG'); length(unique(VecOfSynFourFoldDegenerateSites)) # 32
MUT = MUT[MUT$AncestorCodon %in% VecOfSynFourFoldDegenerateSites & MUT$DescendantCodon %in% VecOfSynFourFoldDegenerateSites,]
nrow(MUT)

######## CALCULATE MATRIX
### CALCULATE MATRIX 1: TAKE ONLY SPECIES WITH MANY SUBSTITUTIONS:
MUT$NumberOfAnyMutPerSpecies = 1
AGG1 = aggregate(MUT$NumberOfAnyMutPerSpecies, by = list(MUT$Species), FUN = sum)
summary(AGG1$x) # 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.0    30.0   100.0   239.5   285.0  2629.0 
ListOfSpeciesWithManySubst = AGG1[AGG1$x >= 30,]$Group.1; length(ListOfSpeciesWithManySubst) # 232
MUT = MUT[MUT$Species %in% ListOfSpeciesWithManySubst,]
nrow(MUT)

###### Tr and Tv

VecOfTransitionToCG = c('A_G','T_C') 
VecOfTransitionFromCG = c('G_A','C_T') 
VecOfTransversionToCG = c('A_C','T_G') 
VecOfTransversionFromCG = c('C_A','G_T') 
VecOfTransversionCGNeutral = c('C_G','G_C','T_A','A_T')
VecAllToCG = c('A_G','T_C','A_C','T_G') 
VecAllFromCG = c('G_A','C_T','C_A','G_T')
VecOfCommonTransitionToCG = c('T_C')
VecOfCommonTransitionFromCG = c('G_A') 

#MUT_1 = MUT[MUT$Subs %in% VecAllToCG,]
#MUT_2 = MUT[MUT$Subs %in% VecAllFromCG,]

MUT_1 = MUT[MUT$Subs %in% VecOfTransitionToCG,]
MUT_2 = MUT[MUT$Subs %in% VecOfTransitionFromCG,]

#MUT_1 = MUT[MUT$Subs %in% VecOfTransversionToCG,]
#MUT_2 = MUT[MUT$Subs %in% VecOfTransversionFromCG,]

#MUT_1 = MUT[MUT$Subs %in% VecAllToCG,]
#MUT_2 = MUT[MUT$Subs %in% VecOfTransversionCGNeutral,]

#MUT_1 = MUT[MUT$Subs %in% VecAllFromCG,]
#MUT_2 = MUT[MUT$Subs %in% VecOfTransversionCGNeutral,]

#MUT_1 = MUT[MUT$Subs %in% VecOfTransversionToCG,]
#MUT_2 = MUT[MUT$Subs %in% VecOfTransversionCGNeutral,]

#MUT_1 = MUT[MUT$Subs %in% VecOfTransversionFromCG,]
#MUT_2 = MUT[MUT$Subs %in% VecOfTransversionCGNeutral,]

#MUT_1 = MUT[MUT$Subs %in% VecOfCommonTransitionToCG,]
#MUT_2 = MUT[MUT$Subs %in% VecOfCommonTransitionFromCG,]

AGG1 = aggregate(MUT_1$NumberOfSynMutPerSpecies, by = list(MUT_1$Species), FUN = sum); names(AGG1) = c('Species','Type1')
AGG2 = aggregate(MUT_2$NumberOfSynMutPerSpecies, by = list(MUT_2$Species), FUN = sum); names(AGG2) = c('Species','Type2')
TrTv = merge(AGG1,AGG2, by = 'Species')
TrTv = merge(TrTv,GT,by = 'Species'); summary(TrTv$GenerationLength_d)
TrTv$Ratio = TrTv$Type1/TrTv$Type2; summary(TrTv$Ratio); TrTv = TrTv[!is.na(TrTv$Ratio),]

hist(TrTv$Ratio, breaks = 100)
cor.test(TrTv$GenerationLength_d,TrTv$Ratio, method = 'spearman')
a<-lm(scale(TrTv$Ratio) ~ 0+scale(TrTv$GenerationLength_d)); summary(a) # all genes: 1.668e-01; CYTB: 2.666e-01; COX1: 2.785e-01;
plot(log2(TrTv$GenerationLength_d),log2(TrTv$Ratio))

# VecAllToCG/VecAllFromCG: rho = 0.1666929; p = 0.0005767
# VecOfTransitionToCG/VecOfTransitionFromCG: rho = 0.1693457; p = 0.0004689
# VecOfTransversionToCG/VecOfTransversionFromCG: rho = 0.1598554; p = 0.003699
# VecAllToCG/VecOfTransversionCGNeutral: rho = 0.29208; p = 2.184e-09
# VecAllFromCG/VecOfTransversionCGNeutral: rho = 0.2026437; p = 3.57e-05 
# VecOfTransversionToCG/VecOfTransversionCGNeutral: rho = 0.03320413; p = 0.5472
# VecOfTransversionFromCG/VecOfTransversionCGNeutral: rho = -0.1825543; p = 0.0003352

###################################
###### 04.04.2018: correlation between substitutions, PCA
###################################

rm(list=ls(all=TRUE))

############# GENERATION LENGTH FOR ALL MAMMALS
#  https://natureconservation.pensoft.net/article/1343/  data for generation length of ALL mammals!!!
GenerTtime = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/GenerationLenghtforAllMammals/GenerationLenghtforMammals.xlsx.txt', header = TRUE, sep = '\t')
head(GenerTtime)
str(GenerTtime)
GenerTtime$Species = gsub(' ','_',GenerTtime$Scientific_name)
length(unique(GenerTtime$Species))
summary(GenerTtime$AdultBodyMass_g)     # 2        21        71    136058       614 154321304 
summary(GenerTtime$GenerationLength_d)  # 129.0   624.4  1101.1  1578.8  2064.3 18980.0 # max = 18980 days => 52 years. ok
GenerTime = GenerTtime[,c(11,13)]
summary(GenerTime$GenerationLength_d)
ListOfMammals = GenerTime$Species; length(ListOfMammals) # 5426 - at least will use it to work with mammals only

# AnAge
AnAge = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/1_RAW_DATA/anage_data.txt', header = TRUE, sep = '\t')
head(AnAge)
AnAge = AnAge[AnAge$Class == 'Mammalia' & ! is.na(AnAge$Metabolic.rate..W.) & ! is.na(AnAge$Adult.weight..g.),];  #  & ! is.na(AnAge$R)
plot(log2(AnAge$Metabolic.rate..W.), log2(AnAge$Adult.weight..g.))
AnAge$BodyMassNormalizedBmr = log2(AnAge$Metabolic.rate..W.)/log2(AnAge$Adult.weight..g.)
AnAge$Species = paste(AnAge$Genus,AnAge$Species, sep='_') 
a = lm(log2(AnAge$Metabolic.rate..W.)~log2(AnAge$Adult.weight..g.))
b = lm(log2(AnAge$Adult.weight..g.)~log2(AnAge$Metabolic.rate..W.))
AnAge$Residuals = a$residuals
summary(AnAge$Residuals)
TooColdSpecies = AnAge[AnAge$Residuals < quantile(AnAge$Residuals,0.25),]$Species; length(TooColdSpecies)
TooHotSpecies = AnAge[AnAge$Residuals > quantile(AnAge$Residuals,0.75),]$Species; length(TooHotSpecies)
AnAge1 = data.frame(AnAge$Species,AnAge$BodyMassNormalizedBmr,AnAge$Residuals,log2(AnAge$Metabolic.rate..W.),log2(AnAge$Adult.weight..g.)); names(AnAge1)=c('Species','BodyMassNormalizedBmr','Residuals','log2MR','log2Weight')

AnAge = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/1_RAW_DATA/anage_data.txt', header = TRUE, sep = '\t')
head(AnAge)
AnAge = AnAge[AnAge$Class == 'Mammalia' & ! is.na(AnAge$Temperature..K.),]; 
AnAge$Species = paste(AnAge$Genus,AnAge$Species, sep='_') 
AnAge2 = data.frame(AnAge$Species,AnAge$Temperature..K.); names(AnAge2)=c('Species','Temperature')

AnAge = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/1_RAW_DATA/anage_data.txt', header = TRUE, sep = '\t')
head(AnAge)
AnAge = AnAge[AnAge$Class == 'Mammalia' & ! is.na(AnAge$Maximum.longevity..yrs.),]; 
AnAge$Species = paste(AnAge$Genus,AnAge$Species, sep='_') 
AnAge3 = data.frame(AnAge$Species,AnAge$Maximum.longevity..yrs.); names(AnAge3)=c('Species','MLS')

AnAge = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/1_RAW_DATA/anage_data.txt', header = TRUE, sep = '\t')
head(AnAge)
AnAge = AnAge[AnAge$Class == 'Mammalia' & ! is.na(AnAge$Litter.Clutch.size),]; 
AnAge$Species = paste(AnAge$Genus,AnAge$Species, sep='_') 
AnAge4 = data.frame(AnAge$Species,AnAge$Litter.Clutch.size); names(AnAge4)=c('Species','ClutchSize')

AnAge = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/1_RAW_DATA/anage_data.txt', header = TRUE, sep = '\t')
head(AnAge)
AnAge = AnAge[AnAge$Class == 'Mammalia' & ! is.na(AnAge$Birth.weight..g.),]; 
AnAge$Species = paste(AnAge$Genus,AnAge$Species, sep='_') 
AnAge5 = data.frame(AnAge$Species,AnAge$Birth.weight..g.); names(AnAge5)=c('Species','BirthWeight')

AnAge = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/1_RAW_DATA/anage_data.txt', header = TRUE, sep = '\t')
head(AnAge)
AnAge = AnAge[AnAge$Class == 'Mammalia' & ! is.na(AnAge$Maximum.longevity..yrs.) & ! is.na(AnAge$Adult.weight..g.),];
AnAge$Species = paste(AnAge$Genus,AnAge$Species, sep='_') 
plot(log2(AnAge$Adult.weight..g.),log2(AnAge$Maximum.longevity..yrs.))
a = lm(log2(AnAge$Maximum.longevity..yrs.) ~ log2(AnAge$Adult.weight..g.))
b = lm(log2(AnAge$Adult.weight..g.) ~ log2(AnAge$Maximum.longevity..yrs.))
AnAge$Residuals = a$residuals
AnAge6 = data.frame(AnAge$Species,log2(AnAge$Adult.weight..g.),log2(AnAge$Maximum.longevity..yrs.),AnAge$Residuals); names(AnAge6)=c('Species','Log2Weight','log2MLS','Residuals')

AnAge = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/1_RAW_DATA/anage_data.txt', header = TRUE, sep = '\t')
head(AnAge)
AnAge = AnAge[AnAge$Class == 'Mammalia' & ! is.na(AnAge$Birth.weight..g.) & ! is.na(AnAge$Adult.weight..g.),];
AnAge$Species = paste(AnAge$Genus,AnAge$Species, sep='_')
AnAge$PropagulaRelativeSize = AnAge$Birth.weight..g./AnAge$Adult.weight..g.; summary(AnAge$PropagulaRelativeSize)
AnAge7 = data.frame(AnAge$Species,AnAge$PropagulaRelativeSize); names(AnAge7)=c('Species','PropagulaRelativeSize')

AnAge = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/1_RAW_DATA/anage_data.txt', header = TRUE, sep = '\t')
head(AnAge)
AnAge = AnAge[AnAge$Class == 'Mammalia' & ! is.na(AnAge$Growth.rate..1.days.),];
AnAge$Species = paste(AnAge$Genus,AnAge$Species, sep='_')
AnAge8 = data.frame(AnAge$Species,AnAge$Growth.rate..1.days.); names(AnAge8)=c('Species','Growth.rate..1.days.')

# hibernation
Hib = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/hibernation/tornew5.csv', header = TRUE, sep = ',')
Hib$Species = gsub(" ",'_',Hib$Species)
HibOnlySpecies = Hib[Hib$Type == 'HIB',]$Species; length(HibOnlySpecies)
DtOnlySpecies = Hib[Hib$Type == 'DT',]$Species; length(DtOnlySpecies)
HibAndDtSpecies = Hib[Hib$Type == 'HIB' || Hib$Type == 'DT',]$Species; length(HibAndDtSpecies)

# little (30 g) but not hibernating or DT:
AnAge = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/1_RAW_DATA/anage_data.txt', header = TRUE, sep = '\t')
head(AnAge)
AnAge = AnAge[AnAge$Class == 'Mammalia' & ! is.na(AnAge$Adult.weight..g.),]; 
AnAge$Species = paste(AnAge$Genus,AnAge$Species, sep='_')
summary(AnAge[AnAge$Species %in% HibAndDtSpecies,]$Adult.weight..g.) # 2.1     15.0     32.5   4701.6    175.0 277500.0
summary(AnAge[AnAge$Species %in% DtOnlySpecies,]$Adult.weight..g.) # 2.10    11.45    25.25   499.57    54.55 10000.00
summary(AnAge[AnAge$Species %in% HibOnlySpecies,]$Adult.weight..g.) # 4.6     21.0     48.2   8824.3    380.0 277500.0
AnAge = AnAge[!AnAge$Species %in% HibAndDtSpecies,]; 
LittleSpeciesButNotHibOrDt = AnAge[AnAge$Adult.weight..g. < 175,]$Species; length(LittleSpeciesButNotHibOrDt)

AnAge = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/1_RAW_DATA/anage_data.txt', header = TRUE, sep = '\t')
head(AnAge)
AnAge = AnAge[AnAge$Class == 'Mammalia',];
AnAge$Species = paste(AnAge$Genus,AnAge$Species, sep='_')
Marsupials = c(AnAge[AnAge$Order == 'Diprotodontia' | AnAge$Order == 'Didelphimorphia' | AnAge$Order == 'Dasyuromorphia',]$Species)
Placental = c(AnAge[AnAge$Order == 'Artiodactyla' | AnAge$Order == 'Cetacea' | AnAge$Order == 'Carnivora',]$Species)


############ MUTATION SPECTRUM 
# METHOD = 'PARSIMONY'
METHOD = 'MAXLIKELIHOOD'

if (METHOD == 'MAXLIKELIHOOD') {MUT = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/2_DERIVED/Table_fixed.txt', header = TRUE)}
if (METHOD == 'PARSIMONY')     {MUT = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/2_DERIVED/Table_fixed_parsymony.txt', header = TRUE)}
nrow(MUT)
VecOfNormalSubstitutions = c('A_C','C_A','A_G','G_A','C_G','G_C','C_T','T_C','G_T','T_G','T_A','A_T')

##### FILTER: normal Substitutions
MUT = MUT[MUT$Subs %in% VecOfNormalSubstitutions,]
nrow(MUT)

##### FILTER: Synonymous Substitutions
MUT = MUT[MUT$AncestralAA == MUT$DescendantAA,]
nrow(MUT)

##### FILTER: Gene
table(MUT$Gene)
MUT = MUT[MUT$Gene == 'CytB',]
# MUT = MUT[MUT$Gene == 'COX1' | MUT$Gene == 'COX2' | MUT$Gene == 'COX3',]

##### NORMALIZATION by nucleotide count in the third position of four-fold synonymous substitutions:
NUC = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/2_DERIVED/ATGC_counts_in_SYN_codons_wit_full_gene.txt', header = TRUE)
NUC$Gene = gsub("(.*)\\.",'',NUC$Species)
NUC$Species = gsub("\\.(.*)",'',NUC$Species)
MUT = merge(MUT,NUC, by = c("Species","Gene"))  # compare CountA.x and CountA.y  - they should be identical.
nrow(MUT) # A bit less !!!! WHY????

EXTRACT = function(x) {first = unlist(strsplit(as.character(x),'_'))[1]; return(first);}; MUT$AncestralNuc = apply(as.matrix(MUT$Subs), 1, EXTRACT)
MUT$NumberOfSynMutPerSpecies = 1
MUT_A = MUT[MUT$AncestralNuc == 'A',]; MUT_T = MUT[MUT$AncestralNuc == 'T',]; MUT_G = MUT[MUT$AncestralNuc == 'G',]; MUT_C = MUT[MUT$AncestralNuc == 'C',] # 64145+123587+97657+128195=413584 
MUT_A$NumberOfSynMutPerSpecies = MUT_A$NumberOfSynMutPerSpecies/MUT_A$CountA_Syn;
MUT_T$NumberOfSynMutPerSpecies = MUT_T$NumberOfSynMutPerSpecies/MUT_T$CountT_Syn;
MUT_G$NumberOfSynMutPerSpecies = MUT_G$NumberOfSynMutPerSpecies/MUT_G$CountG_Syn;
MUT_C$NumberOfSynMutPerSpecies = MUT_C$NumberOfSynMutPerSpecies/MUT_C$CountC_Syn;
MUT = rbind(MUT_A,MUT_T,MUT_G,MUT_C)

##### FILTER: fourfold degenerate sites:
VecOfSynFourFoldDegenerateSites = c('CTT','CTC','CTA','CTG','GTT','GTC','GTA','GTG','TCT','TCC','TCA','TCG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','GCT','GCC','GCA','GCG','CGT','CGC','CGA','CGG','GGT','GGC','GGA','GGG'); length(unique(VecOfSynFourFoldDegenerateSites)) # 32
MUT = MUT[MUT$AncestorCodon %in% VecOfSynFourFoldDegenerateSites & MUT$DescendantCodon %in% VecOfSynFourFoldDegenerateSites,]
nrow(MUT)

### VECTOR OF SPECIES WITH MANY SUBSTITUTIONS:
MUT$NumberOfAnyMutPerSpecies = 1
AGG1 = aggregate(MUT$NumberOfAnyMutPerSpecies, by = list(MUT$Species), FUN = sum)
summary(AGG1$x) # 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.0    30.0   100.0   239.5   285.0  2629.0 
ListOfSpeciesWithManySubst = AGG1[AGG1$x >= 30,]$Group.1; length(ListOfSpeciesWithManySubst) # 232

##### COUNT THE TOTAL NUMBER OF NORMALIZED MUTATIONS PER SPECIES

AggTotalMutSpectrumPerSpecies = aggregate(MUT$NumberOfSynMutPerSpecies, by = list(MUT$Species), FUN = sum); 
names(AggTotalMutSpectrumPerSpecies) = c('Species','TotalMutRate')
MutTypes = data.frame(VecOfNormalSubstitutions); names(MutTypes)=c('MutType')
nrow(AggTotalMutSpectrumPerSpecies) # 1710
AggTotalMutSpectrumPerSpecies = merge(AggTotalMutSpectrumPerSpecies,MutTypes)
nrow(AggTotalMutSpectrumPerSpecies) # 20520/12==1710!
AggTotalMutSpectrumPerSpeciesPerMutType = aggregate(MUT$NumberOfSynMutPerSpecies, by = list(MUT$Species,MUT$Subs), FUN = sum); 
names(AggTotalMutSpectrumPerSpeciesPerMutType) = c('Species','MutType','MutTypeRate')
nrow(AggTotalMutSpectrumPerSpeciesPerMutType) # 12357
ALL = merge(AggTotalMutSpectrumPerSpecies,AggTotalMutSpectrumPerSpeciesPerMutType, by = c('Species','MutType'), all.x=TRUE)
ALL[is.na(ALL)]<-0
ALL$Fraction = ALL$MutTypeRate/ALL$TotalMutRate; summary(ALL$Fraction)
ALL = ALL[ALL$Species %in% ListOfSpeciesWithManySubst,]

# check of these species have all 12 nonzero 12 types of substitutions:
ALL$number = 1
AGG = aggregate(ALL$number, by = list(ALL$Species), FUN = sum); summary(AGG$x)
VecOfSpeciesWith12Subst = AGG$Group.1

# create matrix for PCA:

AT = ALL[ALL$MutType == 'A_T',]; AT = AT[c(1,5)]; names(AT) = c('Species','AT'); AT = AT[order(AT$Species),]
AG = ALL[ALL$MutType == 'A_G',]; AG = AG[c(1,5)]; names(AG) = c('Species','AG'); AG = AG[order(AG$Species),]
AC = ALL[ALL$MutType == 'A_C',]; AC = AC[c(1,5)]; names(AC) = c('Species','AC'); AC = AC[order(AC$Species),]
TA = ALL[ALL$MutType == 'T_A',]; TA = TA[c(1,5)]; names(TA) = c('Species','TA'); TA = TA[order(TA$Species),]
TG = ALL[ALL$MutType == 'T_G',]; TG = TG[c(1,5)]; names(TG) = c('Species','TG'); TG = TG[order(TG$Species),]
TC = ALL[ALL$MutType == 'T_C',]; TC = TC[c(1,5)]; names(TC) = c('Species','TC'); TC = TC[order(TC$Species),]
CA = ALL[ALL$MutType == 'C_A',]; CA = CA[c(1,5)]; names(CA) = c('Species','CA'); CA = CA[order(CA$Species),]
CG = ALL[ALL$MutType == 'C_G',]; CG = CG[c(1,5)]; names(CG) = c('Species','CG'); CG = CG[order(CG$Species),]
CT = ALL[ALL$MutType == 'C_T',]; CT = CT[c(1,5)]; names(CT) = c('Species','CT'); CT = CT[order(CT$Species),]
GA = ALL[ALL$MutType == 'G_A',]; GA = GA[c(1,5)]; names(GA) = c('Species','GA'); GA = GA[order(GA$Species),]
GC = ALL[ALL$MutType == 'G_C',]; GC = GC[c(1,5)]; names(GC) = c('Species','GC'); GC = GC[order(GC$Species),]
GT = ALL[ALL$MutType == 'G_T',]; GT = GT[c(1,5)]; names(GT) = c('Species','GT'); GT = GT[order(GT$Species),]
MATRIX = cbind(AT,AG[,2],AC[,2],TA[,2],TG[,2],TC[,2],CA[,2],CG[,2],CT[,2],GA[,2],GC[,2],GT[,2]); names(MATRIX) = c('Species','AT','AG','AC','TA','TG','TC','CA','CG','CT','GA','GC','GT')
MATRIX = MATRIX[MATRIX$Species %in% ListOfMammals,] # 361

##### scale 12 substitutions. PCa1 strongly depends on GA - probably this is just because GA is very common and thus fluctuations are very important by absolute values - not by relative...
for (i in 2:13)
{ # i =2
  summary(MATRIX[,i])
  MATRIX[,i] = as.numeric(scale(MATRIX[,i]))
  summary(as.numeric(MATRIX[,i]))
}

###### PCA 
row.names(MATRIX)=MATRIX$Species
head(MATRIX)
matrix = MATRIX[,c(2:13)]
PCA = prcomp(matrix)
print(PCA) 
summary(PCA)
MATRIX$Pca1 = PCA$x[,1]
MATRIX$Pca2 = PCA$x[,2]
MATRIX$Pca3 = PCA$x[,3]

### list of species with extreme PC1
MATRIX[MATRIX$Pca1 < quantile(MATRIX$Pca1,0.05),]$Species 
# Ailuropoda_melanoleuca   Ateles_geoffroyi         Dasypus_novemcinctus     Dasyurus_hallucatus      Dryomys_nitedula         Genetta_servalina        Hylomys_suillus         
# Martes_flavigula         Niviventer_cremoriventer Nycticebus_coucang       Ochotona_rufescens       Oecomys_concolor         Oligoryzomys_fulvescens  Perodicticus_potto      
# Pteronotus_gymnonotus    Rhinophylla_pumilio      Tamias_striatus          Tamiops_swinhoei        

MATRIX[MATRIX$Pca1 > quantile(MATRIX$Pca1,0.95),]$Species  # Sorex, Hemiechinus...
# Apodemus_argenteus      Artibeus_lituratus      Galago_senegalensis     Handleyomys_rostratus   Hemiechinus_auritus     Leptonychotes_weddellii Macaca_mulatta         
# Metachirus_nudicaudatus Microcebus_murinus      Myodes_gapperi          Myotis_nattereri        Nycticebus_bengalensis  Plagiodontia_aedium     Rousettus_leschenaultii
# Sorex_araneus           Tamias_dorsalis         Tamias_minimus          Tursiops_truncatus     

### test 1: 

Test = merge(MATRIX,AnAge1) 
cor.test(Test$Pca1,Test$BodyMassNormalizedBmr, method = 'spearman')
cor.test(Test$Pca2,Test$BodyMassNormalizedBmr, method = 'spearman')
cor.test(Test$Pca3,Test$BodyMassNormalizedBmr, method = 'spearman')

cor.test(Test$Pca1,Test$Residuals, method = 'spearman')
cor.test(Test$Pca2,Test$Residuals, method = 'spearman')
cor.test(Test$Pca3,Test$Residuals, method = 'spearman')

cor.test(Test$Pca1,Test$log2MR, method = 'spearman')
cor.test(Test$Pca2,Test$log2MR, method = 'spearman')
cor.test(Test$Pca3,Test$log2MR, method = 'spearman')

cor.test(Test$Pca1,Test$log2Weight, method = 'spearman')
cor.test(Test$Pca2,Test$log2Weight, method = 'spearman')
cor.test(Test$Pca3,Test$log2Weight, method = 'spearman')

#### test 2:
Test = merge(MATRIX,GenerTime) 
cor.test(Test$Pca1,Test$GenerationLength_d, method = 'spearman') # a bit positive
cor.test(Test$Pca2,Test$GenerationLength_d, method = 'spearman') # super positive
cor.test(Test$Pca3,Test$GenerationLength_d, method = 'spearman') 

#### test 3: hibernating species have a bit less G>A! 
# control for body mass? they are small!
# other tests - decreased BMR as compared to body mass - compare strong outliers (cold and hot)...
# animals, which live more than should according to bosy mass (naked mole rat)

par(mfrow=c(1,2))
boxplot(MATRIX[MATRIX$Species %in% HibOnlySpecies,]$Pca1,MATRIX[MATRIX$Species %in% LittleSpeciesButNotHibOrDt,]$Pca1, notch = TRUE, names = c('HibOnlySpecies','LittleSpeciesButNotHibOrDt'), ylab = 'PC1');  # !!!! HibSpecies have a bit lower Pc1
wilcox.test(MATRIX[MATRIX$Species %in% HibOnlySpecies,]$Pca1,MATRIX[MATRIX$Species %in% LittleSpeciesButNotHibOrDt,]$Pca1)

boxplot(MATRIX[MATRIX$Species %in% TooColdSpecies,]$Pca1,MATRIX[MATRIX$Species %in% TooHotSpecies,]$Pca1, notch = TRUE, names = c('TooCold','TooHot'), ylab = 'PC1');  # TooCold have a bit lower PC1 than TooHot
wilcox.test(MATRIX[MATRIX$Species %in% TooColdSpecies,]$Pca1,MATRIX[MATRIX$Species %in% TooHotSpecies,]$Pca1, alternative = 'less') 

boxplot(MATRIX[MATRIX$Species %in% DtOnlySpecies,]$GA,MATRIX[MATRIX$Species %in% LittleSpeciesButNotHibOrDt,]$GA, notch = TRUE)

#### test 4:
Test = merge(MATRIX,AnAge2) 
cor.test(Test$Pca1,Test$Temperature, method = 'spearman')
cor.test(Test$Pca2,Test$Temperature, method = 'spearman')
cor.test(Test$Pca3,Test$Temperature, method = 'spearman')

#### test 5:
Test = merge(MATRIX,AnAge3) 
cor.test(Test$Pca1,Test$MLS, method = 'spearman')  # a bit positive
cor.test(Test$Pca2,Test$MLS, method = 'spearman')  # positive
cor.test(Test$Pca3,Test$MLS, method = 'spearman')  # a bit negative

#### test 6:
Test = merge(MATRIX,AnAge4) 
cor.test(Test$Pca1,Test$ClutchSize, method = 'spearman')  #
cor.test(Test$Pca2,Test$ClutchSize, method = 'spearman')  # negative
cor.test(Test$Pca3,Test$ClutchSize, method = 'spearman')  # 

#### test 7:
Test = merge(MATRIX,AnAge5) 
cor.test(Test$Pca1,Test$BirthWeight, method = 'spearman')  #
cor.test(Test$Pca2,Test$BirthWeight, method = 'spearman')  # positive
cor.test(Test$Pca3,Test$BirthWeight, method = 'spearman')  # 

#### test 8:
Test = merge(MATRIX,AnAge6) 
cor.test(Test$Pca1,Test$Residuals, method = 'spearman')  #
cor.test(Test$Pca2,Test$Residuals, method = 'spearman')  #
cor.test(Test$Pca3,Test$Residuals, method = 'spearman')  # 

#### test 9:
Test = merge(MATRIX,AnAge7) 
cor.test(Test$Pca1,Test$PropagulaRelativeSize, method = 'spearman')  #
cor.test(Test$Pca2,Test$PropagulaRelativeSize, method = 'spearman')  #
cor.test(Test$Pca3,Test$PropagulaRelativeSize, method = 'spearman')  # 

#### test 10:
Test = merge(MATRIX,AnAge8) 
cor.test(Test$Pca1,Test$Growth.rate..1.days., method = 'spearman')  #
cor.test(Test$Pca2,Test$Growth.rate..1.days., method = 'spearman')  # a bit negative
cor.test(Test$Pca3,Test$Growth.rate..1.days., method = 'spearman')  # 

#### test 11:
boxplot(MATRIX[MATRIX$Species %in% Marsupials,]$Pca1,MATRIX[MATRIX$Species %in% Placental,]$Pca1, notch = TRUE)
wilcox.test(MATRIX[MATRIX$Species %in% Marsupials,]$Pca1,MATRIX[MATRIX$Species %in% Placental,]$Pca1) # 0.09
boxplot(MATRIX[MATRIX$Species %in% Marsupials,]$Pca2,MATRIX[MATRIX$Species %in% Placental,]$Pca2, notch = TRUE)
wilcox.test(MATRIX[MATRIX$Species %in% Marsupials,]$Pca2,MATRIX[MATRIX$Species %in% Placental,]$Pca2) # significant
boxplot(MATRIX[MATRIX$Species %in% Marsupials,]$Pca3,MATRIX[MATRIX$Species %in% Placental,]$Pca3, notch = TRUE)

#### figures:

MATRIX = merge(MATRIX,GenerTime)
MATRIX = MATRIX[order(MATRIX$GenerationLength_d),]
MATRIX$Col = c(rep('green',100),rep('gray',161),rep('red',100))
head(MATRIX)
summary(PCA)
print(PCA)


pdf('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/4_FIGURES/PCA.pdf', width = 14, height = 14)
par(mfrow=c(2,3))
summary(PCA)
plot(PCA)
plot(MATRIX$Pca1,MATRIX$Pca2, col = MATRIX$Col)
plot(MATRIX$Pca2,MATRIX$Pca3, col = MATRIX$Col)
# plot(PCA$x[,1],MATRIX$GenerationLength_d); cor.test(PCA$x[,1],MATRIX$GenerationLength_d, method = 'spearman') # nothing  - First mutagen signature! Body mass normalized BMR!
plot(MATRIX$Pca2,log2(MATRIX$GenerationLength_d)); cor.test(MATRIX$Pca2,MATRIX$GenerationLength_d, method = 'spearman') 
biplot(PCA, col = c('grey','black'), cex = 0.5)
biplot(PCA, choices=c(2,3), col = c('grey','black'), cex = 0.5) #  biplot(princomp(USArrests),choices=c(1,3))
dev.off()


#### several analyses:
flag = 0
GoodOrder = c('T_C','G_A','A_G','C_T','A_C','C_G','G_C','T_A','T_G','A_T','C_A','G_T'); length(GoodOrder)
for (i in 1:length(GoodOrder))
{# i = 1
  Temp1 = ALL[ALL$MutType == GoodOrder[i],]; 
  Temp1 = Temp1[,c(1,5)]
  for (j in i:length(GoodOrder))
  {
    if (i != j)
    {
    Temp2 = ALL[ALL$MutType == GoodOrder[j],]
    Temp2 = Temp2[,c(1,5)]
    Temp = merge(Temp1,Temp2, by = 'Species')
    a = cor.test(Temp[,2],Temp[,3],method = 'spearman')
    Pval = as.numeric(a[3])
    Rho = as.numeric(a[4])
    OneLine = data.frame(GoodOrder[i],GoodOrder[j],Pval,Rho)
    if (flag == 1) {Final = rbind(Final,OneLine)}
    if (flag == 0) {Final = OneLine; flag = 1}
    }
  }
}
  
Final=Final[order(Final$Pval),]

### Princ Comp!!!

########################################
##### 6.04.2018: MUT SPECTRUM FOR DIFFERENT CLASSES
########################################

rm(list=ls(all=TRUE))

# AnAge
AnAge = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/1_RAW_DATA/anage_data.txt', header = TRUE, sep = '\t')
head(AnAge)
AnAge$Species = paste(AnAge$Genus,AnAge$Species, sep='_') 
AnAge = AnAge[,c(2:8)]

############ MUTATION SPECTRUM 
# METHOD = 'PARSIMONY'
METHOD = 'MAXLIKELIHOOD'

if (METHOD == 'MAXLIKELIHOOD') {MUT = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/2_DERIVED/Table_fixed.txt', header = TRUE)}
if (METHOD == 'PARSIMONY')     {MUT = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/2_DERIVED/Table_fixed_parsymony.txt', header = TRUE)}
nrow(MUT)
VecOfNormalSubstitutions = c('A_C','C_A','A_G','G_A','C_G','G_C','C_T','T_C','G_T','T_G','T_A','A_T')

##### FILTER: normal Substitutions
MUT = MUT[MUT$Subs %in% VecOfNormalSubstitutions,]
nrow(MUT)

##### FILTER: Synonymous Substitutions
MUT = MUT[MUT$AncestralAA == MUT$DescendantAA,]
nrow(MUT)

##### FILTER: Gene
table(MUT$Gene)
MUT = MUT[MUT$Gene == 'CytB',]
# MUT = MUT[MUT$Gene == 'COX1' | MUT$Gene == 'COX2' | MUT$Gene == 'COX3',]

##### NORMALIZATION by nucleotide count in the third position of four-fold synonymous substitutions:
NUC = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/2_DERIVED/ATGC_counts_in_SYN_codons_wit_full_gene.txt', header = TRUE)
NUC$Gene = gsub("(.*)\\.",'',NUC$Species)
NUC$Species = gsub("\\.(.*)",'',NUC$Species)
MUT = merge(MUT,NUC, by = c("Species","Gene"))  # compare CountA.x and CountA.y  - they should be identical.
nrow(MUT) # A bit less !!!! WHY????

EXTRACT = function(x) {first = unlist(strsplit(as.character(x),'_'))[1]; return(first);}; MUT$AncestralNuc = apply(as.matrix(MUT$Subs), 1, EXTRACT)
MUT$NumberOfSynMutPerSpecies = 1
MUT_A = MUT[MUT$AncestralNuc == 'A',]; MUT_T = MUT[MUT$AncestralNuc == 'T',]; MUT_G = MUT[MUT$AncestralNuc == 'G',]; MUT_C = MUT[MUT$AncestralNuc == 'C',] # 64145+123587+97657+128195=413584 
MUT_A$NumberOfSynMutPerSpecies = MUT_A$NumberOfSynMutPerSpecies/MUT_A$CountA_Syn;
MUT_T$NumberOfSynMutPerSpecies = MUT_T$NumberOfSynMutPerSpecies/MUT_T$CountT_Syn;
MUT_G$NumberOfSynMutPerSpecies = MUT_G$NumberOfSynMutPerSpecies/MUT_G$CountG_Syn;
MUT_C$NumberOfSynMutPerSpecies = MUT_C$NumberOfSynMutPerSpecies/MUT_C$CountC_Syn;
MUT = rbind(MUT_A,MUT_T,MUT_G,MUT_C)

##### FILTER: fourfold degenerate sites:
VecOfSynFourFoldDegenerateSites = c('CTT','CTC','CTA','CTG','GTT','GTC','GTA','GTG','TCT','TCC','TCA','TCG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','GCT','GCC','GCA','GCG','CGT','CGC','CGA','CGG','GGT','GGC','GGA','GGG'); length(unique(VecOfSynFourFoldDegenerateSites)) # 32
MUT = MUT[MUT$AncestorCodon %in% VecOfSynFourFoldDegenerateSites & MUT$DescendantCodon %in% VecOfSynFourFoldDegenerateSites,]
nrow(MUT)

### VECTOR OF SPECIES WITH MANY SUBSTITUTIONS:
MUT$NumberOfAnyMutPerSpecies = 1
AGG1 = aggregate(MUT$NumberOfAnyMutPerSpecies, by = list(MUT$Species), FUN = sum)
summary(AGG1$x) # 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.0    30.0   100.0   239.5   285.0  2629.0 
ListOfSpeciesWithManySubst = AGG1[AGG1$x >= 30,]$Group.1; length(ListOfSpeciesWithManySubst) # 232

##### COUNT THE TOTAL NUMBER OF NORMALIZED MUTATIONS PER SPECIES

AggTotalMutSpectrumPerSpecies = aggregate(MUT$NumberOfSynMutPerSpecies, by = list(MUT$Species), FUN = sum); 
names(AggTotalMutSpectrumPerSpecies) = c('Species','TotalMutRate')
MutTypes = data.frame(VecOfNormalSubstitutions); names(MutTypes)=c('MutType')
nrow(AggTotalMutSpectrumPerSpecies) # 1710
AggTotalMutSpectrumPerSpecies = merge(AggTotalMutSpectrumPerSpecies,MutTypes)
nrow(AggTotalMutSpectrumPerSpecies) # 20520/12==1710!
AggTotalMutSpectrumPerSpeciesPerMutType = aggregate(MUT$NumberOfSynMutPerSpecies, by = list(MUT$Species,MUT$Subs), FUN = sum); 
names(AggTotalMutSpectrumPerSpeciesPerMutType) = c('Species','MutType','MutTypeRate')
nrow(AggTotalMutSpectrumPerSpeciesPerMutType) # 12357
ALL = merge(AggTotalMutSpectrumPerSpecies,AggTotalMutSpectrumPerSpeciesPerMutType, by = c('Species','MutType'), all.x=TRUE)
ALL[is.na(ALL)]<-0
ALL$Fraction = ALL$MutTypeRate/ALL$TotalMutRate; summary(ALL$Fraction)
ALL = ALL[ALL$Species %in% ListOfSpeciesWithManySubst,]

# check of these species have all 12 nonzero 12 types of substitutions:
ALL$number = 1
AGG = aggregate(ALL$number, by = list(ALL$Species), FUN = sum); summary(AGG$x)
VecOfSpeciesWith12Subst = AGG$Group.1

# create matrix for PCA:

AT = ALL[ALL$MutType == 'A_T',]; AT = AT[c(1,5)]; names(AT) = c('Species','AT'); AT = AT[order(AT$Species),]
AG = ALL[ALL$MutType == 'A_G',]; AG = AG[c(1,5)]; names(AG) = c('Species','AG'); AG = AG[order(AG$Species),]
AC = ALL[ALL$MutType == 'A_C',]; AC = AC[c(1,5)]; names(AC) = c('Species','AC'); AC = AC[order(AC$Species),]
TA = ALL[ALL$MutType == 'T_A',]; TA = TA[c(1,5)]; names(TA) = c('Species','TA'); TA = TA[order(TA$Species),]
TG = ALL[ALL$MutType == 'T_G',]; TG = TG[c(1,5)]; names(TG) = c('Species','TG'); TG = TG[order(TG$Species),]
TC = ALL[ALL$MutType == 'T_C',]; TC = TC[c(1,5)]; names(TC) = c('Species','TC'); TC = TC[order(TC$Species),]
CA = ALL[ALL$MutType == 'C_A',]; CA = CA[c(1,5)]; names(CA) = c('Species','CA'); CA = CA[order(CA$Species),]
CG = ALL[ALL$MutType == 'C_G',]; CG = CG[c(1,5)]; names(CG) = c('Species','CG'); CG = CG[order(CG$Species),]
CT = ALL[ALL$MutType == 'C_T',]; CT = CT[c(1,5)]; names(CT) = c('Species','CT'); CT = CT[order(CT$Species),]
GA = ALL[ALL$MutType == 'G_A',]; GA = GA[c(1,5)]; names(GA) = c('Species','GA'); GA = GA[order(GA$Species),]
GC = ALL[ALL$MutType == 'G_C',]; GC = GC[c(1,5)]; names(GC) = c('Species','GC'); GC = GC[order(GC$Species),]
GT = ALL[ALL$MutType == 'G_T',]; GT = GT[c(1,5)]; names(GT) = c('Species','GT'); GT = GT[order(GT$Species),]
MATRIX = cbind(AT,AG[,2],AC[,2],TA[,2],TG[,2],TC[,2],CA[,2],CG[,2],CT[,2],GA[,2],GC[,2],GT[,2]); names(MATRIX) = c('Species','AT','AG','AC','TA','TG','TC','CA','CG','CT','GA','GC','GT')
MATRIX = merge(MATRIX,AnAge)
table(MATRIX$Class)

pdf('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/4_FIGURES/GAMutSpectrumForDifferenClasses.pdf')
boxplot(MATRIX[MATRIX$Class == 'Actinopterygii',]$GA,MATRIX[MATRIX$Class == 'Amphibia',]$GA,MATRIX[MATRIX$Class == 'Reptilia',]$GA,MATRIX[MATRIX$Class == 'Mammalia',]$GA,MATRIX[MATRIX$Class == 'Aves',]$GA, notch = TRUE, names = c('Actinopterygii','Amphibia','Reptilia','Mammalia','Aves'), ylab = 'Proportion of G>A')
dev.off()

###################################
###### 06.03.2018: tr/tv analysis versus GT; barplots for 12 types; boxplots with asymmetry
###################################

############# GENERATION LENGTH FOR ALL MAMMALS
#  https://natureconservation.pensoft.net/article/1343/  data for generation length of ALL mammals!!!
rm(list=ls(all=TRUE))
GT = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/GenerationLenghtforAllMammals/GenerationLenghtforMammals.xlsx.txt', header = TRUE, sep = '\t')
head(GT)
str(GT)
GT$Species = gsub(' ','_',GT$Scientific_name)
length(unique(GT$Species))
summary(GT$AdultBodyMass_g)     # 2        21        71    136058       614 154321304 
summary(GT$GenerationLength_d)  # 129.0   624.4  1101.1  1578.8  2064.3 18980.0 # max = 18980 days => 52 years. ok
GT = GT[,c(11,13)]
summary(GT$GenerationLength_d)

############ MUTATION SPECTRUM 
# METHOD = 'PARSIMONY'
 METHOD = 'MAXLIKELIHOOD'

if (METHOD == 'MAXLIKELIHOOD') {MUT = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/2_DERIVED/Table_fixed.txt', header = TRUE)}
if (METHOD == 'PARSIMONY')     {MUT = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/2_DERIVED/Table_fixed_parsymony.txt', header = TRUE)}
nrow(MUT)
VecOfNormalSubstitutions = c('A_C','C_A','A_G','G_A','C_G','G_C','C_T','T_C','G_T','T_G','T_A','A_T')

##### FILTER: normal Substitutions
MUT = MUT[MUT$Subs %in% VecOfNormalSubstitutions,]
nrow(MUT)

##### FILTER: Synonymous Substitutions
MUT = MUT[MUT$AncestralAA == MUT$DescendantAA,]
nrow(MUT)

##### FILTER: Gene
table(MUT$Gene)
# MUT = MUT[MUT$Gene == 'CytB',]
MUT = MUT[MUT$Gene == 'COX1' | MUT$Gene == 'COX2' | MUT$Gene == 'COX3',]

##### NORMALIZATION by nucleotide count in the third position of four-fold synonymous substitutions:
NUC = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/2_DERIVED/ATGC_counts_in_SYN_codons_wit_full_gene.txt', header = TRUE)
NUC$Gene = gsub("(.*)\\.",'',NUC$Species)
NUC$Species = gsub("\\.(.*)",'',NUC$Species)
MUT = merge(MUT,NUC, by = c("Species","Gene"))  # compare CountA.x and CountA.y  - they should be identical.
nrow(MUT) # A bit less !!!! WHY????

EXTRACT = function(x) {first = unlist(strsplit(as.character(x),'_'))[1]; return(first);}; MUT$AncestralNuc = apply(as.matrix(MUT$Subs), 1, EXTRACT)
MUT$NumberOfSynMutPerSpecies = 1
MUT_A = MUT[MUT$AncestralNuc == 'A',]; MUT_T = MUT[MUT$AncestralNuc == 'T',]; MUT_G = MUT[MUT$AncestralNuc == 'G',]; MUT_C = MUT[MUT$AncestralNuc == 'C',] # 64145+123587+97657+128195=413584 
MUT_A$NumberOfSynMutPerSpecies = MUT_A$NumberOfSynMutPerSpecies/MUT_A$CountA_Syn;
MUT_T$NumberOfSynMutPerSpecies = MUT_T$NumberOfSynMutPerSpecies/MUT_T$CountT_Syn;
MUT_G$NumberOfSynMutPerSpecies = MUT_G$NumberOfSynMutPerSpecies/MUT_G$CountG_Syn;
MUT_C$NumberOfSynMutPerSpecies = MUT_C$NumberOfSynMutPerSpecies/MUT_C$CountC_Syn;
MUT = rbind(MUT_A,MUT_T,MUT_G,MUT_C)

##### FILTER: fourfold degenerate sites:
VecOfSynFourFoldDegenerateSites = c('CTT','CTC','CTA','CTG','GTT','GTC','GTA','GTG','TCT','TCC','TCA','TCG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','GCT','GCC','GCA','GCG','CGT','CGC','CGA','CGG','GGT','GGC','GGA','GGG'); length(unique(VecOfSynFourFoldDegenerateSites)) # 32
MUT = MUT[MUT$AncestorCodon %in% VecOfSynFourFoldDegenerateSites & MUT$DescendantCodon %in% VecOfSynFourFoldDegenerateSites,]
nrow(MUT)

### VECTOR OF SPECIES WITH MANY SUBSTITUTIONS:
MUT$NumberOfAnyMutPerSpecies = 1
AGG1 = aggregate(MUT$NumberOfAnyMutPerSpecies, by = list(MUT$Species), FUN = sum)
summary(AGG1$x) # 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.0    30.0   100.0   239.5   285.0  2629.0 
ListOfSpeciesWithManySubst = AGG1[AGG1$x >= 30,]$Group.1; length(ListOfSpeciesWithManySubst) # 232
ListOfSpeciesWithSuperManySubst = AGG1[AGG1$x >= 30,]$Group.1; length(ListOfSpeciesWithManySubst) # 232

##### COUNT THE TOTAL NUMBER OF NORMALIZED MUTATIONS PER SPECIES

AggTotalMutSpectrumPerSpecies = aggregate(MUT$NumberOfSynMutPerSpecies, by = list(MUT$Species), FUN = sum); 
names(AggTotalMutSpectrumPerSpecies) = c('Species','TotalMutRate')
MutTypes = data.frame(VecOfNormalSubstitutions); names(MutTypes)=c('MutType')
nrow(AggTotalMutSpectrumPerSpecies) # 1710
AggTotalMutSpectrumPerSpecies = merge(AggTotalMutSpectrumPerSpecies,MutTypes)
nrow(AggTotalMutSpectrumPerSpecies) # 20520/12==1710!
AggTotalMutSpectrumPerSpeciesPerMutType = aggregate(MUT$NumberOfSynMutPerSpecies, by = list(MUT$Species,MUT$Subs), FUN = sum); 
names(AggTotalMutSpectrumPerSpeciesPerMutType) = c('Species','MutType','MutTypeRate')
nrow(AggTotalMutSpectrumPerSpeciesPerMutType) # 12357
ALL = merge(AggTotalMutSpectrumPerSpecies,AggTotalMutSpectrumPerSpeciesPerMutType, by = c('Species','MutType'), all.x=TRUE)
ALL[is.na(ALL)]<-0
ALL$Fraction = ALL$MutTypeRate/ALL$TotalMutRate; summary(ALL$Fraction)
ALL = ALL[ALL$Species %in% ListOfSpeciesWithSuperManySubst,]
ALL = merge(ALL,GT); length(unique(ALL$Species)) # 140

#### several analyses:

GoodOrder = c('T_C','G_A','A_G','C_T','A_C','C_G','G_C','T_A','T_G','A_T','C_A','G_T'); length(GoodOrder)

for (i in 1:length(GoodOrder))
{ # i = 1
  TEMP = ALL[ALL$MutType == GoodOrder[i],]
  Pval = as.numeric(cor.test(TEMP$GenerationLength_d,TEMP$Fraction, method = 'spearman')[3])*12
  Rho = as.numeric(cor.test(TEMP$GenerationLength_d,TEMP$Fraction, method = 'spearman')[4])
  OneLine = data.frame(GoodOrder[i],Pval,Rho)
  if (i == 1) {FINAL = OneLine}
  if (i >  1) {FINAL = rbind(FINAL,OneLine)}
}
FINAL
# GoodOrder.i.         Pval         Rho 
# 1           G_A 8.306945e-01 -0.01129201
# 2           T_C 1.232595e-06  0.25202204 PAPER
# 3           A_G 6.070599e-03  0.14416292
# 4           C_T 3.209208e-02 -0.11283411
# 5           C_A 8.046965e-07 -0.25623635 PAPER
# 6           A_C 2.297772e-03 -0.15998190
# 7           C_G 3.837695e-01 -0.04597583
# 8           G_C 9.972836e-03 -0.13546359
# 9           G_T 6.658911e-07 -0.25808409 PAPER
# 10          T_G 1.841122e-04 -0.19560912 # 1.841122e-04 * 12 = 0.002209346 -> we don't include it and thus our threshold can be bonferroni corrected p value < 0.001
# 11          T_A 5.702422e-03 -0.14522730
# 12          A_T 2.767910e-06 -0.24381436 PAPER

### multiple regression: 

VecOfSubst = unique(ALL$MutType) # C_T G_T T_C C_A G_A C_G A_C A_T T_A A_G G_C T_G
TR = data.frame(ALL[ALL$MutType == 'C_T',]$Species, ALL[ALL$MutType == 'C_T',]$GenerationLength_d, 
               ALL[ALL$MutType == 'C_T',]$Fraction, ALL[ALL$MutType == 'G_T',]$Fraction,  ALL[ALL$MutType == 'T_C',]$Fraction, 
               ALL[ALL$MutType == 'C_A',]$Fraction, ALL[ALL$MutType == 'G_A',]$Fraction,  ALL[ALL$MutType == 'C_G',]$Fraction, 
               ALL[ALL$MutType == 'A_C',]$Fraction, ALL[ALL$MutType == 'A_T',]$Fraction,  ALL[ALL$MutType == 'T_A',]$Fraction, 
               ALL[ALL$MutType == 'A_G',]$Fraction, ALL[ALL$MutType == 'G_C',]$Fraction,  ALL[ALL$MutType == 'T_G',]$Fraction)
names(TR)=c('Species','GenerationLength_d','C_T','G_T','T_C','C_A','G_A','C_G','A_C','A_T','T_A','A_G','G_C','T_G')

mean(TR$G_A) # 0.612 PAPER
mean(TR$T_C) # 0.164 PAPER
mean(TR$C_T) # 0.085 PAPER
mean(TR$G_C) # 0.028 PAPER
mean(TR$T_A) # 0.027 PAPER
mean(TR$A_G) # 0.024 PAPER
mean(TR$G_T) # 0.022 PAPER
mean(TR$C_A) # 0.020 PAPER
mean(TR$A_T) # 0.006 PAPER
mean(TR$A_C) # 0.006 PAPER
mean(TR$T_G) # 0.003 PAPER
mean(TR$C_G) # 0.002 PAPER

# 0.612+0.164 - 0.774; 3/4 = 0.75 PAPER

a<-lm(log2(GenerationLength_d) ~ T_C + A_T + C_A + G_T, data = TR)
summary(a)
a<-lm(log2(GenerationLength_d) ~ G_A + T_C + A_T + C_A + G_T, data = TR)
summary(a)
a<-lm(log2(GenerationLength_d) ~ scale(T_C) + scale(A_T) + scale(C_A) + scale(G_T), data = TR)  #### PAPER
summary(a)
a<-lm(log2(GenerationLength_d) ~ scale(G_A) + scale(T_C) + scale(A_T) + scale(C_A) + scale(G_T), data = TR)  #### PAPER
summary(a)

TR$GC_TA = TR$G_T + TR$C_A
a<-lm(log2(GenerationLength_d) ~ scale(T_C) + scale(A_T) + scale(TR$GC_TA), data = TR)
summary(a)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     10.34591    0.05927 174.561  < 2e-16 ***
#  scale(T_C)       0.24628    0.06065   4.061 6.02e-05 ***
#  scale(A_T)      -0.18224    0.05963  -3.056 0.002411 ** 
#  scale(TR$GC_TA) -0.23376    0.06063  -3.856 0.000137 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 1.126 on 357 degrees of freedom
# Multiple R-squared:  0.1277,	Adjusted R-squared:  0.1204 
# F-statistic: 17.43 on 3 and 357 DF,  p-value: 1.399e-10

TR$AllTs = TR$A_G + TR$G_A + TR$C_T + TR$T_C
TR$AllTv = TR$C_A + TR$A_C + TR$C_G + TR$G_C + TR$G_T + TR$T_G + TR$T_A + TR$A_T 
TR$TsTv = TR$AllTs/TR$AllTv

if (METHOD == 'PARSIMONY') {pdf("/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/4_FIGURES/TsTvVsGtFourScatterplots.PARSIMONY.pdf")}
if (METHOD == 'MAXLIKELIHOOD') {pdf("/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/4_FIGURES/TsTvVsGtFourScatterplots.MAXLIKELIHOOD.pdf")}  
par(mfrow=c(2,2))
plot(log2(TR$T_C),log2(TR$GenerationLength_d), main = 'T>C')
plot(log2(TR$GC_TA),log2(TR$GenerationLength_d), main = 'GC>TA')
plot(log2(TR$A_T),log2(TR$GenerationLength_d), main = 'A>T')
plot(log2(TR$TsTv),log2(TR$GenerationLength_d), main = 'Ts/Tv')
dev.off()
 
### single strand 

TR$TcAg = TR$T_C / TR$A_G; TcAg = TR[!is.na(TR$TcAg) & TR$TcAg > 0 & TR$TcAg < Inf,]$TcAg;
TR$GaCt = TR$G_A / TR$C_T; GaCt = TR[!is.na(TR$GaCt) & TR$GaCt > 0 & TR$GaCt < Inf,]$GaCt;
TR$GcCg = TR$G_C / TR$C_G; GcCg = TR[!is.na(TR$GcCg) & TR$GcCg > 0 & TR$GcCg < Inf,]$GcCg;
TR$TaAt = TR$T_A / TR$A_T; TaAt = TR[!is.na(TR$TaAt) & TR$TaAt > 0 & TR$TaAt < Inf,]$TaAt;
TR$GtCa = TR$G_T / TR$C_A; GtCa = TR[!is.na(TR$GtCa) & TR$GtCa > 0 & TR$GtCa < Inf,]$GtCa;
TR$AcTg = TR$A_C / TR$T_G; AcTg = TR[!is.na(TR$AcTg) & TR$AcTg > 0 & TR$AcTg < Inf,]$AcTg;

if (METHOD == 'PARSIMONY') {pdf('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/4_FIGURES/StrandBias6Ratios.PARSIMONY.pdf')} # PAPER
if (METHOD == 'MAXLIKELIHOOD') {pdf('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/4_FIGURES/StrandBias6Ratios.MAXLIKELIHOOD.pdf')} # PAPER
boxplot(TcAg,GaCt,GcCg,TaAt,GtCa,AcTg, names = c('TC/AG','GA/CT','GC/CG','TA/AT','GT/CA','AC/TG'), outline = FALSE, notch = TRUE)
abline(h = 1, col = 'red')
dev.off()

wilcox.test(AcTg, mu = 1) # p = 0.14           ; PARS: 0.8506
wilcox.test(GtCa, mu = 1) # p-value = 9.709e-14; PARS: 1.173e-05
wilcox.test(TaAt, mu = 1) # p-value < 2.2e-16  ; PARS: 9.851e-15
wilcox.test(GcCg, mu = 1) # p-value < 2.2e-16  ; PARS: 9.851e-15
wilcox.test(GaCt, mu = 1) # p-value < 2.2e-16
wilcox.test(TcAg, mu = 1) # p-value < 2.2e-16

head(TR)
##################### strand bias and generation time (nothing significant): PAPER
cor.test(TR[!is.na(TR$TcAg) & TR$TcAg > 0 & TR$TcAg < Inf,]$GenerationLength_d,TR[!is.na(TR$TcAg) & TR$TcAg > 0 & TR$TcAg < Inf,]$TcAg, method = 'spearman'); nrow(TR[!is.na(TR$TcAg) & TR$TcAg > 0 & TR$TcAg < Inf,]) # 351
cor.test(TR[!is.na(TR$GaCt) & TR$GaCt > 0 & TR$GaCt < Inf,]$GenerationLength_d,TR[!is.na(TR$GaCt) & TR$GaCt > 0 & TR$GaCt < Inf,]$GaCt, method = 'spearman'); nrow(TR[!is.na(TR$GaCt) & TR$GaCt > 0 & TR$GaCt < Inf,]) # 359
cor.test(TR[!is.na(TR$GcCg) & TR$GcCg > 0 & TR$GcCg < Inf,]$GenerationLength_d,TR[!is.na(TR$GcCg) & TR$GcCg > 0 & TR$GcCg < Inf,]$GcCg, method = 'spearman'); nrow(TR[!is.na(TR$GcCg) & TR$GcCg > 0 & TR$GcCg < Inf,]) # 108
cor.test(TR[!is.na(TR$TaAt) & TR$TaAt > 0 & TR$TaAt < Inf,]$GenerationLength_d,TR[!is.na(TR$TaAt) & TR$TaAt > 0 & TR$TaAt < Inf,]$TaAt, method = 'spearman'); nrow(TR[!is.na(TR$TaAt) & TR$TaAt > 0 & TR$TaAt < Inf,]) # 267
cor.test(TR[!is.na(TR$GtCa) & TR$GtCa > 0 & TR$GtCa < Inf,]$GenerationLength_d,TR[!is.na(TR$GtCa) & TR$GtCa > 0 & TR$GtCa < Inf,]$GtCa, method = 'spearman'); nrow(TR[!is.na(TR$GtCa) & TR$GtCa > 0 & TR$GtCa < Inf,]) # 128
cor.test(TR[!is.na(TR$AcTg) & TR$AcTg > 0 & TR$AcTg < Inf,]$GenerationLength_d,TR[!is.na(TR$AcTg) & TR$AcTg > 0 & TR$AcTg < Inf,]$AcTg, method = 'spearman'); nrow(TR[!is.na(TR$AcTg) & TR$AcTg > 0 & TR$AcTg < Inf,]) # 97

dataset = data.frame(ALL[ALL$MutType == 'C_T',]$Species,ALL[ALL$MutType == 'C_T',]$GenerationLength_d,ALL[ALL$MutType == 'T_C',]$Fraction,ALL[ALL$MutType == 'C_T',]$Fraction,ALL[ALL$MutType == 'G_A',]$Fraction,ALL[ALL$MutType == 'A_G',]$Fraction,ALL[ALL$MutType == 'G_T',]$Fraction,ALL[ALL$MutType == 'C_A',]$Fraction,ALL[ALL$MutType == 'A_T',]$Fraction,ALL[ALL$MutType == 'T_A',]$Fraction,ALL[ALL$MutType == 'T_G',]$Fraction,ALL[ALL$MutType == 'A_C',]$Fraction); 
names(dataset)=c('Species','GenerationLength_d','T_C','C_T','G_A','A_G','G_T','C_A','A_T','T_A','T_G','A_C')

dataset$TcAssym = dataset$T_C/(dataset$T_C + dataset$A_G); summary(dataset[dataset$T_C > 0 & dataset$A_G > 0,]$TcAssym) # median = 0.87
dataset$GaAssym = dataset$G_A/(dataset$G_A + dataset$T_C); summary(dataset[dataset$G_A > 0 & dataset$T_C > 0,]$GaAssym) # median = 0.79

dataset$GtAssym = dataset$G_T/(dataset$G_T + dataset$C_A); summary(dataset[dataset$G_T > 0 & dataset$C_A > 0,]$GtAssym) # median = 0.79
dataset$CaAssym = dataset$C_A/(dataset$C_A + dataset$G_T); summary(dataset[dataset$G_T > 0 & dataset$C_A > 0,]$CaAssym) # median = 0.79
dataset$AtAssym = dataset$A_T/(dataset$A_T + dataset$T_A); summary(dataset[dataset$A_T > 0 & dataset$T_A > 0,]$AtAssym) # median = 0.79
dataset$TgAssym = dataset$T_G/(dataset$T_G + dataset$A_C); summary(dataset[dataset$T_G > 0 & dataset$A_C > 0,]$TgAssym) #

summary(dataset$GenerationLength_d) # 624 1384 2496

if (METHOD == 'PARSIMONY') {pdf('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/4_FIGURES/BoxplotsAsymmetry.PARSIMONY.pdf')}
if (METHOD == 'MAXLIKELIHOOD') {pdf('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/4_FIGURES/BoxplotsAsymmetry.MAXLIKELIHOOD.pdf')}
boxplot(dataset[dataset$T_C > 0 & dataset$A_G > 0,]$TcAssym,dataset[dataset$G_A > 0 & dataset$T_C > 0,]$GaAssym,dataset[dataset$G_T > 0 & dataset$C_A > 0,]$GtAssym,dataset[dataset$G_T > 0 & dataset$C_A > 0,]$CaAssym,dataset[dataset$A_T > 0 & dataset$T_A > 0,]$AtAssym,dataset[dataset$T_G > 0 & dataset$A_C > 0,]$TgAssym, notch = TRUE, outline = FALSE, names = c('T>C','G>A','G>T','C>A','A>T','T>G'))
abline(h = 0.5, col = 'red')
dev.off()

dataset1 = dataset[dataset$GenerationLength_d <= 1384,]
dataset2 = dataset[dataset$GenerationLength_d > 1384,]
dataset = dataset1
boxplot(dataset[dataset$T_C > 0 & dataset$A_G > 0,]$TcAssym,dataset[dataset1$G_A > 0 & dataset$T_C > 0,]$GaAssym,dataset[dataset$G_T > 0 & dataset$C_A > 0,]$GtAssym,dataset[dataset$G_T > 0 & dataset$C_A > 0,]$CaAssym,dataset[dataset$A_T > 0 & dataset$T_A > 0,]$AtAssym,dataset[dataset$T_G > 0 & dataset$A_C > 0,]$TgAssym, notch = TRUE, outline = FALSE, names = c('T>C','G>A','G>T','C>A','A>T','T>G'))
abline(h = 0.5, col = 'red')
dataset = dataset2
boxplot(dataset[dataset$T_C > 0 & dataset$A_G > 0,]$TcAssym,dataset[dataset1$G_A > 0 & dataset$T_C > 0,]$GaAssym,dataset[dataset$G_T > 0 & dataset$C_A > 0,]$GtAssym,dataset[dataset$G_T > 0 & dataset$C_A > 0,]$CaAssym,dataset[dataset$A_T > 0 & dataset$T_A > 0,]$AtAssym,dataset[dataset$T_G > 0 & dataset$A_C > 0,]$TgAssym, notch = TRUE, outline = FALSE, names = c('T>C','G>A','G>T','C>A','A>T','T>G'))
abline(h = 0.5, col = 'red')

cor.test(dataset$TcAssym,dataset$GenerationLength_d,method = 'spearman') # nothing
cor.test(dataset$GaAssym,dataset$GenerationLength_d,method = 'spearman') # negative!!!! in short=lived is more asymmetrical. saturation...
cor.test(dataset$GtAssym,dataset$GenerationLength_d,method = 'spearman') # negative!
cor.test(dataset$CaAssym,dataset$GenerationLength_d,method = 'spearman') # positive!
cor.test(dataset$AtAssym,dataset$GenerationLength_d,method = 'spearman') # negative!
cor.test(dataset$TgAssym,dataset$GenerationLength_d,method = 'spearman') # negative!

#### barplot:

ALL = ALL[order(ALL$GenerationLength_d),]; VecOfSpecies = unique(ALL$Species)

for (i in 1:length(VecOfSpecies))
{ # i = 1
  TEMP = ALL[ALL$Species == VecOfSpecies[i],]
  TEMP = TEMP[order(match(TEMP$MutType,GoodOrder)),]
  if (i == 1) {FINAL = data.frame(TEMP$Fraction)}
  if (i >  1) {FINAL = cbind(FINAL,TEMP$Fraction)}
}
names(FINAL)=VecOfSpecies
row.names(FINAL)=c(GoodOrder)

if (METHOD == 'PARSIMONY') {pdf('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/4_FIGURES/BarPlot12MutTypes.PARSIMONY.pdf', width = 60)}
if (METHOD == 'MAXLIKELIHOOD') {pdf('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/4_FIGURES/BarPlot12MutTypes.MAXLIKELIHOOD.pdf', width = 60)}
par(las=2)
par(mai=c(2.4,0.82,0.82,0.42))
# default par(mai=c(1.02,0.82,0.82,0.42))
FinalForBarPLot = as.matrix(FINAL)
VecOfColors = c('green4','green','green3','greenyellow','steelblue1','steelblue3','royalblue1','royalblue4','orchid1','orchid4','purple','purple4') # 'mediumblue','navy')
barplot(FinalForBarPLot,col=VecOfColors, legend.text=TRUE, ylab = 'Proportion of substitutions')
dev.off()

ncol(FINAL) # 361 = take first 50 and last 50 
a = ncol(FINAL)-30;
b = ncol(FINAL);
FinalShort = FINAL[,c(1:30,a:b)]

pdf('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/4_FIGURES/BarPlot12MutTypesShortLong.pdf', width = 30)
par(las=2)
par(mai=c(2.4,0.82,0.82,0.42))
# default par(mai=c(1.02,0.82,0.82,0.42))
FinalForBarPLot = as.matrix(FinalShort)
VecOfColors = c('green4','green','green3','greenyellow','steelblue1','steelblue3','royalblue1','royalblue4','orchid1','orchid4','purple','purple4') # 'mediumblue','navy')
barplot(FinalForBarPLot,col=VecOfColors, legend.text=TRUE)
# legend(c('C_T'), col = (VecOfColors))
dev.off()


###### Tr and Tv
MUT = MUT[MUT$Species %in% ListOfSpeciesWithManySubst,]
nrow(MUT)

VecOfTransitionSubstitutions = c('A_G','G_A','C_T','T_C') # all tr
VecOfTransversionSubstitutions = c('C_A','A_C','C_G','G_C','G_T','T_G','T_A','A_T') # ALL Tv 

MUT_TR = MUT[MUT$Subs %in% VecOfTransitionSubstitutions,]
MUT_TV = MUT[MUT$Subs %in% VecOfTransversionSubstitutions,]

AGG1 = aggregate(MUT_TR$NumberOfSynMutPerSpecies, by = list(MUT_TR$Species), FUN = sum); names(AGG1) = c('Species','Tr')
AGG2 = aggregate(MUT_TV$NumberOfSynMutPerSpecies, by = list(MUT_TV$Species), FUN = sum); names(AGG2) = c('Species','Tv')
TrTv = merge(AGG1,AGG2, by = 'Species')
TrTv = merge(TrTv,GT,by = 'Species'); summary(TrTv$GenerationLength_d)
TrTv$TrTv = TrTv$Tr/TrTv$Tv; summary(TrTv$TrTv); TrTv = TrTv[!is.na(TrTv$TrTv),]

hist(TrTv$TrTv, breaks = 100)
cor.test(TrTv$GenerationLength_d,TrTv$TrTv, method = 'spearman') # rho = 0.2541812 , p = p-value = 1.063e-06
a<-lm(log2(TrTv$GenerationLength_d) ~ scale(TrTv$TrTv)); summary(a) # all genes: 1.668e-01; CYTB: 2.666e-01; COX1: 2.785e-01;
pdf('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/4_FIGURES/ScatterplotTrTv.pdf')
plot(log2(TrTv$GenerationLength_d),log2(TrTv$TrTv), xlab = c('Generation Time (log2(days))'), ylab = 'log2(Tr/Tv)', main = 'N = 359 species')
dev.off()

#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      10.33388    0.06126 168.685  < 2e-16 ***
#  scale(TrTv$TrTv)  0.28266    0.06135   4.607 5.68e-06 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 1.161 on 357 degrees of freedom
# Multiple R-squared:  0.05613,	Adjusted R-squared:  0.05348 
# F-statistic: 21.23 on 1 and 357 DF,  p-value: 5.68e-06

###################################
###### 06.03.2018: Gene specificity 
###################################

############# GENERATION LENGTH FOR ALL MAMMALS
#  https://natureconservation.pensoft.net/article/1343/  data for generation length of ALL mammals!!!
rm(list=ls(all=TRUE))
GT = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/GenerationLenghtforAllMammals/GenerationLenghtforMammals.xlsx.txt', header = TRUE, sep = '\t')
head(GT)
str(GT)
GT$Species = gsub(' ','_',GT$Scientific_name)
length(unique(GT$Species))
summary(GT$AdultBodyMass_g)     # 2        21        71    136058       614 154321304 
summary(GT$GenerationLength_d)  # 129.0   624.4  1101.1  1578.8  2064.3 18980.0 # max = 18980 days => 52 years. ok
GT = GT[,c(11,13)]
summary(GT$GenerationLength_d)

############ MUTATION SPECTRUM 
MUT = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/2_DERIVED/Table_fixed.txt', header = TRUE)
nrow(MUT)
VecOfNormalSubstitutions = c('A_C','C_A','A_G','G_A','C_G','G_C','C_T','T_C','G_T','T_G','T_A','A_T')

##### FILTER: normal Substitutions
MUT = MUT[MUT$Subs %in% VecOfNormalSubstitutions,]
nrow(MUT)

##### FILTER: Synonymous Substitutions
MUT = MUT[MUT$AncestralAA == MUT$DescendantAA,]
nrow(MUT)

##### FILTER: Gene
table(MUT$Gene)
VecOfGenes = unique(MUT$Gene)
MUT_ALL = MUT
for (gene in 1:length(VecOfGenes))
  { # gene = 9
  MUT = MUT_ALL[MUT_ALL$Gene == VecOfGenes[gene],]
  
  ##### NORMALIZATION: the third position of four-fold synonymous substitutions:
  NUC = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/2_DERIVED/ATGC_counts_in_SYN_codons_wit_full_gene.txt', header = TRUE)
  NUC$Gene = gsub("(.*)\\.",'',NUC$Species)
  NUC$Species = gsub("\\.(.*)",'',NUC$Species)
  MUT = merge(MUT,NUC, by = c("Species","Gene"))  # compare CountA.x and CountA.y  - they should be identical.
  nrow(MUT) # A bit less !!!! WHY????

  EXTRACT = function(x) {first = unlist(strsplit(as.character(x),'_'))[1]; return(first);}; MUT$AncestralNuc = apply(as.matrix(MUT$Subs), 1, EXTRACT)
  MUT$NumberOfSynMutPerSpecies = 1
  MUT_A = MUT[MUT$AncestralNuc == 'A',]; MUT_T = MUT[MUT$AncestralNuc == 'T',]; MUT_G = MUT[MUT$AncestralNuc == 'G',]; MUT_C = MUT[MUT$AncestralNuc == 'C',] # 64145+123587+97657+128195=413584 
  MUT_A$NumberOfSynMutPerSpecies = MUT_A$NumberOfSynMutPerSpecies/MUT_A$CountA_Syn;
  MUT_T$NumberOfSynMutPerSpecies = MUT_T$NumberOfSynMutPerSpecies/MUT_T$CountT_Syn;
  MUT_G$NumberOfSynMutPerSpecies = MUT_G$NumberOfSynMutPerSpecies/MUT_G$CountG_Syn;
  MUT_C$NumberOfSynMutPerSpecies = MUT_C$NumberOfSynMutPerSpecies/MUT_C$CountC_Syn;
  MUT = rbind(MUT_A,MUT_T,MUT_G,MUT_C)

  ##### FILTER: fourfold degenerate sites:
  VecOfSynFourFoldDegenerateSites = c('CTT','CTC','CTA','CTG','GTT','GTC','GTA','GTG','TCT','TCC','TCA','TCG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','GCT','GCC','GCA','GCG','CGT','CGC','CGA','CGG','GGT','GGC','GGA','GGG'); length(unique(VecOfSynFourFoldDegenerateSites)) # 32
  MUT = MUT[MUT$AncestorCodon %in% VecOfSynFourFoldDegenerateSites & MUT$DescendantCodon %in% VecOfSynFourFoldDegenerateSites,]
  nrow(MUT)

  ######## CALCULATE MATRIX
  ### CALCULATE MATRIX 1: TAKE ONLY SPECIES WITH MANY SUBSTITUTIONS:
  MUT$NumberOfAnyMutPerSpecies = 1
  AGG1 = aggregate(MUT$NumberOfAnyMutPerSpecies, by = list(MUT$Species), FUN = sum)
  summary(AGG1$x) # 
  ListOfSpeciesWithManySubst = AGG1[AGG1$x >= 20,]$Group.1; length(ListOfSpeciesWithManySubst) # 232
  MUT = MUT[MUT$Species %in% ListOfSpeciesWithManySubst,]
  nrow(MUT)

  ###### Tr and Tv
  # VecOfTransitionSubstitutions = c('A_G','G_A','C_T','T_C') # all tr
  # VecOfTransitionSubstitutions = c('G_A','T_C') # transition in SS DNA
  VecOfTransitionSubstitutions = c('T_C') # transition in SS DNA
  VecOfTransversionSubstitutions = c('C_A','A_C','C_G','G_C','G_T','T_G','T_A','A_T') # ALL Tv 

  MUT$NumberOfSynMutPerSpecies = 1

  MUT_TR = MUT[MUT$Subs %in% VecOfTransitionSubstitutions,]
  MUT_TV = MUT[MUT$Subs %in% VecOfTransversionSubstitutions,]

  AGG1 = aggregate(MUT_TR$NumberOfSynMutPerSpecies, by = list(MUT_TR$Species), FUN = sum); names(AGG1) = c('Species','Tr')
  AGG2 = aggregate(MUT_TV$NumberOfSynMutPerSpecies, by = list(MUT_TV$Species), FUN = sum); names(AGG2) = c('Species','Tv')
  TrTv = merge(AGG1,AGG2, by = 'Species')
  TrTv$TrTv = TrTv$Tr/TrTv$Tv; summary(TrTv$TrTv); TrTv = TrTv[!is.na(TrTv$TrTv),]
  OneLine = TrTv
  OneLine$Gene = VecOfGenes[gene]
  if (gene == 1) {Final = OneLine}
  if (gene >  1) {Final = rbind(Final,OneLine)}
}

COX1 = Final[Final$Gene == 'COX1',]; COX1 = COX1[,c(1,4)]; names(COX1)=c('Species','TrTvCox1')
COX2 = Final[Final$Gene == 'COX2',]; COX2 = COX2[,c(1,4)]; names(COX2)=c('Species','TrTvCox2')
COX3 = Final[Final$Gene == 'COX3',]; COX3 = COX3[,c(1,4)]; names(COX3)=c('Species','TrTvCox3')
ATP8 = Final[Final$Gene == 'ATP8',]; ATP8 = ATP8[,c(1,4)]; names(ATP8)=c('Species','TrTvATP8')
ATP6 = Final[Final$Gene == 'ATP6',]; ATP6 = ATP6[,c(1,4)]; names(ATP6)=c('Species','TrTvATP6')
ND1  = Final[Final$Gene == 'ND1',];  ND1  =  ND1[,c(1,4)]; names(ND1) =c('Species','TrTvND1')
ND2  = Final[Final$Gene == 'ND2',];  ND2  =  ND2[,c(1,4)]; names(ND2) =c('Species','TrTvND2')
CytB = Final[Final$Gene == 'CytB',]; CytB = CytB[,c(1,4)]; names(CytB)=c('Species','TrTvCytB')

ShortSSD = Final[Final$Gene %in% c('COX1','COX2','COX3','ATP6','ATP8','ND4','ND4L'),]
ShortSSDAgg = aggregate(ShortSSD$TrTv, by = list(ShortSSD$Species), FUN = mean); names(ShortSSDAgg) = c('Species','TrTv')

COMP = merge(ShortSSDAgg,CytB)
# COMP = merge(COX1,CytB)
COMP$TrTvRatio = COMP$TrTvCytB/COMP$TrTv
summary(COMP$TrTvRatio)
COMP = merge(COMP,GT)
cor.test(COMP$TrTvRatio,COMP$GenerationLength_d, method = 'spearman') # negative... the bigger generation time the less difference between COX and CYTB. 
summary(COMP[COMP$GenerationLength_d <= quantile(COMP$GenerationLength_d,0.25),]$TrTvRatio) # 0.2643  0.7922  1.1111  1.3743  1.4252  4.7222 
summary(COMP[COMP$GenerationLength_d <= quantile(COMP$GenerationLength_d,0.5),]$TrTvRatio)  # 0.2643  0.7197  1.0256  1.2266  1.2500  4.7222
summary(COMP[COMP$GenerationLength_d >= quantile(COMP$GenerationLength_d,0.75),]$TrTvRatio) #### !!!!! ???? In long lived tr/tv(CYTB) is not higher than Tr/Tv of other genes????
summary(COMP[COMP$GenerationLength_d >= quantile(COMP$GenerationLength_d,0.5),]$TrTvRatio) #### !!!!! ???? In long lived tr/tv(CYTB) is not higher than Tr/Tv of other genes????
boxplot(COMP[COMP$GenerationLength_d <= quantile(COMP$GenerationLength_d,0.25),]$TrTvRatio,COMP[COMP$GenerationLength_d >= quantile(COMP$GenerationLength_d,0.75),]$TrTvRatio, notch = TRUE)
boxplot(COMP[COMP$GenerationLength_d <= quantile(COMP$GenerationLength_d,0.5),]$TrTvRatio,COMP[COMP$GenerationLength_d >= quantile(COMP$GenerationLength_d,0.5),]$TrTvRatio, notch = TRUE)

### saturation curve (especially for G>A and T>C) is more flat (or even goes down) for long-lived mammals than for short-lived (do I see it also for AT/GC skew in whole genome data...?):
### saturation in big animals in CYTB!!!!!! YEYEYEYE
### OR THERE iS a trade off between rate of mutations and bias...
### OR alternative OH, OL






















##### 
##### the strongest correlations:
# VecOfTransitionSubstitutions = c('C_T','T_C') vs VecOfTransversionSubstitutions = c('A_C','C_A','C_G','G_C','G_T','T_G','T_A','A_T'): rho = 0.37
# VecOfTransitionSubstitutions = c('T_C') vs VecOfTransversionSubstitutions = c('A_C','C_A','C_G','G_C','G_T','T_G','T_A','A_T'): rho = 0.38
# VecOfTransitionSubstitutions = c('T_C') vs VecOfTransversionSubstitutions = c('T_A','A_T'): rho = 0.40
# VecOfTransitionSubstitutions = c('T_C') vs VecOfTransversionSubstitutions = c('A_T'): rho = 0.39  

# VecOfTransitionSubstitutions = c('A_G','G_A','C_T','T_C') # all tr
# VecOfTransitionSubstitutions = c('G_A','T_C') # Ts, which depends on the time of beeng single-stranded 
VecOfTransitionSubstitutions = c('G_A','T_C') # they are driven by DNA polymerase?
# VecOfTransitionSubstitutions = c('T_C') # 

# VecOfTransversionSubstitutions = c('C_A','A_C','C_G','G_C','G_T','T_G','T_A','A_T') # ALL Tv 
# VecOfTransversionSubstitutions = c('C_T','A_G')
# VecOfTransversionSubstitutions = c('A_G') #
VecOfTransversionSubstitutions = c('C_A','G_T')  # Tv, induced by ROS
# VecOfTransversionSubstitutions = c('G_T')  # no correlation between tr/tv and GT
# VecOfTransversionSubstitutions = c('C_A')    # there is a correlation between tr/tv and GT: p-value = 0.0007141, rho = 0.255
# VecOfTransversionSubstitutions = c('A_C','C_G','G_C','T_G','T_A','A_T')  # Tv not induced by ROS
# VecOfTransversionSubstitutions = c('A_T') # ALL Tv
# VecOfTransversionSubstitutions = c('C_G','G_C','T_G') # ALL Tv


ECO = unique(data.frame(MUT$Species, MUT$Female.maturity..days.)); names(ECO)=c('Species','GenerationTime')
ECO = ECO[!is.na(ECO$GenerationTime),] # 197!!! FUCK, WHY LESS THAN BEFORE!!!!
TrTv = merge(TrTv,ECO, by = 'Species')  # FUCK, AGAIN LESS!!!!
TrTv = merge(TrTv,NumberOfSubst, by = 'Species')
# TrTv = merge(TrTv,GC, by = 'Species')

A = cor.test(TrTv$TrTv,TrTv$GenerationTime, method = 'spearman') # positive!!!
# B = cor.test(TrTv$TrTv,TrTv$GcContent, method = 'spearman') # positive but week!!!
C<-lm(log2(TrTv$TrTv) ~ log2(TrTv$GenerationTime) + log2(TrTv$GcContent)); summary(C)  # GT is more significant than GC

PVal = round(as.numeric(A[3]),10)
Rho  = round(as.numeric(A[4]),5)
N = nrow(TrTv)
pdf('TrTvDirectVersusGenTimeFor188Mammals.pdf')
# pdf('TrTvNormalizedVersusGenTimeFor188Mammals.pdf')
title = paste('Spearman Rho = ',Rho,', Pval = ',PVal,', (N=',N,')', sep='')
plot(log2(TrTv$GenerationTime),log2(TrTv$TrTv), ylab = 'log2(Tr/Tv)', xlab = 'log2(generation time in days)', main = title)
dev.off()

write.table(TrTv, file = '204Species.AllSynon.NoNormalisationByNuclCount.txt', quote = FALSE, row.names = FALSE, sep = '\t')




















###################################
###### 21.02.2018: analysis of Other Chordata from Alina's file "/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/cds_at_gc_skew.txt"
###################################

###################### Actinopterygii
rm(list=ls(all=TRUE))
Skew = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/cds_at_gc_skew.txt', header = TRUE)
names(Skew) = c('Species','Gene','AtSkew','GcSkew')

AnAge = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/1_RAW_DATA/anage_data.txt', header = TRUE, sep = '\t')
table(AnAge$Class) # 
#Actinopterygii           Amphibia               Aves Cephalaspidomorphi     Chondrichthyes           Mammalia           Reptilia      Sarcopterygii 
#824                174               1191                 16                116               1327                544                  4 
AnAge = AnAge[AnAge$Class == 'Actinopterygii' & ! is.na(AnAge$Female.maturity..days.),]; 
AnAge$Species = paste(AnAge$Genus,AnAge$Species, sep='_') 
AnAge = data.frame(AnAge$Species,AnAge$Female.maturity..days); names(AnAge)=c('Species','FemaleMaturityDays')

Skew = merge(Skew,AnAge, by = 'Species')
table(Skew$Gene) # ATP6 ATP8 COX1 COX2 COX3 CytB  ND1  ND2  ND3  ND4 ND4L  ND5  ND6

Skew = Skew[Skew$Gene %in% c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','CytB','ND1','ND2'),] # no ND6, ND1 and ND2 are located at the end approximately

Gene = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','CytB', 'ND1','ND2')
Timing = seq(1:12)
NewData = data.frame(Gene,Timing)

Skew = merge(Skew,NewData)
Skew$Timing = as.factor(Skew$Timing)
boxplot(Skew$AtSkew ~ Skew$Timing, outline = FALSE, notch = TRUE, names = Gene)
boxplot(Skew$GcSkew ~ Skew$Timing, outline = FALSE, notch = TRUE, names = Gene)

summary(Skew$FemaleMaturityDays)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 60     784    1277    1949    2153    9855 

SkewShort = Skew[Skew$FemaleMaturityDays <= 784,];  SkewShort$GT = 'Short'
SkewLong = Skew[Skew$FemaleMaturityDays >= 2153,];  SkewLong$GT = 'Long'
SkewExtremes = rbind(SkewShort,SkewLong)

SkewShort = Skew[Skew$FemaleMaturityDays <= 784,];  SkewShort$GT = 'Short'
SkewLong = Skew[Skew$FemaleMaturityDays >= 2153,];  SkewLong$GT = 'Long'
SkewExtremes = rbind(SkewShort,SkewLong)

pdf('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/4_FIGURES/AtGcSkewActinopterigii.pdf', width = 14)
par(cex.axis=0.7)
par(mfrow = c(2,1))

par(mfrow = c(2,1))
boxplot(Skew$AtSkew ~ Skew$Timing, outline = FALSE, notch = TRUE, ylim = c(-0.1,1), main = 'AT skew for genes with different time being single stranded (All Actinopterygii)', names = Gene)
boxplot(Skew$GcSkew ~ Skew$Timing, outline = FALSE, notch = TRUE, ylim = c(-1,0.1), main = 'GC skew for genes with different time being single stranded (All Actinopterygii)', names = Gene)

GeneVec = c('L.Cox1','S.Cox1','L.Cox2','S.Cox2','L.Atp8','S.Atp8','L.Atp6','S.Atp6','L.Cox3','S.Cox3','L.Nd3','S.Nd3','L.Nd4l','S.Nd4l','L.Nd4','S.Nd4','L.Nd5','S.Nd5','L.Cytb','S.Cytb','L.Nd1','S.Nd1','L.Nd2','S.Nd2')
par(mfrow = c(2,1))
boxplot(SkewExtremes$AtSkew ~ SkewExtremes$GT*SkewExtremes$Timing, outline = FALSE, notch = TRUE, ylim = c(-0.1,1), col = c('green','red'), main = paste('AT skew for genes with different time being single stranded','(Red - Short GT < 784 days; Green - Long GT: > 2153 days)',sep = '\n'), names = GeneVec)
abline(h = 0, lwd = 2, lt = 2, col = 'red')
boxplot(SkewExtremes$GcSkew ~ SkewExtremes$GT*SkewExtremes$Timing, outline = FALSE, notch = TRUE, ylim = c(-1,0.1), col = c('green','red'), main = paste('GC skew for genes with different time being single stranded','(Red - Short GT < 784 days; Green - Long GT: > 2153 days)',sep = '\n'), names = GeneVec)
abline(h = 0, lwd = 2, lt = 2, col = 'red')
dev.off()

################### Aves
rm(list=ls(all=TRUE))
Skew = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/cds_at_gc_skew.txt', header = TRUE)
names(Skew) = c('Species','Gene','AtSkew','GcSkew')

AnAge = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/1_RAW_DATA/anage_data.txt', header = TRUE, sep = '\t')
table(AnAge$Class) # 
#Actinopterygii           Amphibia               Aves Cephalaspidomorphi     Chondrichthyes           Mammalia           Reptilia      Sarcopterygii 
#824                174               1191                 16                116               1327                544                  4 
AnAge = AnAge[AnAge$Class == 'Aves' & ! is.na(AnAge$Female.maturity..days.),]; 
AnAge$Species = paste(AnAge$Genus,AnAge$Species, sep='_') 
AnAge = data.frame(AnAge$Species,AnAge$Female.maturity..days); names(AnAge)=c('Species','FemaleMaturityDays')

Skew = merge(Skew,AnAge, by = 'Species')
table(Skew$Gene) # ATP6 ATP8 COX1 COX2 COX3 CytB  ND1  ND2  ND3  ND4 ND4L  ND5  ND6

Skew = Skew[Skew$Gene %in% c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','CytB','ND1','ND2'),] # no ND6, ND1 and ND2 are located at the end approximately

Gene = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','CytB', 'ND1','ND2')
Timing = seq(1:12)
NewData = data.frame(Gene,Timing)

Skew = merge(Skew,NewData)
Skew$Timing = as.factor(Skew$Timing)
boxplot(Skew$AtSkew ~ Skew$Timing, outline = FALSE, notch = TRUE, names = Gene)
boxplot(Skew$GcSkew ~ Skew$Timing, outline = FALSE, notch = TRUE, names = Gene)

summary(Skew$FemaleMaturityDays)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 60.0   365.0   365.0   634.3   730.0  4015.0 

SkewShort = Skew[Skew$FemaleMaturityDays <= 365.0,];  SkewShort$GT = 'Short'
SkewLong = Skew[Skew$FemaleMaturityDays >= 730.0,];  SkewLong$GT = 'Long'
SkewExtremes = rbind(SkewShort,SkewLong)

SkewShort = Skew[Skew$FemaleMaturityDays <= 365.0,];  SkewShort$GT = 'Short'
SkewLong = Skew[Skew$FemaleMaturityDays >= 730.0,];  SkewLong$GT = 'Long'
SkewExtremes = rbind(SkewShort,SkewLong)

pdf('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/4_FIGURES/AtGcSkewAves.pdf', width = 14)
par(cex.axis=0.7)
par(mfrow = c(2,1))

par(mfrow = c(2,1))
boxplot(Skew$AtSkew ~ Skew$Timing, outline = FALSE, notch = TRUE, ylim = c(-0.1,1), main = 'AT skew for genes with different time being single stranded (All Aves)', names = Gene)
boxplot(Skew$GcSkew ~ Skew$Timing, outline = FALSE, notch = TRUE, ylim = c(-1,0.1), main = 'GC skew for genes with different time being single stranded (All Aves)', names = Gene)

GeneVec = c('L.Cox1','S.Cox1','L.Cox2','S.Cox2','L.Atp8','S.Atp8','L.Atp6','S.Atp6','L.Cox3','S.Cox3','L.Nd3','S.Nd3','L.Nd4l','S.Nd4l','L.Nd4','S.Nd4','L.Nd5','S.Nd5','L.Cytb','S.Cytb','L.Nd1','S.Nd1','L.Nd2','S.Nd2')
par(mfrow = c(2,1))
boxplot(SkewExtremes$AtSkew ~ SkewExtremes$GT*SkewExtremes$Timing, outline = FALSE, notch = TRUE, ylim = c(-0.1,1), col = c('green','red'), main = paste('AT skew for genes with different time being single stranded','(Red - Short GT < 365 days; Green - Long GT: > 730 days)',sep = '\n'), names = GeneVec)
abline(h = 0, lwd = 2, lt = 2, col = 'red')
boxplot(SkewExtremes$GcSkew ~ SkewExtremes$GT*SkewExtremes$Timing, outline = FALSE, notch = TRUE, ylim = c(-1,0.1), col = c('green','red'), main = paste('GC skew for genes with different time being single stranded','(Red - Short GT < 365 days; Green - Long GT: > 730 days)',sep = '\n'), names = GeneVec)
abline(h = 0, lwd = 2, lt = 2, col = 'red')
dev.off()

################### Aves
rm(list=ls(all=TRUE))
Skew = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/cds_at_gc_skew.txt', header = TRUE)
names(Skew) = c('Species','Gene','AtSkew','GcSkew')

AnAge = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/1_RAW_DATA/anage_data.txt', header = TRUE, sep = '\t')
table(AnAge$Class) # 
#Actinopterygii           Amphibia               Aves Cephalaspidomorphi     Chondrichthyes           Mammalia           Reptilia      Sarcopterygii 
#824                174               1191                 16                116               1327                544                  4 
AnAge = AnAge[AnAge$Class == 'Reptilia' & ! is.na(AnAge$Female.maturity..days.),]; 
AnAge$Species = paste(AnAge$Genus,AnAge$Species, sep='_') 
AnAge = data.frame(AnAge$Species,AnAge$Female.maturity..days); names(AnAge)=c('Species','FemaleMaturityDays')

Skew = merge(Skew,AnAge, by = 'Species')
table(Skew$Gene) # ATP6 ATP8 COX1 COX2 COX3 CytB  ND1  ND2  ND3  ND4 ND4L  ND5  ND6

Skew = Skew[Skew$Gene %in% c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','CytB','ND1','ND2'),] # no ND6, ND1 and ND2 are located at the end approximately

Gene = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','CytB', 'ND1','ND2')
Timing = seq(1:12)
NewData = data.frame(Gene,Timing)

Skew = merge(Skew,NewData)
Skew$Timing = as.factor(Skew$Timing)
boxplot(Skew$AtSkew ~ Skew$Timing, outline = FALSE, notch = TRUE, names = Gene)
boxplot(Skew$GcSkew ~ Skew$Timing, outline = FALSE, notch = TRUE, names = Gene)

summary(Skew$FemaleMaturityDays)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 365    1277    2750    2586    3650    4560 

SkewShort = Skew[Skew$FemaleMaturityDays <= 2750.0,];  SkewShort$GT = 'Short'
SkewLong = Skew[Skew$FemaleMaturityDays > 2750.0,];  SkewLong$GT = 'Long'
SkewExtremes = rbind(SkewShort,SkewLong)

SkewShort = Skew[Skew$FemaleMaturityDays <= 2750.0,];  SkewShort$GT = 'Short'
SkewLong = Skew[Skew$FemaleMaturityDays > 2750.0,];  SkewLong$GT = 'Long'
SkewExtremes = rbind(SkewShort,SkewLong)

pdf('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/4_FIGURES/AtGcSkewReptilia.pdf', width = 14)
par(cex.axis=0.7)
par(mfrow = c(2,1))

par(mfrow = c(2,1))
boxplot(Skew$AtSkew ~ Skew$Timing, outline = FALSE, notch = TRUE, ylim = c(-0.1,1), main = 'AT skew for genes with different time being single stranded (All Aves)', names = Gene)
boxplot(Skew$GcSkew ~ Skew$Timing, outline = FALSE, notch = TRUE, ylim = c(-1,0.1), main = 'GC skew for genes with different time being single stranded (All Aves)', names = Gene)

GeneVec = c('L.Cox1','S.Cox1','L.Cox2','S.Cox2','L.Atp8','S.Atp8','L.Atp6','S.Atp6','L.Cox3','S.Cox3','L.Nd3','S.Nd3','L.Nd4l','S.Nd4l','L.Nd4','S.Nd4','L.Nd5','S.Nd5','L.Cytb','S.Cytb','L.Nd1','S.Nd1','L.Nd2','S.Nd2')
par(mfrow = c(2,1))
boxplot(SkewExtremes$AtSkew ~ SkewExtremes$GT*SkewExtremes$Timing, outline = FALSE, notch = TRUE, ylim = c(-0.1,1), col = c('green','red'), main = paste('AT skew for genes with different time being single stranded','(Red - Short GT < 2750 days; Green - Long GT: > 2750 days)',sep = '\n'), names = GeneVec)
abline(h = 0, lwd = 2, lt = 2, col = 'red')
boxplot(SkewExtremes$GcSkew ~ SkewExtremes$GT*SkewExtremes$Timing, outline = FALSE, notch = TRUE, ylim = c(-1,0.1), col = c('green','red'), main = paste('GC skew for genes with different time being single stranded','(Red - Short GT < 2750 days; Green - Long GT: > 2750 days)',sep = '\n'), names = GeneVec)
abline(h = 0, lwd = 2, lt = 2, col = 'red')
dev.off()

###################################
###### 20.02.2018: analysis of Mammals from Alina's file "/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/cds_at_gc_skew.txt": PAPER
###################################

rm(list=ls(all=TRUE))
# Skew = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/cds_at_gc_skew.txt', header = TRUE)
# Skew = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/cds_at_gc_skew_with_sites_numbers.txt', header = TRUE)
Skew = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/cds_at_gc_skew_with_sites_number_python.csv', header = TRUE)

names(Skew) = c('Parasha','Species','Gene','GcSkew','AtSkew')

#### want to remove species with small number of sites.... it is better now. BUT think more!!!
Temp = Skew[is.na(Skew$AtSkew) | is.na(Skew$GcSkew),]
NaVec = unique(Temp$Species); length(NaVec)  # 1505!!!
Skew = Skew[!Skew$Species %in% NaVec,]

#### want to remove species completely saturated (skew = 1)
Temp = Skew[abs(Skew$AtSkew) == 1 | abs(Skew$GcSkew) == 1,]
VecOfSaturatedSpecies = unique(Temp$Species); length(VecOfSaturatedSpecies)  # 1505!!!
Skew = Skew[!Skew$Species %in% VecOfSaturatedSpecies,]

##### GENERATION LENGTH FOR ALL MAMMALS
#  https://natureconservation.pensoft.net/article/1343/  data for generation length of ALL mammals!!!
GT = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/GenerationLenghtforAllMammals/GenerationLenghtforMammals.xlsx.txt', header = TRUE, sep = '\t')
head(GT)
str(GT)
GT$Species = gsub(' ','_',GT$Scientific_name)
length(unique(GT$Species))
summary(GT$AdultBodyMass_g)     # 2        21        71    136058       614 154321304 
summary(GT$GenerationLength_d)  # 129.0   624.4  1101.1  1578.8  2064.3 18980.0 # max = 18980 days => 52 years. ok
GT = GT[,c(11,13)]
summary(GT$GenerationLength_d)

#AnAge = read.table('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/COMPARATIVE_GENOMICS/1_RAW_DATA/anage_data.txt', header = TRUE, sep = '\t')
#head(AnAge)
# AnAge = AnAge[AnAge$Class == 'Mammalia' & ! is.na(AnAge$Female.maturity..days.),];
#AnAge = AnAge[AnAge$Class == 'Mammalia' & ! is.na(AnAge$Metabolic.rate..W.)  & ! is.na(AnAge$Female.maturity..days.) ,]; 
#AnAge$Species = paste(AnAge$Genus,AnAge$Species, sep='_') 
# AnAge = data.frame(AnAge$Species,AnAge$Female.maturity..days); names(AnAge)=c('Species','FemaleMaturityDays')

Skew = merge(Skew,GT, by = 'Species')
#Skew = merge(Skew,AnAge, by = 'Species')
table(Skew$Gene) # ATP6 ATP8 COX1 COX2 COX3 CytB  ND1  ND2  ND3  ND4 ND4L  ND5  ND6

############# SKEW FOR ND6

ND6 = Skew[Skew$Gene == 'ND6',]
ND6$AtSkew = ND6$AtSkew*(-1)
ND6$GcSkew = ND6$GcSkew*(-1)

Skew = Skew[Skew$Gene != 'ND6',]
Skew = rbind(Skew,ND6)

Skew = Skew[Skew$Gene %in% c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','ND6','CytB','ND1','ND2'),] # only genes with simple reconstruction of DSSH from short to long
# Skew = Skew[Skew$Gene %in% c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','ND6','CytB'),] # only genes with simple reconstruction of DSSH from short to long
# Skew = Skew[Skew$Gene %in% c('COX1','COX2','ATP6_8','COX3','ND3_4L','ND4','ND5','CytB'),] # only genes with simple reconstruction of DSSH from short to long

Gene = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','ND6','CytB', 'ND1','ND2')
# Gene = c('COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','CytB')
# Gene = c('COX1','COX2','ATP6_8','COX3','ND3_4L','ND4','ND5','CytB')
 
  Timing = seq(1:13)
# Timing = seq(1:8)
NewData = data.frame(Gene,Timing)
Skew = merge(Skew,NewData)
Skew$Timing = as.factor(Skew$Timing)
boxplot(Skew$AtSkew ~ Skew$Timing, outline = FALSE, notch = TRUE, names = Gene)
boxplot(Skew$GcSkew ~ Skew$Timing, outline = FALSE, notch = TRUE, names = Gene)

####### segments

summary(Skew$GenerationLength_d) # 1671 3650
#summary(Skew$Metabolic.rate..W.) # 1671 3650

VecOfLongLivedSpecies = unique(Skew[Skew$GenerationLength_d > quantile(Skew$GenerationLength_d, 0.5),]$Species); length(VecOfLongLivedSpecies)
VecOfShortLivedSpecies = unique(Skew[Skew$GenerationLength_d < quantile(Skew$GenerationLength_d, 0.5),]$Species); length(VecOfShortLivedSpecies)

AverageShortLived = Skew[Skew$Species %in% VecOfShortLivedSpecies,]
AverageShortLivedAgg1 = aggregate(AverageShortLived$AtSkew, by = list(AverageShortLived$Gene), FUN = mean); names(AverageShortLivedAgg1)=c('Gene','AtSkew')
AverageShortLivedAgg2 = aggregate(AverageShortLived$GcSkew, by = list(AverageShortLived$Gene), FUN = mean); names(AverageShortLivedAgg2)=c('Gene','GcSkew')
AverageShortLivedAgg = merge(AverageShortLivedAgg1, AverageShortLivedAgg2)
AverageShortLivedAgg = merge(AverageShortLivedAgg, NewData)
AverageShortLivedAgg = AverageShortLivedAgg[order(AverageShortLivedAgg$Timing),]

AverageLongLived = Skew[Skew$Species %in% VecOfLongLivedSpecies,]
AverageLongLivedAgg1 = aggregate(AverageLongLived$AtSkew, by = list(AverageLongLived$Gene), FUN = mean);  names(AverageLongLivedAgg1)=c('Gene','AtSkew')
AverageLongLivedAgg2 = aggregate(AverageLongLived$GcSkew, by = list(AverageLongLived$Gene), FUN = mean);  names(AverageLongLivedAgg2)=c('Gene','GcSkew')
AverageLongLivedAgg = merge(AverageLongLivedAgg1, AverageLongLivedAgg2)
AverageLongLivedAgg = merge(AverageLongLivedAgg, NewData)
AverageLongLivedAgg = AverageLongLivedAgg[order(AverageLongLivedAgg$Timing),]

pdf('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/4_FIGURES/AtGcSkew50Percentiles.pdf', width = 20, height = 14)
# plot(NA, xlim=c(1,10), ylim=c(-1,1), xlab='', ylab="skew", main = 'skew')
plot(NA, xlim=c(1,13), ylim=c(-1,1), xlab='', ylab="skew", main = 'skew')
for (i in 1:length(VecOfShortLivedSpecies))
{ # i = 1
  Temp = Skew[Skew$Species == VecOfShortLivedSpecies[i],]
  Temp = Temp[order(Temp$Timing),]
   if (nrow(Temp) == 13) # 10
  {
    for (count in 1:(nrow(Temp)-1))
    {
      segments(count, Temp$AtSkew[count], count+1, Temp$AtSkew[count+1], col = rgb(0.1,0.1,0.1,0.1), lwd = 3) # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
      segments(count, Temp$GcSkew[count], count+1, Temp$GcSkew[count+1], col = rgb(0.1,0.1,0.1,0.1), lwd = 3) # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
    }
  }
}
for (i in 1:length(VecOfLongLivedSpecies))
{ # i = 1
  Temp = Skew[Skew$Species == VecOfLongLivedSpecies[i],]
  Temp = Temp[order(Temp$Timing),]
  if (nrow(Temp) == 13) # 10
    {
    for (count in 1:(nrow(Temp)-1))
    {
    segments(count, Temp$AtSkew[count], count+1, Temp$AtSkew[count+1], col = rgb(1,0.1,0.1,0.1), lwd = 3) # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
    segments(count, Temp$GcSkew[count], count+1, Temp$GcSkew[count+1], col = rgb(1,0.1,0.1,0.1), lwd = 3) # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
    } 
  }
} 
### averages: 
Temp = AverageShortLivedAgg
for (count in 1:(nrow(Temp)-1))
{
  segments(count, Temp$AtSkew[count], count+1, Temp$AtSkew[count+1], col = rgb(0.1,0.1,0.1,1), lwd = 6) # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
  segments(count, Temp$GcSkew[count], count+1, Temp$GcSkew[count+1], col = rgb(0.1,0.1,0.1,1), lwd = 6) # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
}
Temp = AverageLongLivedAgg
for (count in 1:(nrow(Temp)-1))
{
  segments(count, Temp$AtSkew[count], count+1, Temp$AtSkew[count+1], col = rgb(1,0,0,1), lwd = 6) # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
  segments(count, Temp$GcSkew[count], count+1, Temp$GcSkew[count+1], col = rgb(1,0,0,1), lwd = 6) # rgb(red, green, blue, alpha, names = NULL, maxColorValue = 1)
}
dev.off()

#### 


#### boxplots

pdf('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/4_FIGURES/AtGcSkewMammals.pdf', width = 14)
par(cex.axis=0.7)
par(mfrow = c(2,1))
boxplot(Skew$AtSkew ~ Skew$Timing, outline = FALSE, notch = TRUE, ylim = c(-0.1,1), main = 'AT skew for 10 genes as a function of their Single Strand Time (All Mammals)', names = Gene)
boxplot(Skew$GcSkew ~ Skew$Timing, outline = FALSE, notch = TRUE, ylim = c(-1,0.1), main = 'GC skew for 10 genes as a function of their Single Strand Time (All Mammals)', names = Gene)

GeneVec = c('L.Cox1','S.Cox1','L.Cox2','S.Cox2','L.Atp8','S.Atp8','L.Atp6','S.Atp6','L.Cox3','S.Cox3','L.Nd3','S.Nd3','L.Nd4l','S.Nd4l','L.Nd4','S.Nd4','L.Nd5','S.Nd5','L.Cytb','S.Cytb','L.Nd1','S.Nd1','L.Nd2','S.Nd2')
par(mfrow = c(2,1))
boxplot(SkewExtremes$AtSkew ~ SkewExtremes$GT*SkewExtremes$Timing, outline = FALSE, notch = TRUE, ylim = c(-0.1,1), col = c('green','red'), main = paste('AT skew for genes with different time being single stranded','(Red - Short GT < 335 days; Green - Long GT: >1125 days)',sep = '\n'), names = GeneVec)
abline(h = 0, lwd = 2, lt = 2, col = 'red')
boxplot(SkewExtremes$GcSkew ~ SkewExtremes$GT*SkewExtremes$Timing, outline = FALSE, notch = TRUE, ylim = c(-1,0.1), col = c('green','red'), main = paste('GC skew for genes with different time being single stranded','(Red - Short GT < 335 days; Green - Long GT: >1125 days)',sep = '\n'), names = GeneVec)
abline(h = 0, lwd = 2, lt = 2, col = 'red')

### 10% instead of 25%
quantile(Skew$FemaleMaturityDays, 0.1) # 213
quantile(Skew$FemaleMaturityDays, 0.9) # 1826

SkewShort = Skew[Skew$FemaleMaturityDays <= 213,];  SkewShort$GT = 'Short'
SkewLong = Skew[Skew$FemaleMaturityDays >= 1826,];  SkewLong$GT = 'Long'
SkewExtremes = rbind(SkewShort,SkewLong)

SkewShort = Skew[Skew$FemaleMaturityDays <= 213,];  SkewShort$GT = 'Short'
SkewLong = Skew[Skew$FemaleMaturityDays >= 1826,];  SkewLong$GT = 'Long'
SkewExtremes = rbind(SkewShort,SkewLong)

GeneVec = c('L.Cox1','S.Cox1','L.Cox2','S.Cox2','L.Atp8','S.Atp8','L.Atp6','S.Atp6','L.Cox3','S.Cox3','L.Nd3','S.Nd3','L.Nd4l','S.Nd4l','L.Nd4','S.Nd4','L.Nd5','S.Nd5','L.Cytb','S.Cytb','L.Nd1','S.Nd1','L.Nd2','S.Nd2')
par(mfrow = c(2,1))
boxplot(SkewExtremes$AtSkew ~ SkewExtremes$GT*SkewExtremes$Timing, outline = FALSE, notch = TRUE, ylim = c(-0.1,1), col = c('green','red'), main = paste('AT skew for genes with different time being single stranded','(Red - Short GT < 213 days; Green - Long GT: >1826 days)',sep = '\n'), names = GeneVec)
abline(h = 0, lwd = 2, lt = 2, col = 'red')
boxplot(SkewExtremes$GcSkew ~ SkewExtremes$GT*SkewExtremes$Timing, outline = FALSE, notch = TRUE, ylim = c(-1,0.1), col = c('green','red'), main = paste('GC skew for genes with different time being single stranded','(Red - Short GT < 213 days; Green - Long GT: >1826 days)',sep = '\n'), names = GeneVec)
abline(h = 0, lwd = 2, lt = 2, col = 'red')

dev.off()

###### Genome-wide GC content (and Generation time) negatively correlate with genome-wide GC skew, and according to my analyses above (with GC skew from neutral positions) it is the same = GC is more negative for elephants than mice

AGG = aggregate(Skew$GcSkew, by = list(Skew$Species,Skew$FemaleMaturityDays), FUN = mean)
names(AGG) = c('Species','FemaleMaturityDays','MeanGcSkew')
cor.test(AGG$FemaleMaturityDays,AGG$MeanGcSkew, method = 'spearman') # negative significant => GC skew totally is higher (less negative) in short-lived mammals.
plot(AGG$FemaleMaturityDays,AGG$MeanGcSkew)


###################################
###### 1: Derive From MAtrix File paired node names: MoreDeepNode, MoreShallowNode
###################################

rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

# source("https://bioconductor.org/biocLite.R")
# biocLite("Biostrings")


library(Biostrings)
USER = 'KOSTYA';
# USER = 'ALYA';

if (USER == 'KOSTYA') {setwd('/home/kostya/konstantin/SCIENCE_PROJECTS_HEAD/MITOCHONDRIA/MutSpectrum/Results/re/');}

### READ THREE FILES FOR EACH GENE-SPECIES:

Terminal <- readDNAStringSet("Acanthogobius_hasta.ATP6.terminals.fa")
Terminal = data.frame(Terminal); names(Terminal) = c('TerminalSequence')
Terminal$TerminalNodeNumber = rownames(Terminal)
str(Terminal)

Ancestral  <- read.table("Acanthogobius_hasta.ATP6.ancestors.fa", header = FALSE, sep = ' '); names(Ancestral)=c('InternalNodeNumber', 'InternalSequence')
Ancestral$InternalNodeNumber = as.character(Ancestral$InternalNodeNumber)
Ancestral$InternalSequence = as.character(Ancestral$InternalSequence)
str(Ancestral)

Matrix <- read.table("Acanthogobius_hasta.ATP6.tree.matrix", header=TRUE, sep="\t")

### FIND ALL PAIRS OF ANCESTOR-DERIVED BRANCHES AND SAVE THEM AS MoreDeepNode AND MoreShallowNode

Final = data.frame('PleaseDeleteMe','PleaseDeleteMe'); names(Final) = c('MoreDeepNode','MoreShallowNode')
row.names(Matrix) = Matrix$X  # name each raw as first column 
Matrix = Matrix[-1]           # delete first column
OneBeforeLast = nrow(Matrix); 
TempMatrix = Matrix;
TempMatrix = TempMatrix[-OneBeforeLast,]; TempMatrix = TempMatrix[,-OneBeforeLast]
Max = nrow(TempMatrix); 
# TempMatrix is square matrix which we cut down step by step: 
for (OneBeforeLast in Max:1)
  { # OneBeforeLast = 6
  LINE = TempMatrix[OneBeforeLast,];
  for (j in 1:ncol(LINE))
    {#  j = 2
    if (LINE[j] == 1) {OneLine = data.frame(rownames(LINE),colnames(LINE)[j]); names(OneLine) = c('MoreDeepNode','MoreShallowNode'); Final = rbind(Final,OneLine);}
    }
  TempMatrix = TempMatrix[-OneBeforeLast,]; TempMatrix = TempMatrix[,-OneBeforeLast];
  IfOnlyTerminalBranchesLeft = TempMatrix[grep("RN_|OUTGRP",row.names(TempMatrix), invert=TRUE),];  # OUTGRP'
  if (nrow(IfOnlyTerminalBranchesLeft) == 0) 
      {break;}
  }
Final = Final[Final$MoreShallowNode != 'PleaseDeleteMe',]
Final$MoreShallowNode = gsub('X','',Final$MoreShallowNode)  # delete 'X'
Final$MoreDeepNode = gsub('X','',Final$MoreDeepNode)        # delete 'X'

Final$BranchPosition = 'Internal';
for (i in 1:nrow(Final))
  { # i = 1
  if (length(as.character(grep('RN_',Final$MoreShallowNode[i])))) {Final$BranchPosition[i] = 'External'}
  }

write.table(Final, 'FinalExample.txt', quote = FALSE, row.names = FALSE)

###################################
###### 2: print out pairwise alignments, corresponding to each pair of nodes
###################################

library("Biostrings")

###################################
###### 13.01.2018: Analysis of the privet_alya.txt: Tr/Tv for each species and correlation with ecology
###################################

rm(list=ls(all=TRUE))
# setwd('/home/kostya/konstantin/SCIENCE_PROJECTS_HEAD/MITOCHONDRIA/MutSpectrum/')
setwd('/hdd/SCIENCE_PROJECTS_BODY/MITOCHONDRIA/MutSpectrum/1_RAW/Results/')

MUT = read.table('privet_alya.txt', header = TRUE)
table(MUT$Subs)
VecOfNormalSubstitutions = c('A_C','C_A','A_G','G_A','C_G','G_C','C_T','T_C','G_T','T_G','T_A','A_T')
MUT = MUT[MUT$Subs %in% VecOfNormalSubstitutions,]

MUT = MUT[MUT$TAXON == 'Mammalia',]
# MUT = MUT[!is.na(MUT$Maximum.longevity..yrs.),]
MUT = MUT[MUT$AncestralAA == MUT$DescendantAA,]

### normalization by nucletoties in the third position of four-fold synonymous substitutions:
NORM = read.table('ATGC_counts_in_SYN_codons_wit_full_gene.txt', header = TRUE)
NORM$Gene = gsub("(.*)\\.",'',NORM$Species)
NORM$Species = gsub("\\.(.*)",'',NORM$Species)
MUT = merge(MUT,NORM, by = c("Species","Gene"))  # compare CountA.x and CountA.y  - they should be identical.

### GC content genome-wide for each species
#GC = MUT
#GC$GcContent = (GC$CountG + GC$CountC)/(GC$CountA + GC$CountT + GC$CountG + GC$CountC)
#GC = aggregate(GC$GcContent, by = list(GC$Species), FUN = mean)
#names(GC)=c('Species','GcContent')

### fourfold degenerate sites:
VecOfSynFourFoldDegenerateSites = c('CTT','CTC','CTA','CTG','GTT','GTC','GTA','GTG','TCT','TCC','TCA','TCG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','GCT','GCC','GCA','GCG','CGT','CGC','CGA','CGG','GGT','GGC','GGA','GGG'); length(unique(VecOfSynFourFoldDegenerateSites)) # 32
MUT = MUT[MUT$AncestorCodon %in% VecOfSynFourFoldDegenerateSites & MUT$DescendantCodon %in% VecOfSynFourFoldDegenerateSites,]

### 

table(MUT$Gene)
# ATP6  ATP8  COX1  COX2  COX3  CytB   ND1   ND2   ND3   ND4  ND4L   ND5   ND6 
# 2949   703  6798  3923  3265 25081  5302  5538  1796  8516  1661  8416    52
# MUT = MUT[MUT$Gene != 'ND6',]
# MUT = MUT[MUT$Gene == 'CytB',]

MUT$NumberOfSynMutPerSpecies = 1
AGG1 = aggregate(MUT$NumberOfSynMutPerSpecies, by = list(MUT$Species), FUN = sum)
summary(AGG1$x) # 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.0    30.0   100.0   239.5   285.0  2629.0 
ListOfSpeciesWithManySubst = AGG1[AGG1$x >= 30,]$Group.1; length(ListOfSpeciesWithManySubst) # 232
MUT = MUT[MUT$Species %in% ListOfSpeciesWithManySubst,]
NumberOfSubst = AGG1; names(NumberOfSubst) = c('Species','NumberOfAllSynonymousSubstitutions')
NumberOfSubst = NumberOfSubst[NumberOfSubst$NumberOfAllSynonymousSubstitutions >= 30,]

#### normalize by nucleotide count (miss this part if don't want to normilaze by nucleotide context) [results are very similar]
EXTRACT = function(x) {first = unlist(strsplit(as.character(x),'_'))[1]; return(first);}; MUT$AncestralNuc = apply(as.matrix(MUT$Subs), 1, EXTRACT)
MUT_A = MUT[MUT$AncestralNuc == 'A',]; MUT_T = MUT[MUT$AncestralNuc == 'T',]; MUT_G = MUT[MUT$AncestralNuc == 'G',]; MUT_C = MUT[MUT$AncestralNuc == 'C',]
MUT_A$NumberOfSynMutPerSpecies = MUT_A$NumberOfSynMutPerSpecies/MUT_A$CountA_Syn;
MUT_T$NumberOfSynMutPerSpecies = MUT_T$NumberOfSynMutPerSpecies/MUT_T$CountT_Syn;
MUT_G$NumberOfSynMutPerSpecies = MUT_G$NumberOfSynMutPerSpecies/MUT_G$CountG_Syn;
MUT_C$NumberOfSynMutPerSpecies = MUT_C$NumberOfSynMutPerSpecies/MUT_C$CountC_Syn;
MUT = rbind(MUT_A,MUT_T,MUT_G,MUT_C)

##### 
##### the strongest correlations:
# VecOfTransitionSubstitutions = c('C_T','T_C') vs VecOfTransversionSubstitutions = c('A_C','C_A','C_G','G_C','G_T','T_G','T_A','A_T'): rho = 0.37
# VecOfTransitionSubstitutions = c('T_C') vs VecOfTransversionSubstitutions = c('A_C','C_A','C_G','G_C','G_T','T_G','T_A','A_T'): rho = 0.38
# VecOfTransitionSubstitutions = c('T_C') vs VecOfTransversionSubstitutions = c('T_A','A_T'): rho = 0.40
# VecOfTransitionSubstitutions = c('T_C') vs VecOfTransversionSubstitutions = c('A_T'): rho = 0.39  

# VecOfTransitionSubstitutions = c('A_G','G_A','C_T','T_C') # all tr
# VecOfTransitionSubstitutions = c('G_A','T_C') # Ts, which depends on the time of beeng single-stranded 
VecOfTransitionSubstitutions = c('G_A','T_C') # they are driven by DNA polymerase?
# VecOfTransitionSubstitutions = c('T_C') # 

# VecOfTransversionSubstitutions = c('C_A','A_C','C_G','G_C','G_T','T_G','T_A','A_T') # ALL Tv 
# VecOfTransversionSubstitutions = c('C_T','A_G')
# VecOfTransversionSubstitutions = c('A_G') #
VecOfTransversionSubstitutions = c('C_A','G_T')  # Tv, induced by ROS
# VecOfTransversionSubstitutions = c('G_T')  # no correlation between tr/tv and GT
# VecOfTransversionSubstitutions = c('C_A')    # there is a correlation between tr/tv and GT: p-value = 0.0007141, rho = 0.255
# VecOfTransversionSubstitutions = c('A_C','C_G','G_C','T_G','T_A','A_T')  # Tv not induced by ROS
# VecOfTransversionSubstitutions = c('A_T') # ALL Tv
# VecOfTransversionSubstitutions = c('C_G','G_C','T_G') # ALL Tv

MUT_TR = MUT[MUT$Subs %in% VecOfTransitionSubstitutions,]
MUT_TV = MUT[MUT$Subs %in% VecOfTransversionSubstitutions,]

AGG1 = aggregate(MUT_TR$NumberOfSynMutPerSpecies, by = list(MUT_TR$Species), FUN = sum); names(AGG1) = c('Species','Tr')
AGG2 = aggregate(MUT_TV$NumberOfSynMutPerSpecies, by = list(MUT_TV$Species), FUN = sum); names(AGG2) = c('Species','Tv')
TrTv = merge(AGG1,AGG2, by = 'Species')
TrTv$TrTv = TrTv$Tr/TrTv$Tv

hist(TrTv$TrTv, breaks = 100)

ECO = unique(data.frame(MUT$Species, MUT$Female.maturity..days.)); names(ECO)=c('Species','GenerationTime')
ECO = ECO[!is.na(ECO$GenerationTime),] # 197!!! FUCK, WHY LESS THAN BEFORE!!!!
TrTv = merge(TrTv,ECO, by = 'Species')  # FUCK, AGAIN LESS!!!!
TrTv = merge(TrTv,NumberOfSubst, by = 'Species')
# TrTv = merge(TrTv,GC, by = 'Species')

A = cor.test(TrTv$TrTv,TrTv$GenerationTime, method = 'spearman') # positive!!!
# B = cor.test(TrTv$TrTv,TrTv$GcContent, method = 'spearman') # positive but week!!!
C<-lm(log2(TrTv$TrTv) ~ log2(TrTv$GenerationTime) + log2(TrTv$GcContent)); summary(C)  # GT is more significant than GC

PVal = round(as.numeric(A[3]),10)
Rho  = round(as.numeric(A[4]),5)
N = nrow(TrTv)
pdf('TrTvDirectVersusGenTimeFor188Mammals.pdf')
# pdf('TrTvNormalizedVersusGenTimeFor188Mammals.pdf')
title = paste('Spearman Rho = ',Rho,', Pval = ',PVal,', (N=',N,')', sep='')
plot(log2(TrTv$GenerationTime),log2(TrTv$TrTv), ylab = 'log2(Tr/Tv)', xlab = 'log2(generation time in days)', main = title)
dev.off()

write.table(TrTv, file = '204Species.AllSynon.NoNormalisationByNuclCount.txt', quote = FALSE, row.names = FALSE, sep = '\t')
  
###################################
###### 16.01.2018: Analysis of the privet_alya.txt: Three Nucleotide Motivs around Synon Four-Fold Substitutions
###################################

rm(list=ls(all=TRUE))
setwd('/home/kostya/konstantin/SCIENCE_PROJECTS_HEAD/MITOCHONDRIA/MutSpectrum/')
MUT = read.table('privet_alya.txt', header = TRUE)
table(MUT$Subs)
VecOfNormalSubstitutions = c('A_C','C_A','A_G','G_A','C_G','G_C','C_T','T_C','G_T','T_G','T_A','A_T')
MUT = MUT[MUT$Subs %in% VecOfNormalSubstitutions,]

MUT = MUT[MUT$TAXON == 'Mammalia',]
# MUT = MUT[!is.na(MUT$Maximum.longevity..yrs.),]
MUT = MUT[MUT$AncestralAA == MUT$DescendantAA,]
### fourfold degenerate sites:
VecOfSynFourFoldDegenerateSites = c('CTT','CTC','CTA','CTG','GTT','GTC','GTA','GTG','TCT','TCC','TCA','TCG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','GCT','GCC','GCA','GCG','CGT','CGC','CGA','CGG','GGT','GGC','GGA','GGG'); length(unique(VecOfSynFourFoldDegenerateSites)) # 32
MUT = MUT[MUT$AncestorCodon %in% VecOfSynFourFoldDegenerateSites & MUT$DescendantCodon %in% VecOfSynFourFoldDegenerateSites,]

# synonimous four fold substitutions occur only at the third position => to get three nucleotide motif we have to take second and third positions in codon and add first position from the next
EXTRACT1 = function(x) {temp = unlist(strsplit(as.character(x),''))[1]; return(temp);}; 
EXTRACT2 = function(x) {temp = unlist(strsplit(as.character(x),''))[2]; return(temp);};
EXTRACT3 = function(x) {temp = unlist(strsplit(as.character(x),''))[3]; return(temp);};

MUT$MotivBefore = paste(apply(as.matrix(MUT$AncestorCodon), 1, EXTRACT2),apply(as.matrix(MUT$AncestorCodon), 1, EXTRACT3),apply(as.matrix(MUT$NextAncCodon), 1, EXTRACT1), sep = '')
MUT$MotivAfter = paste(apply(as.matrix(MUT$DescendantCodon), 1, EXTRACT2),apply(as.matrix(MUT$DescendantCodon), 1, EXTRACT3),apply(as.matrix(MUT$NextDesCodon), 1, EXTRACT1), sep = '')
MUT$MutationWithContext = paste(MUT$MotivBefore,MUT$MotivAfter, sep = '_')
MUT$NeighborsBefore = paste(apply(as.matrix(MUT$AncestorCodon), 1, EXTRACT2),apply(as.matrix(MUT$NextAncCodon), 1, EXTRACT1), sep = '')
MUT$NeighborsAfter = paste(apply(as.matrix(MUT$DescendantCodon), 1, EXTRACT2),apply(as.matrix(MUT$NextDesCodon), 1, EXTRACT1), sep = '')
MUT = MUT[MUT$NeighborsBefore == MUT$NeighborsAfter,]
MUT = MUT[!grepl('-',MUT$NeighborsBefore),] # delete all neighbors with at least one '-'
MUT = MUT[!grepl('\\?',MUT$NeighborsBefore),] # delete all neighbors with at least one '?'
# MUT = MUT[MUT$Branch == 'External',]

MUT$NumberOfMut = 1
AGG = aggregate(MUT$NumberOfMut, by = list(MUT$MutationWithContext), FUN = sum)
names(AGG)=c('MutationWithContext','FreqAll')
AGG = AGG[order(AGG$Freq, decreasing = TRUE),]
AllMammals = AGG

### tha same for short- and long- lived mammals:
MUT = MUT[!is.na(MUT$Female.maturity..days.),]
summary(MUT$Female.maturity..days.) # 548
MUT1 = MUT[MUT$Female.maturity..days. < 548,]
AGG = aggregate(MUT1$NumberOfMut, by = list(MUT1$MutationWithContext), FUN = sum)
names(AGG)=c('MutationWithContext','FreqShortLived')
AGG = AGG[order(AGG$Freq, decreasing = TRUE),]
ShortLivedMammals = AGG

MUT1 = MUT[MUT$Female.maturity..days. >= 548,]
AGG = aggregate(MUT1$NumberOfMut, by = list(MUT1$MutationWithContext), FUN = sum)
names(AGG)=c('MutationWithContext','FreqLongLived')
AGG = AGG[order(AGG$Freq, decreasing = TRUE),]
LongLivedMammals = AGG

RES = merge(AllMammals,ShortLivedMammals, by = 'MutationWithContext', all = TRUE)
RES = merge(RES,LongLivedMammals, by = 'MutationWithContext', all = TRUE)
RES = RES[order(RES$FreqAll, decreasing = TRUE),]

EXTRACT1 = function(x) {temp = unlist(strsplit(as.character(x),''))[1]; return(temp);}; 
EXTRACT7 = function(x) {temp = unlist(strsplit(as.character(x),''))[7]; return(temp);};
RES$Context = paste(apply(as.matrix(RES$MutationWithContext), 1, EXTRACT1),apply(as.matrix(RES$MutationWithContext), 1, EXTRACT7), sep = '')

EXTRACT2 = function(x) {temp = unlist(strsplit(as.character(x),''))[2]; return(temp);}; 
EXTRACT6 = function(x) {temp = unlist(strsplit(as.character(x),''))[6]; return(temp);};
RES$Subst = paste(apply(as.matrix(RES$MutationWithContext), 1, EXTRACT2),apply(as.matrix(RES$MutationWithContext), 1, EXTRACT6), sep = '_')

VecOfTransitionSubstitutions = c('A_G','G_A','C_T','T_C') # all tr
VecOfTransversionSubstitutions = c('C_A','A_C','C_G','G_C','G_T','T_G','T_A','A_T') # ALL Tv

TR = RES[RES$Subst %in% VecOfTransitionSubstitutions,]; TR$MutType = 'Tr'
TV = RES[RES$Subst %in% VecOfTransversionSubstitutions,]; TV$MutType = 'Tv'
RES = rbind(TR,TV)
RES = RES[order(RES$FreqAll, decreasing = TRUE),]

write.table(RES, 'Motivs.txt', sep = '\t', quote = FALSE, row.names = FALSE)
# FromTorC = c('T_A','T_G','T_C','C_T','C_A','C_G') # 
FromAorG = c('A_T','A_G','A_C','G_T','G_A','G_C') # 

# ResFromTorC = RES[RES$Subst %in%  FromTorC,]
ResFromTorC = RES[RES$Subst %in%  FromAorG,]
ResFromTorC = ResFromTorC[!is.na(ResFromTorC$FreqLongLived),]
ResFromTorC = ResFromTorC[!is.na(ResFromTorC$FreqShortLived),]
ResFromTorC$FreqShortLived = ResFromTorC$FreqShortLived/(sum(ResFromTorC$FreqShortLived))
ResFromTorC$FreqLongLived = ResFromTorC$FreqLongLived/(sum(ResFromTorC$FreqLongLived))

Short = ResFromTorC[c(1,3,5,6,7)]; names(Short) = c('MutationWithContext' ,'Freq' ,'Context', 'Subst', 'MutType'); Short$Species = 'ShortLived'
Long = ResFromTorC[c(1,4,5,6,7)]; names(Long) = c('MutationWithContext' ,'Freq' ,'Context', 'Subst', 'MutType'); Long$Species = 'LongLived'
All = rbind(Short,Long)
All = All[order(All$Subst,All$Context,All$Species),]
pdf('Motiv.pdf', width=20, height=7) # pdf(width=0.5, height=0.3)
barplot(All$Freq,names = All$MutationWithContext, las = 2 , col = c('blue', 'red'))
dev.off()

###################################
###### 17.01.2018: Analysis of the privet_alya.txt: Three Nucleotide Motivs around Synon Four-Fold Substitutions
###################################

rm(list=ls(all=TRUE))
setwd('/home/kostya/konstantin/SCIENCE_PROJECTS_HEAD/MITOCHONDRIA/MutSpectrum/')
MUT = read.table('privet_alya.txt', header = TRUE)
table(MUT$Subs)
VecOfNormalSubstitutions = c('A_C','C_A','A_G','G_A','C_G','G_C','C_T','T_C','G_T','T_G','T_A','A_T')
MUT = MUT[MUT$Subs %in% VecOfNormalSubstitutions,]

MUT = MUT[MUT$TAXON == 'Mammalia',]
MUT = MUT[MUT$AncestralAA == MUT$DescendantAA,]
### fourfold degenerate sites:
VecOfSynFourFoldDegenerateSites = c('CTT','CTC','CTA','CTG','GTT','GTC','GTA','GTG','TCT','TCC','TCA','TCG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','GCT','GCC','GCA','GCG','CGT','CGC','CGA','CGG','GGT','GGC','GGA','GGG'); length(unique(VecOfSynFourFoldDegenerateSites)) # 32
MUT = MUT[MUT$AncestorCodon %in% VecOfSynFourFoldDegenerateSites & MUT$DescendantCodon %in% VecOfSynFourFoldDegenerateSites,]
# if the list of fourfild generate sites there is no A at the second place!!!! Consequently, there are no motfs, which start with A!!!!

###### longevity!!!!!
MUT = MUT[!is.na(MUT$Maximum.longevity..yrs.),]  # 
MUT = MUT[MUT$Maximum.longevity..yrs. <= median(MUT$Maximum.longevity..yrs.),]


# synonimous four fold substitutions occur only at the third position => to get three nucleotide motif we have to take second and third positions in codon and add first position from the next
EXTRACT1 = function(x) {temp = unlist(strsplit(as.character(x),''))[1]; return(temp);}; 
EXTRACT2 = function(x) {temp = unlist(strsplit(as.character(x),''))[2]; return(temp);};
EXTRACT3 = function(x) {temp = unlist(strsplit(as.character(x),''))[3]; return(temp);};

MUT$MotivBefore = paste(apply(as.matrix(MUT$AncestorCodon), 1, EXTRACT2),apply(as.matrix(MUT$AncestorCodon), 1, EXTRACT3),apply(as.matrix(MUT$NextAncCodon), 1, EXTRACT1), sep = '')
MUT$MotivAfter = paste(apply(as.matrix(MUT$DescendantCodon), 1, EXTRACT2),apply(as.matrix(MUT$DescendantCodon), 1, EXTRACT3),apply(as.matrix(MUT$NextDesCodon), 1, EXTRACT1), sep = '')
table(MUT$MotivBefore)  # there are no Motivs started from A!!!!!! FUCK!!!!!!
MUT$MutationWithContext = paste(MUT$MotivBefore,MUT$MotivAfter, sep = '_')
MUT$NeighborsBefore = paste(apply(as.matrix(MUT$AncestorCodon), 1, EXTRACT2),apply(as.matrix(MUT$NextAncCodon), 1, EXTRACT1), sep = '')
MUT$NeighborsAfter = paste(apply(as.matrix(MUT$DescendantCodon), 1, EXTRACT2),apply(as.matrix(MUT$NextDesCodon), 1, EXTRACT1), sep = '')
MUT = MUT[MUT$NeighborsBefore == MUT$NeighborsAfter,]
MUT = MUT[!grepl('-',MUT$NeighborsBefore),] # delete all neighbors with at least one '-'
MUT = MUT[!grepl('\\?',MUT$NeighborsBefore),] # delete all neighbors with at least one '?'
# MUT = MUT[MUT$Branch == 'External',]

MUT$NumberOfMut = 1
AGG = aggregate(MUT$NumberOfMut, by = list(MUT$MutationWithContext), FUN = sum)
names(AGG)=c('MutationWithContext','FreqAll')
AGG = AGG[order(AGG$Freq, decreasing = TRUE),]
EXTRACT2 = function(x) {temp = unlist(strsplit(as.character(x),''))[2]; return(temp);}; 
EXTRACT6 = function(x) {temp = unlist(strsplit(as.character(x),''))[6]; return(temp);};
AGG$Subst = paste(apply(as.matrix(AGG$MutationWithContext), 1, EXTRACT2),apply(as.matrix(AGG$MutationWithContext), 1, EXTRACT6), sep = '_')

EXTRACT1_COMPL = function(x) {temp = unlist(strsplit(as.character(x),''))[1]; if (temp == 'A') {res = 'T'}; if (temp == 'T') {res = 'A'}; if (temp == 'G') {res = 'C'}; if (temp == 'C') {res = 'G'}; return(res);}; 
EXTRACT2_COMPL = function(x) {temp = unlist(strsplit(as.character(x),''))[2]; if (temp == 'A') {res = 'T'}; if (temp == 'T') {res = 'A'}; if (temp == 'G') {res = 'C'}; if (temp == 'C') {res = 'G'}; return(res);};
EXTRACT3_COMPL = function(x) {temp = unlist(strsplit(as.character(x),''))[3]; if (temp == 'A') {res = 'T'}; if (temp == 'T') {res = 'A'}; if (temp == 'G') {res = 'C'}; if (temp == 'C') {res = 'G'}; return(res);}; 
EXTRACT5_COMPL = function(x) {temp = unlist(strsplit(as.character(x),''))[5]; if (temp == 'A') {res = 'T'}; if (temp == 'T') {res = 'A'}; if (temp == 'G') {res = 'C'}; if (temp == 'C') {res = 'G'}; return(res);}; 
EXTRACT6_COMPL = function(x) {temp = unlist(strsplit(as.character(x),''))[6]; if (temp == 'A') {res = 'T'}; if (temp == 'T') {res = 'A'}; if (temp == 'G') {res = 'C'}; if (temp == 'C') {res = 'G'}; return(res);}; 
EXTRACT7_COMPL = function(x) {temp = unlist(strsplit(as.character(x),''))[7]; if (temp == 'A') {res = 'T'}; if (temp == 'T') {res = 'A'}; if (temp == 'G') {res = 'C'}; if (temp == 'C') {res = 'G'}; return(res);}; 
AGG$ReverseMutationWithContext = paste(apply(as.matrix(AGG$MutationWithContext), 1, EXTRACT3_COMPL),apply(as.matrix(AGG$MutationWithContext), 1, EXTRACT2_COMPL),apply(as.matrix(AGG$MutationWithContext), 1, EXTRACT1_COMPL),'_',apply(as.matrix(AGG$MutationWithContext), 1, EXTRACT7_COMPL),apply(as.matrix(AGG$MutationWithContext), 1, EXTRACT6_COMPL),apply(as.matrix(AGG$MutationWithContext), 1, EXTRACT5_COMPL), sep = '')

AGG$ReverseSubst = paste(apply(as.matrix(AGG$ReverseMutationWithContext), 1, EXTRACT2),apply(as.matrix(AGG$ReverseMutationWithContext), 1, EXTRACT6), sep = '_')
TEMP = AGG[c(1,2)]; names(TEMP) = c('ReverseMutationWithContext','ReverseFreqAll') # 144

AllLong = merge(AGG,TEMP, by = 'ReverseMutationWithContext', all = TRUE) # the are many NA - empty places, altogether 180 lines (should be 192 = 96 + 96)

AllShort = merge(AGG,TEMP, by = 'ReverseMutationWithContext', all = FALSE) # the are many NA - empty places, altogether 180 lines (should be 192 = 96 + 96)
FromTorC = c('T_A','T_G','T_C','C_T','C_A','C_G') # 

write.table(AllShort, 'AllShort.txt')

ShortFromTorC = AllShort[AllShort$Subst %in% FromTorC,]

AGG1 = aggregate(ShortFromTorC$FreqAll, by = list(ShortFromTorC$Subst), FUN = sum)
AGG2 = aggregate(ShortFromTorC$ReverseFreqAll, by = list(ShortFromTorC$ReverseSubst), FUN = sum)
AGG1
#     C_A  1539
#     C_G   297
#     C_T 10860
#     T_A  1479
#     T_C 11420
#     T_G   262
AGG2
#     A_C   823
#     A_G  5454
#     A_T  1256
#     G_A 10057
#     G_C   322
#     G_T   281



#################################### Weighted substitutions:

require(wCorr)
> x<-c(8.944444444,6.238095238,20.2,15.66666667,3.823529412,7.833333333,11.5,16.66666667,13.28846154,17,14.44444444,14.7,24.4,16.625,23.6,9,25.53333333,17.75,19.08196721,24.93548387,16.87878788,64,11.5,15,7.333333333,17.15384615,24,23.16666667,27.48571429,20.94339623,18.83333333,21.73529412,12.26086957,29.55555556,1.851851852,26.25,16.2,18.50549451,24.90243902,11.78740157,47.57142857,31.5,5.909090909,7.9,29,14,8.857142857,19.6,7.6,7,6.5,5.888888889,3.984848485,12.22222222,12.6,5.166666667,47.75,19,30.58,19.83333333,17,14.66666667,33.73333333,47,71.5,17.66666667,22.75,36,17.66666667,37,13.8,4.526315789,26,11.75,19.75,28.5,17.16666667,8.105263158,4.115384615,42,70.5,46,16.8,5.484848485,7.9375,10.25,16.5,15.57142857,13.46153846,9.888888889,51,16.85714286,23,5.214285714,9.555555556,13.42857143,10.72857143,4.8,25,35,7.565217391,27.45783133,20.18918919,146,20.875,21.66666667,19.5,29.5,9.2,44,36,54,10.17460317,4.125,29.38461538,15.90740741,3.901408451,10.66666667,29,28.85714286,31.75,8.85046729,12.25,9.769230769,11.66666667,9.350393701,17.57142857,16.9,10.3,9.571428571,13,10.46153846,10.125,16.75,21.5,6.352941176,7.266666667,10,9.133333333,15.5,22.16666667,20,101.3333333,11.2,5.727272727,19.66666667,12.88372093,38.28571429,34.14285714,20.52777778,24,16.16666667,17.42105263,23.22916667,15.8,20.28571429,9,6.913043478,17.66666667,30,11,24.35294118,28.2,8.1875,62.64705882,21.23809524,24.94117647,11.75,18.65,94,19.5,20,8.303571429,7.25,6.395348837,50.5,12.75,2.269230769,16.625,15.82352941,32.85714286,25.48,61,11.34615385,13.4469697,22.33333333,13.3,8.4,14.86956522,9.375,7.466666667,4.333333333,10.74038462,1.818181818,27.52941176,34.3,27.25,35.27272727,19,36.13636364,12.13793103,32.7,21.73333333,3.857142857)
> y<-c(58,2192,550,730,751,1167,1475,335,76,1825,289,1461,304,8212,2701,912,730,738,730,548,502,730,730,1204,1095,1278,1095,274,669,406,797,413,258,639,66,1310,973,1673,852,365,1034,730,730,1461,304,411,75,1095,365,730,213,365,315,882,285,184,3287,708,914,1157,1009,253,956,650,646,1780,304,240,730,730,476,646,2470,2829,710,880,213,301,228,768,183,2555,354,122,327,771,898,800,365,595,600,595,1365,308,228,236,266,380,4018,365,1186,1238,1231,1125,1490,213,1277,365,548,365,456,1643,365,91,3287,243,40,122,540,639,272,42,95,162,120,47,62,210,502,1095,1460,730,106,2190,334,90,304,578,365,347,309,1826,3780,548,730,639,403,548,3194,1095,730,937,1268,3376,1514,329,91,547,487,1279,1095,365,1096,400,365,365,1186,912,730,912,1460,662,90,388,365,674,1003,730,296,1162,198,2983,2834,36,334,2190,182,365,350,1461,1460,504,496,315,2831,345,1278,1313,1734,1095,605,304,304,274)
> w<-c(179,152,106,100,246,53,125,106,743,522,417,157,127,141,738,70,398,300,1225,804,590,65,150,32,100,236,100,145,997,1163,119,773,305,275,77,109,86,1775,1062,1624,340,65,76,178,30,30,69,103,43,80,480,186,329,119,204,37,195,360,1579,250,108,47,521,48,145,56,285,74,56,152,74,105,135,204,83,177,109,346,133,43,143,235,178,214,572,45,35,116,188,98,104,250,96,174,1140,404,821,58,364,108,197,2362,784,147,175,136,41,61,51,45,518,55,704,41,395,913,348,385,90,627,131,2108,159,140,304,2629,130,358,113,74,42,149,89,71,45,125,124,132,152,396,1390,63,307,61,74,62,597,825,246,775,150,412,350,1163,84,894,90,364,56,93,36,431,146,147,1082,467,441,357,1179,95,41,126,521,33,318,103,880,85,141,286,237,662,62,321,1907,280,429,47,730,83,127,32,1221,31,485,353,339,798,60,817,381,337,341,34)
> weightedCorr(x, y, method = c("Spearman"), weights = w, ML = FALSE, fast = TRUE)
[1] 0.5247036
> weightedCorr(x, y, method = c("Pearson"), weights = w, ML = FALSE, fast = TRUE)
[1] 0.301082
> lm(y~x)

Call:
  lm(formula = y ~ x)

Coefficients:
  (Intercept)            x  
612.95        12.06  

> lm(y~x,weights=w)

Call:
  lm(formula = y ~ x, weights = w)

Coefficients:
  (Intercept)            x  
430.85        18.75  

> require(weights)
> wtd.cor(x,y,weight=w)
correlation    std.err  t.value      p.value
Y    0.301082 0.06709494 4.487403 1.209575e-05

> require(directlabels)
> require(ggplot2)
> my_data <- data.frame(x,y,w)
> c<-coef(lm(y~x))
> v<-coef(lm(y~x,weights=w))
> c
(Intercept)           x 
612.95219    12.06348 
> v
(Intercept)           x 
430.85083    18.75158 
> p <- ggplot(my_data, aes(x = x, y = y, size = w)) + geom_point(shape = 21, show_guide=FALSE) + geom_abline(intercept = c[1], slope = c[2], colour="red") + geom_abline(intercept = v[1], slope = v[2], colour="blue") + theme_bw()
> p




###################################
###### 22.01.2018: Analysis of the privet_alya.txt: 12 substitutions in 4-fold, normilazed by frequency of the first one
###################################

rm(list=ls(all=TRUE))
setwd('/home/kostya/konstantin/SCIENCE_PROJECTS_HEAD/MITOCHONDRIA/MutSpectrum/')
MUT = read.table('privet_alya.txt', header = TRUE)
table(MUT$Subs)
VecOfNormalSubstitutions = c('A_C','C_A','A_G','G_A','C_G','G_C','C_T','T_C','G_T','T_G','T_A','A_T')
MUT = MUT[MUT$Subs %in% VecOfNormalSubstitutions,]

VecOfTransitionSubstitutions = c('A_G','G_A','C_T','T_C') # all Tr
VecOfTransversionSubstitutions = c('C_A','A_C','C_G','G_C','G_T','T_G','T_A','A_T') # all Tv

###### FILTER:  synonymous
MUT = MUT[MUT$AncestralAA == MUT$DescendantAA,]

###### FILTER:  taxa  
MUT = MUT[MUT$TAXON == 'Mammalia',]

###### FILTER:  gene
MUT = MUT[MUT$Gene == 'CytB',]

###### FILTER: fourfold degenerate sites:
VecOfSynFourFoldDegenerateSites = c('CTT','CTC','CTA','CTG','GTT','GTC','GTA','GTG','TCT','TCC','TCA','TCG','CCT','CCC','CCA','CCG','ACT','ACC','ACA','ACG','GCT','GCC','GCA','GCG','CGT','CGC','CGA','CGG','GGT','GGC','GGA','GGG'); length(unique(VecOfSynFourFoldDegenerateSites)) # 32
MUT = MUT[MUT$AncestorCodon %in% VecOfSynFourFoldDegenerateSites & MUT$DescendantCodon %in% VecOfSynFourFoldDegenerateSites,]
# if the list of fourfild generate sites there is no A at the second place!!!! Consequently, there are no motfs, which start with A!!!!

###### FILTER: longevity: 
MUT = MUT[!is.na(MUT$Female.maturity..days.),]  # 
MUT = MUT[MUT$Female.maturity..days. >= quantile(MUT$Female.maturity..days.,0.75),]

###### FILTER: external branches or all
# MUT = MUT[MUT$Branch == 'External',]

####### count NumberOfMut, normilazed by the number of the first nucleotide (I want to derive a probability that a given nucleotide will mutate assuming that A = T = G = C = 25%)
##### but now it is not exactly what I want
MUT$NumberOfMut = 1
MUT$AncestralNuc = gsub("_(.*)",'',MUT$Subs)
MUT$DerivedNuc = gsub("(.*)_",'',MUT$Subs)
MUT_A = MUT[MUT$AncestralNuc == 'A',]; nrow(MUT_A) # 10197
MUT_T = MUT[MUT$AncestralNuc == 'T',]; nrow(MUT_T) # 16386
MUT_G = MUT[MUT$AncestralNuc == 'G',]; nrow(MUT_G) # 14367
MUT_C = MUT[MUT$AncestralNuc == 'C',]; nrow(MUT_C) # 16267

## G < C It can reflect that G are rare on L chain.
## But A < T, however we know that A more frequent on L chain. Mut bias?

#### Is this formula correct => the higher the number, the higher the probability to mutate!!!
MUT_A$NumberOfMut = MUT_A$NumberOfMut/MUT_A$PercentageA; summary(MUT_A$NumberOfMut)
MUT_T$NumberOfMut = MUT_T$NumberOfMut/MUT_T$PercentageT; summary(MUT_T$NumberOfMut)
MUT_G$NumberOfMut = MUT_G$NumberOfMut/MUT_G$PercentageG; summary(MUT_G$NumberOfMut)
MUT_C$NumberOfMut = MUT_C$NumberOfMut/MUT_C$PercentageC; summary(MUT_C$NumberOfMut)

#### rbind again 
nrow(MUT)  # 17520
MUT = rbind(MUT_A,MUT_T,MUT_G,MUT_C)
nrow(MUT)  # 17520

AGG = aggregate(MUT$NumberOfMut, by = list(MUT$Subs), FUN = sum)
names(AGG)=c('Subs','Freq')
AGG = AGG[order(AGG$Freq, decreasing = TRUE),]

G_A_T_C = c('G_A','T_C')

## on CYTB and external branches the assymetry is increasing (either because of CYTB (the longestreplication time), or the most parsimonious subset of mutations): 
##  (CYTB, external, short lived)
7   G_A 4218.14354  # 4218.14354 / (1731.24796+4218.14354) = 71%
11  T_C 1880.64310  # 1880.64310 / (512.84337 + 1880.64310) = 78%
6   C_T 1731.24796
2   A_G  512.84337
4   C_A  474.33865
10  T_A  356.55324
3   A_T  231.87422
1   A_C  213.86528
8   G_C  161.34057
9   G_T  160.23656
12  T_G   63.08326
5   C_G   57.27102
Tr_Tv = sum(AGG[AGG$Subs %in% VecOfTransitionSubstitutions,]$Freq)/sum(AGG[AGG$Subs %in% VecOfTransversionSubstitutions,]$Freq) # 4.85
Tr_Tv = sum(AGG[AGG$Subs == 'G_A',]$Freq)/sum(AGG[AGG$Subs %in% VecOfTransversionSubstitutions,]$Freq) # 2.45
Tr_Tv = sum(AGG[AGG$Subs == 'T_C',]$Freq)/sum(AGG[AGG$Subs %in% VecOfTransversionSubstitutions,]$Freq) # 1.09
Tr_Tv = sum(AGG[AGG$Subs %in% G_A_T_C,]$Freq)/sum(AGG[AGG$Subs %in% VecOfTransversionSubstitutions,]$Freq) # 3.548771
Tr_Tv = sum(AGG[AGG$Subs == 'C_T',]$Freq)/sum(AGG[AGG$Subs %in% VecOfTransversionSubstitutions,]$Freq) # 1.007381
Tr_Tv = sum(AGG[AGG$Subs == 'A_G',]$Freq)/sum(AGG[AGG$Subs %in% VecOfTransversionSubstitutions,]$Freq) # 0.2984141

##  (CYTB, external, long lived)
7   G_A 4805.73454  # 4805.73454  / (1747.76804 + 4805.73454) = 73%
11  T_C 1956.84914  # 1956.84914 / (532.44091 + 1956.84914) = 78%
6   C_T 1747.76804
2   A_G  532.44091
4   C_A  387.91641
10  T_A  296.35918
8   G_C  195.97444
1   A_C  181.82652
9   G_T  155.34979
3   A_T  138.49975
5   C_G   96.22694
12  T_G   39.16844
Tr_Tv = sum(AGG[AGG$Subs %in% VecOfTransitionSubstitutions,]$Freq)/sum(AGG[AGG$Subs %in% VecOfTransversionSubstitutions,]$Freq) # 6.06 / (4.85 + 6.06) = 55%
Tr_Tv = sum(AGG[AGG$Subs == 'G_A',]$Freq)/sum(AGG[AGG$Subs %in% VecOfTransversionSubstitutions,]$Freq) # 3.22 / (3.22 + 2.45) 56.7%
Tr_Tv = sum(AGG[AGG$Subs == 'T_C',]$Freq)/sum(AGG[AGG$Subs %in% VecOfTransversionSubstitutions,]$Freq) # 1.31 / (1.31 + 1.09) 54.5%
Tr_Tv = sum(AGG[AGG$Subs %in% G_A_T_C,]$Freq)/sum(AGG[AGG$Subs %in% VecOfTransversionSubstitutions,]$Freq) # 4.53 / (4.53+3.548771) = 56%
Tr_Tv = sum(AGG[AGG$Subs == 'C_T',]$Freq)/sum(AGG[AGG$Subs %in% VecOfTransversionSubstitutions,]$Freq) # 1.171959 / (1.171959 + 1.007381) = 53.7%
Tr_Tv = sum(AGG[AGG$Subs == 'A_G',]$Freq)/sum(AGG[AGG$Subs %in% VecOfTransversionSubstitutions,]$Freq) # 0.3570263 / (0.2984141 + 0.3570263) = 54.4%

So, the difference between long and short-lived mammals is driven by both G_A and T_C as compared to transversions, but the strongest effect is coming from G_A

# LEVEL OF ASSYMETRY BETWEEN SHORT AND LONG- LIVED MAMMALS IS VERY SIMILAR!?
## LONGLIVED level of assymetry: G_A/(C_T+G_A)
## LONGLIVED level of assymetry: T_C/(A_G+T_C)
## SHORTLIVED level of assymetry: G_A/(C_T+G_A)
## SHORTLIVED level of assymetry: T_C/(A_G+T_C)
### - symmetry is the result of ML approach?? PAUP!!!
## 2442.42116/(2442.42116 + 1385.42880)
## 3149.00169 / (3149.00169 + 1563.12330)





