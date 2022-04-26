rm(list=ls(all=TRUE))

library("rsample") # install.packages("rsample"); install.packages("utf8")

#################### 1: SLOPES FROM THE PAPER (Only Major Arc):
# FIG 2:
ChThSlopeMus = 7.42*10^-8
AhGhSlopeMus = 1.23*10^-8
# FIG 5A,B:
ChThSlopeHomo = 4.86*10^-7
AhGhSlopeHomo = 1.84*10^-7

### Homo versus Mus absolute slopes (A>G is becoming much steeper with Age):

ChThSlopeHomo/ChThSlopeMus # 6.5
AhGhSlopeHomo/AhGhSlopeMus # 15

### Homo versus Mus relative slopes of AhGh/ChTh (A>G is more similar to C>T in homo vs mus:)
AhGhSlopeHomo/ChThSlopeHomo # 0.37
AhGhSlopeMus/ChThSlopeMus   # 0.16
# 0.37/0.16 = 2


############################################################################
################### 2: COMPARE MUS YOUNG AND OLD:
############################################################################

################ MusYoungMajorArc
MusYoung = read.table("../../Body/1Raw/DuplexSeqData/DataS3.txt", sep = '\t', header = TRUE)
names(MusYoung)
MusYoung$AhGhFreq = MusYoung$T.C_Freq
MusYoung$ChThFreq = MusYoung$G.A_Freq

MusYoung$AllFreq = MusYoung$T.C_Freq+MusYoung$T.A_Freq + MusYoung$T.G_Freq + MusYoung$A.C_Freq+MusYoung$A.T_Freq + MusYoung$A.G_Freq + MusYoung$C.T_Freq+MusYoung$C.A_Freq + MusYoung$C.G_Freq + MusYoung$G.C_Freq+MusYoung$G.A_Freq + MusYoung$G.T_Freq
MusYoung$AhGhRelFreq = MusYoung$T.C_Freq/MusYoung$AllFreq
MusYoung$ChThRelFreq = MusYoung$G.A_Freq/MusYoung$AllFreq

MusYoung$AllCount = MusYoung$T.C_COUNT+MusYoung$T.A_COUNT + MusYoung$T.G_COUNT + MusYoung$A.C_COUNT+MusYoung$A.T_COUNT + MusYoung$A.G_COUNT + MusYoung$C.T_COUNT+MusYoung$C.A_COUNT + MusYoung$C.G_COUNT + MusYoung$G.C_COUNT+MusYoung$G.A_COUNT + MusYoung$G.T_Freq
MusYoung$AhGhCountFreq = MusYoung$T.C_COUNT/MusYoung$AllCount
MusYoung$ChThCountFreq = MusYoung$G.A_COUNT/MusYoung$AllCount

MusYoungMajorArc=MusYoung[MusYoung$BIN_NUM >= 53 & MusYoung$BIN_NUM <= 153,]
MusYoungMajorArc$MajorArcBins = seq(1:101)
plot(MusYoungMajorArc$MajorArcBins,MusYoungMajorArc$ChThFreq)
plot(MusYoungMajorArc$MajorArcBins,MusYoungMajorArc$AhGhFreq)
plot(MusYoungMajorArc$MajorArcBins,MusYoungMajorArc$AllFreq)
plot(MusYoungMajorArc$MajorArcBins,MusYoungMajorArc$AhGhRelFreq)

CytbYoungMouseAhGhRelFreq = mean(MusYoungMajorArc[MusYoungMajorArc$MajorArcBins>50,]$AhGhRelFreq) # 0.041
CytbYoungMouseChThRelFreq = mean(MusYoungMajorArc[MusYoungMajorArc$MajorArcBins>50,]$ChThRelFreq) # 0.56
CytbYoungMouseChThRelFreqPerBins = MusYoungMajorArc$ChThRelFreq
CytbYoungMouseAhGhRelFreqPerBins = MusYoungMajorArc$AhGhRelFreq

CytbYoungMouseAhGhCountFreq = mean(MusYoungMajorArc[MusYoungMajorArc$MajorArcBins>50,]$AhGhCountFreq) # 0.09
CytbYoungMouseChThCountFreq = mean(MusYoungMajorArc[MusYoungMajorArc$MajorArcBins>50,]$ChThCountFreq) # 0.43

MusYoungMajorArc.ChTh = lm(MusYoungMajorArc$ChThFreq ~ MusYoungMajorArc$MajorArcBins)
MusYoungMajorArc.AhGh = lm(MusYoungMajorArc$AhGhFreq ~ MusYoungMajorArc$MajorArcBins)

########### BOOTSTRAP OF YOUNG MICE: 
DataFrame = MusYoungMajorArc
VecOfChTh = c()
VecOfAhGh = c()
Rows = seq(1:nrow(DataFrame))
for (i in 1:1000)
{
  BootstrappedRows = sample(Rows, replace = TRUE) # bootstrapping
  BootstrappedDataFrame = DataFrame[BootstrappedRows,]
  BootstrappedDataFrame.ChTh = lm(BootstrappedDataFrame$ChThFreq ~ BootstrappedDataFrame$MajorArcBins)
  VecOfChTh = c(VecOfChTh,as.numeric(coefficients(BootstrappedDataFrame.ChTh)[2]))
  BootstrappedDataFrame.AhGh = lm(BootstrappedDataFrame$AhGhFreq ~ BootstrappedDataFrame$MajorArcBins)
  VecOfAhGh = c(VecOfAhGh,as.numeric(coefficients(BootstrappedDataFrame.AhGh)[2]))
}
VecOfYoungMiceChTh = VecOfChTh
VecOfYoungMiceAhGh = VecOfAhGh
wilcox.test(VecOfYoungMiceChTh/VecOfYoungMiceAhGh, mu = 1)
summary(VecOfYoungMiceChTh/VecOfYoungMiceAhGh)

################ MusOldMajorArc
MusOld = read.table("../../Body/1Raw/DuplexSeqData/DataS2.txt", sep = '\t', header = TRUE)
names(MusOld)
MusOld$AhGhFreq = MusOld$T.C_Freq
MusOld$ChThFreq = MusOld$G.A_Freq

MusOld$AllFreq = MusOld$T.C_Freq+MusOld$T.A_Freq + MusOld$T.G_Freq + MusOld$A.C_Freq+MusOld$A.T_Freq + MusOld$A.G_Freq + MusOld$C.T_Freq+MusOld$C.A_Freq + MusOld$C.G_Freq + MusOld$G.C_Freq+MusOld$G.A_Freq + MusOld$G.T_Freq
MusOld$AhGhRelFreq = MusOld$T.C_Freq/MusOld$AllFreq
MusOld$ChThRelFreq = MusOld$G.A_Freq/MusOld$AllFreq

MusOld$AllCount = MusOld$T.C_COUNT+MusOld$T.A_COUNT + MusOld$T.G_COUNT + MusOld$A.C_COUNT+MusOld$A.T_COUNT + MusOld$A.G_COUNT + MusOld$C.T_COUNT+MusOld$C.A_COUNT + MusOld$C.G_COUNT + MusOld$G.C_COUNT+MusOld$G.A_COUNT + MusOld$G.T_Freq
MusOld$AhGhCountFreq = MusOld$T.C_COUNT/MusOld$AllCount
MusOld$ChThCountFreq = MusOld$G.A_COUNT/MusOld$AllCount

MusOldMajorArc=MusOld[MusOld$BIN_NUM >= 53 & MusOld$BIN_NUM <= 153,]
MusOldMajorArc$MajorArcBins = seq(1:101)
plot(MusOldMajorArc$MajorArcBins,MusOldMajorArc$ChThFreq)
plot(MusOldMajorArc$MajorArcBins,MusOldMajorArc$AhGhFreq)

CytbOldMouseAhGhRelFreq = mean(MusOldMajorArc[MusOldMajorArc$MajorArcBins>50,]$AhGhRelFreq) # 0.068
CytbOldMouseChThRelFreq = mean(MusOldMajorArc[MusOldMajorArc$MajorArcBins>50,]$ChThRelFreq) # 0.68
CytbOldMouseChThRelFreqPerBins = MusOldMajorArc$ChThRelFreq
CytbOldMouseAhGhRelFreqPerBins = MusOldMajorArc$AhGhRelFreq

CytbOldMouseAhGhCountFreq = mean(MusOldMajorArc[MusOldMajorArc$MajorArcBins>50,]$AhGhCountFreq) # 0.1517
CytbOldMouseChThCountFreq = mean(MusOldMajorArc[MusOldMajorArc$MajorArcBins>50,]$ChThCountFreq) # 0.506

MusOldMajorArc.ChTh = lm(MusOldMajorArc$ChThFreq ~ MusOldMajorArc$MajorArcBins)
MusOldMajorArc.AhGh = lm(MusOldMajorArc$AhGhFreq ~ MusOldMajorArc$MajorArcBins)

########### BOOTSTRAP OF OLD MICE: 
DataFrame = MusOldMajorArc
VecOfChTh = c()
VecOfAhGh = c()
Rows = seq(1:nrow(DataFrame))
for (i in 1:1000)
{
  BootstrappedRows = sample(Rows, replace = TRUE) # bootstrapping
  BootstrappedDataFrame = DataFrame[BootstrappedRows,]
  BootstrappedDataFrame.ChTh = lm(BootstrappedDataFrame$ChThFreq ~ BootstrappedDataFrame$MajorArcBins)
  VecOfChTh = c(VecOfChTh,as.numeric(coefficients(BootstrappedDataFrame.ChTh)[2]))
  BootstrappedDataFrame.AhGh = lm(BootstrappedDataFrame$AhGhFreq ~ BootstrappedDataFrame$MajorArcBins)
  VecOfAhGh = c(VecOfAhGh,as.numeric(coefficients(BootstrappedDataFrame.AhGh)[2]))
}
VecOfOldMiceChTh = VecOfChTh
VecOfOldMiceAhGh = VecOfAhGh
wilcox.test(VecOfOldMiceChTh/VecOfOldMiceAhGh, mu = 1)
summary(VecOfOldMiceChTh/VecOfOldMiceAhGh)

###############: HUMAN DATA:
Homo = read.table("../../Body/1Raw/DuplexSeqData/DataS6.txt", sep = '\t', header = TRUE)
names(Homo)
Homo$AhGhFreq = Homo$T.C_Freq
Homo$ChThFreq = Homo$G.A_Freq

Homo$AllFreq = Homo$T.C_Freq+Homo$T.A_Freq + Homo$T.G_Freq + Homo$A.C_Freq+Homo$A.T_Freq + Homo$A.G_Freq + Homo$C.T_Freq+Homo$C.A_Freq + Homo$C.G_Freq + Homo$G.C_Freq+Homo$G.A_Freq + Homo$G.T_Freq
Homo$AhGhRelFreq = Homo$T.C_Freq/Homo$AllFreq
Homo$ChThRelFreq = Homo$G.A_Freq/Homo$AllFreq

Homo$AllCount = Homo$T.C_COUNT+Homo$T.A_COUNT + Homo$T.G_COUNT + Homo$A.C_COUNT+Homo$A.T_COUNT + Homo$A.G_COUNT + Homo$C.T_COUNT+Homo$C.A_COUNT + Homo$C.G_COUNT + Homo$G.C_COUNT+Homo$G.A_COUNT + Homo$G.T_Freq
Homo$AhGhCountFreq = Homo$T.C_COUNT/Homo$AllCount
Homo$ChThCountFreq = Homo$G.A_COUNT/Homo$AllCount

HomoMajorArc=Homo[Homo$BIN_NUM >= 28 & Homo$BIN_NUM <= 76,] # 53/2 == 26.5 > start from 28 153/2 = 76.5 > end at 76
HomoMajorArc$MajorArcBins = seq(1:nrow(HomoMajorArc))
plot(HomoMajorArc$MajorArcBins,HomoMajorArc$ChThFreq)
plot(HomoMajorArc$MajorArcBins,HomoMajorArc$AhGhFreq)

CytbHomoAhGhRelFreq = mean(HomoMajorArc[HomoMajorArc$MajorArcBins>25,]$AhGhRelFreq) # 0.14
CytbHomoChThRelFreq = mean(HomoMajorArc[HomoMajorArc$MajorArcBins>25,]$ChThRelFreq) # 0.44

CytbHomoChThRelFreqPerBins = HomoMajorArc$ChThRelFreq
CytbHomoAhGhRelFreqPerBins = HomoMajorArc$AhGhRelFreq

CytbHomoAhGhCountFreq = mean(HomoMajorArc[HomoMajorArc$MajorArcBins>25,]$AhGhCountFreq) # 0.206
CytbHomoChThCountFreq = mean(HomoMajorArc[HomoMajorArc$MajorArcBins>25,]$ChThCountFreq) # 0.28

HomoMajorArc.ChTh = lm(HomoMajorArc$ChThFreq ~ HomoMajorArc$MajorArcBins)
HomoMajorArc.AhGh = lm(HomoMajorArc$AhGhFreq ~ HomoMajorArc$MajorArcBins)

as.numeric(coefficients(HomoMajorArc.ChTh)[2])
as.numeric(coefficients(HomoMajorArc.AhGh)[2])

########### BOOTSTRAP OF HUMAN DATA: 
DataFrame = HomoMajorArc
VecOfChTh = c()
VecOfAhGh = c()
Rows = seq(1:nrow(DataFrame))
for (i in 1:1000)
{
  BootstrappedRows = sample(Rows, replace = TRUE) # bootstrapping
  BootstrappedDataFrame = DataFrame[BootstrappedRows,]
  BootstrappedDataFrame.ChTh = lm(BootstrappedDataFrame$ChThFreq ~ BootstrappedDataFrame$MajorArcBins)
  VecOfChTh = c(VecOfChTh,as.numeric(coefficients(BootstrappedDataFrame.ChTh)[2]))
  BootstrappedDataFrame.AhGh = lm(BootstrappedDataFrame$AhGhFreq ~ BootstrappedDataFrame$MajorArcBins)
  VecOfAhGh = c(VecOfAhGh,as.numeric(coefficients(BootstrappedDataFrame.AhGh)[2]))
}
VecOfHomoChTh = VecOfChTh
VecOfHomoAhGh = VecOfAhGh
wilcox.test(VecOfHomoChTh/VecOfHomoAhGh, mu = 1)
summary(VecOfHomoChTh/VecOfHomoAhGh)

#################################################################
############# 3: MAKE A TABLE WITH RESULTS OF 6 LINEAR MODELS (2 mutation types * 3 age groups)
#################################################################
GeneralMultiplicator = 10000000

summary(MusYoungMajorArc.ChTh)
round(as.numeric(coefficients(MusYoungMajorArc.ChTh)[1])*GeneralMultiplicator,3)
round(as.numeric(coefficients(MusYoungMajorArc.ChTh)[2])*GeneralMultiplicator,3)

summary(MusYoungMajorArc.AhGh)
round(as.numeric(coefficients(MusYoungMajorArc.AhGh)[1])*GeneralMultiplicator,3)
round(as.numeric(coefficients(MusYoungMajorArc.AhGh)[2])*GeneralMultiplicator,3)

summary(MusOldMajorArc.ChTh)
round(as.numeric(coefficients(MusOldMajorArc.ChTh)[1])*GeneralMultiplicator,3)
round(as.numeric(coefficients(MusOldMajorArc.ChTh)[2])*GeneralMultiplicator,3)

summary(MusOldMajorArc.AhGh)
round(as.numeric(coefficients(MusOldMajorArc.AhGh)[1])*GeneralMultiplicator,3)
round(as.numeric(coefficients(MusOldMajorArc.AhGh)[2])*GeneralMultiplicator,3)

summary(HomoMajorArc.ChTh)
round(as.numeric(coefficients(HomoMajorArc.ChTh)[1])*GeneralMultiplicator,3)
round(as.numeric(coefficients(HomoMajorArc.ChTh)[2])*GeneralMultiplicator,3)

summary(HomoMajorArc.AhGh)
round(as.numeric(coefficients(HomoMajorArc.AhGh)[1])*GeneralMultiplicator,3)
round(as.numeric(coefficients(HomoMajorArc.AhGh)[2])*GeneralMultiplicator,3)

###################################################
################ TEST OF PROPORTIONALITY IN INTERCEPT AND SLOPES (DID NOT WORK - BUT WE CAN THINK MORE ON BIGGER DATASETS)
###################################################

###### if increase in intercept and slope is linear? changes in slope / changes in intercept ~ 1
# Homo vs MusOld ChTh: changes in slope / changes in intercept = 1.73987
(as.numeric(coefficients(HomoMajorArc.ChTh)[2])/as.numeric(coefficients(MusOldMajorArc.ChTh)[2])) /
  (as.numeric(coefficients(HomoMajorArc.ChTh)[1])/as.numeric(coefficients(MusOldMajorArc.ChTh)[1])) 

# Homo vs MusOld ChTh: changes in slope / changes in intercept = 0.4557467 => so, so. I expect that it should be >> 1 
# but, probably, initial intercept in young mice was too low... 
(as.numeric(coefficients(HomoMajorArc.AhGh)[2])/as.numeric(coefficients(MusOldMajorArc.AhGh)[2])) /
  (as.numeric(coefficients(HomoMajorArc.AhGh)[1])/as.numeric(coefficients(MusOldMajorArc.AhGh)[1])) 


#################################################################
############# 4: COMPARE SLOPES::::: PAPER!!!!!!!!!!
#################################################################

### real values Homo vs MusOld (AhGh differs stronger than ChTh)
as.numeric(coefficients(HomoMajorArc.ChTh)[2])/as.numeric(coefficients(MusOldMajorArc.ChTh)[2]) # 6.9
as.numeric(coefficients(HomoMajorArc.AhGh)[2])/as.numeric(coefficients(MusOldMajorArc.AhGh)[2]) # 11.7

### bootstrapped values Homo vs MusOld (AhGh differs stronger than ChTh)
summary(VecOfHomoChTh/VecOfOldMiceChTh) # 2.758   5.766   6.810   6.953   8.086  13.276
summary(VecOfHomoAhGh/VecOfOldMiceAhGh) # 3.236   9.824  11.649  11.791  13.491  21.118
wilcox.test(VecOfHomoChTh/VecOfOldMiceChTh,VecOfHomoAhGh/VecOfOldMiceAhGh)
boxplot(VecOfHomoChTh/VecOfOldMiceChTh,VecOfHomoAhGh/VecOfOldMiceAhGh)

### real values Homo vs MusYoung (AhGh differs stronger than ChTh)
as.numeric(coefficients(HomoMajorArc.ChTh)[2])/as.numeric(coefficients(MusYoungMajorArc.ChTh)[2]) # 13.7
as.numeric(coefficients(HomoMajorArc.AhGh)[2])/as.numeric(coefficients(MusYoungMajorArc.AhGh)[2]) # 70.86

### bootstrapped values Homo vs MusYOung (AhGh differs stronger than ChTh)
summary(VecOfHomoChTh/VecOfYoungMiceChTh) # 4.529  11.237  13.630  13.821  16.200  30.459 
summary(VecOfHomoAhGh/VecOfYoungMiceAhGh) # 16.14   57.00   71.43   76.38   88.79  245.52
wilcox.test(VecOfHomoChTh/VecOfYoungMiceChTh,VecOfHomoAhGh/VecOfYoungMiceAhGh)
boxplot(VecOfHomoChTh/VecOfYoungMiceChTh,VecOfHomoAhGh/VecOfYoungMiceAhGh)

### real values MusOld versus MusYoung (AhGh differs stronger than ChTh)
as.numeric(coefficients(MusOldMajorArc.ChTh)[2])/as.numeric(coefficients(MusYoungMajorArc.ChTh)[2]) # 1.98
as.numeric(coefficients(MusOldMajorArc.AhGh)[2])/as.numeric(coefficients(MusYoungMajorArc.AhGh)[2]) # 6.05

### bootstrapped values OldMice vs YoungMice (AhGh differs stronger than ChTh)
summary(VecOfOldMiceChTh/VecOfYoungMiceChTh) # 1.140   1.782   1.983   2.015   2.206   3.355 
summary(VecOfOldMiceAhGh/VecOfYoungMiceAhGh) # 3.045   5.196   6.098   6.518   7.332  18.758
wilcox.test(VecOfOldMiceChTh/VecOfYoungMiceChTh,VecOfOldMiceAhGh/VecOfYoungMiceAhGh)
boxplot(VecOfOldMiceChTh/VecOfYoungMiceChTh,VecOfOldMiceAhGh/VecOfYoungMiceAhGh)

### bootstrapped boxplots:
pdf("../../Body/4Figures/DuplexSeqDataAnalyses03.R.pdf", width = 12, height = 14)
boxplot(VecOfOldMiceChTh/VecOfYoungMiceChTh,VecOfOldMiceAhGh/VecOfYoungMiceAhGh,VecOfHomoChTh/VecOfOldMiceChTh,VecOfHomoAhGh/VecOfOldMiceAhGh,VecOfHomoChTh/VecOfYoungMiceChTh,VecOfHomoAhGh/VecOfYoungMiceAhGh, notch = TRUE, outline = FALSE, col = c('gray','red'), names = c("OldMice/YoungMice","OldMice/YoungMice","Homo/OldMice","Homo/OldMice","Homo/YoungMice","Homo/YoungMice"), ylab = 'ratio of slopes')
dev.off()

###### INTERCEPTS
### Homo vs MusOld (AhGh differs stronger than ChTh)
as.numeric(coefficients(HomoMajorArc.ChTh)[1])/as.numeric(coefficients(MusOldMajorArc.ChTh)[1]) # 3.9768
as.numeric(coefficients(HomoMajorArc.AhGh)[1])/as.numeric(coefficients(MusOldMajorArc.AhGh)[1]) # 25.678

### Homo vs MusYoung (AhGh differs stronger than ChTh)
as.numeric(coefficients(HomoMajorArc.ChTh)[1])/as.numeric(coefficients(MusYoungMajorArc.ChTh)[1]) # 16.33
as.numeric(coefficients(HomoMajorArc.AhGh)[1])/as.numeric(coefficients(MusYoungMajorArc.AhGh)[1]) # 57.92

### MusOld versus MusYoung (AhGh differs stronger than ChTh - NO, it is opposite here or may be non significant)
as.numeric(coefficients(MusOldMajorArc.ChTh)[1])/as.numeric(coefficients(MusYoungMajorArc.ChTh)[1]) # 4.1
as.numeric(coefficients(MusOldMajorArc.AhGh)[1])/as.numeric(coefficients(MusYoungMajorArc.AhGh)[1]) # 2.25


###### BOXPLOTES AND ANALYSES
pdf("../../Body/4Figures/DuplexSeqDataAnalyses02.R.pdf")
par(mfrow=c(2,2))
boxplot(CytbYoungMouseAhGhRelFreqPerBins,CytbOldMouseAhGhRelFreqPerBins,CytbHomoAhGhRelFreqPerBins,names=c('young mice','old mice','human'), notch = TRUE, main = 'Ah>Gh freq in spectrum', col = c(rgb(1,0,0,0.2),rgb(1,0,0,0.6),rgb(1,0,0,0.9)))
boxplot(CytbYoungMouseChThRelFreqPerBins,CytbOldMouseChThRelFreqPerBins,CytbHomoChThRelFreqPerBins,names=c('young mice','old mice','human'), notch = TRUE, main = 'Ch>Th freq in spectrum', col = c(rgb(0,0,0,0.2),rgb(0,0,0,0.4),rgb(0,0,0,0.8)))
#boxplot(CytbYoungMouseAhGhRelFreqPerBins[50:100],CytbOldMouseAhGhRelFreqPerBins[50:100],CytbHomoAhGhRelFreqPerBins[25:50],names=c('young mice','old mice','human'), notch = TRUE)
#boxplot(CytbYoungMouseChThRelFreqPerBins[50:100],CytbOldMouseChThRelFreqPerBins[50:100],CytbHomoChThRelFreqPerBins[25:50],names=c('young mice','old mice','human'), notch = TRUE)
dev.off()

wilcox.test(CytbYoungMouseAhGhRelFreqPerBins,CytbOldMouseAhGhRelFreqPerBins) # 1.583e-08
wilcox.test(CytbOldMouseAhGhRelFreqPerBins,CytbHomoAhGhRelFreqPerBins)       # 2.2e-16
wilcox.test(CytbYoungMouseAhGhRelFreqPerBins,CytbHomoAhGhRelFreqPerBins)     # 2.2e-16  

wilcox.test(CytbYoungMouseChThRelFreqPerBins,CytbOldMouseChThRelFreqPerBins) # 2.2e-16
wilcox.test(CytbOldMouseChThRelFreqPerBins,CytbHomoChThRelFreqPerBins)       # 2.2e-16
wilcox.test(CytbYoungMouseChThRelFreqPerBins,CytbHomoChThRelFreqPerBins)     # 1.805e-05  

######### SCATTERPLOTS AND BINS:

pdf("../../Body/4Figures/DuplexSeqDataAnalyses01.R.pdf", width = 7, height = 14)
par(mfrow=c(1,3))

#### Ch>Th
plot(HomoMajorArc$MajorArcBins*2,HomoMajorArc$ChThFreq*10^5, col = rgb(0,0,0,0.9), xlim = c(1,101), ylim = c(0,8), xlab = 'bins of TBSS', ylab = 'mut freq * 10^-5', pch = 0, main = 'Ch>Th'); par(new = TRUE)
A = HomoMajorArc$ChThFreq*10^5; B = HomoMajorArc$MajorArcBins*2; HomoMajorArc.ChTh = lm(A ~ B)
abline(HomoMajorArc.ChTh, col = rgb(0,0,0,0.9));  par(new = TRUE)

plot(MusOldMajorArc$MajorArcBins,MusOldMajorArc$ChThFreq*10^5, col = rgb(0,0,0,0.6), xlim = c(1,101), ylim = c(0,8), xlab = 'bins of TBSS', ylab = 'mut freq * 10^-5', pch = 16);  par(new = TRUE)
A = MusOldMajorArc$ChThFreq*10^5; B = MusOldMajorArc$MajorArcBins; MusOldMajorArc.ChTh = lm(A ~ B)
abline(MusOldMajorArc.ChTh, col = rgb(0,0,0,0.6));  par(new = TRUE)

plot(MusYoungMajorArc$MajorArcBins,MusYoungMajorArc$ChThFreq*10^5, col = rgb(0,0,0,0.2), xlim = c(1,101), ylim = c(0,8), xlab = 'bins of TBSS', ylab = 'mut freq * 10^-5', pch = 6)
A = MusYoungMajorArc$ChThFreq*10^5; B = MusYoungMajorArc$MajorArcBins; MusYoungMajorArc.ChTh = lm(A ~ B)
abline(MusYoungMajorArc.ChTh, col = rgb(0,0,0,0.2))

#### Ah>Gh
plot(HomoMajorArc$MajorArcBins*2,HomoMajorArc$AhGhFreq*10^5, col = rgb(1,0,0,0.9), xlim = c(1,101), ylim = c(0,8), xlab = 'bins of TBSS', ylab = 'mut freq * 10^-5', pch = 0, main = 'Ah>Gh'); par(new = TRUE)
A = HomoMajorArc$AhGhFreq*10^5; B = HomoMajorArc$MajorArcBins*2; HomoMajorArc.AhGh = lm(A ~ B)
abline(HomoMajorArc.AhGh, col = rgb(1,0,0,0.9));  par(new = TRUE)

plot(MusOldMajorArc$MajorArcBins,MusOldMajorArc$AhGhFreq*10^5, col = rgb(1,0,0,0.6), xlim = c(1,101), ylim = c(0,8), xlab = 'bins of TBSS', ylab = 'mut freq * 10^-5', pch = 1);  par(new = TRUE)
A = MusOldMajorArc$AhGhFreq*10^5; B = MusOldMajorArc$MajorArcBins; MusOldMajorArc.AhGh = lm(A ~ B)
abline(MusOldMajorArc.AhGh, col = rgb(1,0,0,0.6));  par(new = TRUE)
  
plot(MusYoungMajorArc$MajorArcBins,MusYoungMajorArc$AhGhFreq*10^5, col = rgb(1,0,0,0.2), xlim = c(1,101), ylim = c(0,8), xlab = 'bins of TBSS', ylab = 'mut freq * 10^-5', pch = 6)
A = MusYoungMajorArc$AhGhFreq*10^5; B = MusYoungMajorArc$MajorArcBins; MusYoungMajorArc.AhGh = lm(A ~ B)
abline(MusYoungMajorArc.AhGh, col = rgb(1,0,0,0.2))

barplot(
c(
  as.numeric(coefficients(MusYoungMajorArc.ChTh)[2]),
  as.numeric(coefficients(MusYoungMajorArc.AhGh)[2]),
  as.numeric(coefficients(MusOldMajorArc.ChTh)[2]),
  as.numeric(coefficients(MusOldMajorArc.AhGh)[2]),
  as.numeric(coefficients(HomoMajorArc.ChTh)[2]),
  as.numeric(coefficients(HomoMajorArc.AhGh)[2])
),
col=c(rgb(0,0,0,0.2),rgb(1,0,0,0.2),rgb(0,0,0,0.4),rgb(1,0,0,0.6),rgb(0,0,0,0.7),rgb(1,0,0,0.9)),
ylab = 'slopes'
)

dev.off()


##### all these numbers just mention in text
### ratio of AhGh / ChTh is increasing with age
c(
as.numeric(coefficients(MusYoungMajorArc.AhGh)[2])/as.numeric(coefficients(MusYoungMajorArc.ChTh)[2]), # 0.058
as.numeric(coefficients(MusOldMajorArc.AhGh)[2])/as.numeric(coefficients(MusOldMajorArc.ChTh)[2]),     # 0.18
as.numeric(coefficients(HomoMajorArc.AhGh)[2])/as.numeric(coefficients(HomoMajorArc.ChTh)[2])          # 0.30
)

### relative increase if we use young mouse as baseline is stronger for AhGh

c(
  as.numeric(coefficients(MusOldMajorArc.ChTh)[2])/as.numeric(coefficients(MusYoungMajorArc.ChTh)[2]),  # 1.98
  as.numeric(coefficients(MusOldMajorArc.AhGh)[2])/as.numeric(coefficients(MusYoungMajorArc.AhGh)[2]),  # 6.05
  
  as.numeric(coefficients(HomoMajorArc.ChTh)[2])/as.numeric(coefficients(MusYoungMajorArc.ChTh)[2]),  # 6.85
  as.numeric(coefficients(HomoMajorArc.AhGh)[2])/as.numeric(coefficients(MusYoungMajorArc.AhGh)[2])  # 35
)

#### BIND ALL DATA INTO ONE BIG DATASET

# MusYoung 
MusYoungMajorArc.1 = MusYoungMajorArc[colnames(MusYoungMajorArc) %in% c('AhGhFreq','MajorArcBins')]
names(MusYoungMajorArc.1)=c('Freq','MajorArcBins')
MusYoungMajorArc.1$Mut = 1 # means AhGh
MusYoungMajorArc.1$Age = 0 # means mus young

MusYoungMajorArc.2 = MusYoungMajorArc[colnames(MusYoungMajorArc) %in% c('ChThFreq','MajorArcBins')]
names(MusYoungMajorArc.2)=c('Freq','MajorArcBins')
MusYoungMajorArc.2$Mut = 0 # means C>T
MusYoungMajorArc.2$Age = 0 # means mus young

MusYoungMajorArc = rbind(MusYoungMajorArc.1,MusYoungMajorArc.2)

# MusOld
MusOldMajorArc.1 = MusOldMajorArc[colnames(MusOldMajorArc) %in% c('AhGhFreq','MajorArcBins')]
names(MusOldMajorArc.1)=c('Freq','MajorArcBins')
MusOldMajorArc.1$Mut = 1 # means A>G 
MusOldMajorArc.1$Age = 1 # means mus old

MusOldMajorArc.2 = MusOldMajorArc[colnames(MusOldMajorArc) %in% c('ChThFreq','MajorArcBins')]
names(MusOldMajorArc.2)=c('Freq','MajorArcBins')
MusOldMajorArc.2$Mut = 0 # means C>T
MusOldMajorArc.2$Age = 1 # means mus old

MusOldMajorArc = rbind(MusOldMajorArc.1,MusOldMajorArc.2)

# Homo
HomoMajorArc.1 = HomoMajorArc[colnames(HomoMajorArc) %in% c('AhGhFreq','MajorArcBins')]
names(HomoMajorArc.1)=c('Freq','MajorArcBins')
HomoMajorArc.1$Mut = 1 # means AhGh
HomoMajorArc.1$Age = 2 # means homo

HomoMajorArc.2 = HomoMajorArc[colnames(HomoMajorArc) %in% c('ChThFreq','MajorArcBins')]
names(HomoMajorArc.2)=c('Freq','MajorArcBins')
HomoMajorArc.2$Mut = 0 # means C>T
HomoMajorArc.2$Age = 2 # means mus young

HomoMajorArc = rbind(HomoMajorArc.1,HomoMajorArc.2)


### bind all three datasets:
AllTogether = rbind(MusYoungMajorArc,MusOldMajorArc,HomoMajorArc)
table(AllTogether$Age)
table(AllTogether$Mut)

### only A>G in mice: old mice have a bit higher intercept and significantly higher slope
AG = AllTogether[AllTogether$Mut == 1 & AllTogether$Age < 2,]
summary(lm(AG$Freq ~ AG$MajorArcBins + AG$Age))
summary(lm(AG$Freq ~ AG$MajorArcBins*AG$Age)) ###
#  (Intercept)            1.736e-07  5.027e-08   3.453 0.000678 ***
#  AG$MajorArcBins        2.092e-09  8.557e-10   2.445 0.015365 *  
#  AG$Age                 2.180e-07  7.109e-08   3.066 0.002473 **  # old mice have higher intercept
#  AG$MajorArcBins:AG$Age 1.058e-08  1.210e-09   8.739 9.98e-16 *** # old mice have higher slope:

### only C>T in mice. old mice have a higher intercept and higher slope
CT = AllTogether[AllTogether$Mut == 0 & AllTogether$Age < 2,]
summary(lm(CT$Freq ~ CT$MajorArcBins + CT$Age))
summary(lm(CT$Freq ~ CT$MajorArcBins*CT$Age))
#  (Intercept)            2.027e-06  3.671e-07   5.523 1.04e-07 ***
#  CT$MajorArcBins        3.552e-08  6.248e-09   5.684 4.66e-08 ***
#  CT$Age                 6.299e-06  5.191e-07  12.135  < 2e-16 *** # old mice have higher intercept
#  CT$MajorArcBins:CT$Age 3.481e-08  8.836e-09   3.939 0.000113 *** # old nice have higher slope 3.481e-08/3.552e-08


### HomoVsYoungMiceAG
HomoVsYoungMiceAG = AllTogether[AllTogether$Mut == 1 & AllTogether$Age != 1,]
HomoVsYoungMiceAG[HomoVsYoungMiceAG$Age == 2,]$Age = 1; table(HomoVsYoungMiceAG$Age) # recoded Age
summary(lm(HomoVsYoungMiceAG$Freq ~ HomoVsYoungMiceAG$MajorArcBins + HomoVsYoungMiceAG$Age))
summary(lm(HomoVsYoungMiceAG$Freq ~ HomoVsYoungMiceAG$MajorArcBins*HomoVsYoungMiceAG$Age))

HomoVsYoungMiceCT = AllTogether[AllTogether$Mut == 0 & AllTogether$Age != 1,]
HomoVsYoungMiceCT[HomoVsYoungMiceCT$Age == 2,]$Age = 1; table(HomoVsYoungMiceCT$Age) # recoded Age
summary(lm(HomoVsYoungMiceCT$Freq ~ HomoVsYoungMiceCT$MajorArcBins + HomoVsYoungMiceCT$Age))
summary(lm(HomoVsYoungMiceCT$Freq ~ HomoVsYoungMiceCT$MajorArcBins*HomoVsYoungMiceCT$Age))


# Freq is increasing with bins, it is higher in C>T and old animals
summary(lm(AllTogether$Freq ~ AllTogether$MajorArcBins + AllTogether$Mut +  AllTogether$Age))

# Freq is increasing with bins, it is higher in C>T and old animals and also, the slope is less in A>G 
summary(lm(AllTogether$Freq ~ AllTogether$MajorArcBins*AllTogether$Mut +  AllTogether$Age))

# Freq is increasing with bins, it is higher in C>T and old animals and also, the slope is less in A>G 
summary(lm(AllTogether$Freq ~ AllTogether$MajorArcBins*AllTogether$Mut + AllTogether$Mut*AllTogether$Age))




A = lm(AllTogether$Freq ~ AllTogether$MajorArcBins*AllTogether$Mut)


############## LM TEST:

x <- runif(100, 1, 100)
A = 2
b = 4
y = A*x+b
LmRes = lm(y~x)
LmRes

y = (A*x+b)*2
LmRes = lm(y~x)
LmRes




  
  
  