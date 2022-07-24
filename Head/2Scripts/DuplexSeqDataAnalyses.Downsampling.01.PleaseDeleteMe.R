##############################################################################################################
################################################### DOWNSAMPLIG OF C>T TOWARDS A>G IN YOUNG MICE
#############################################################################################################

## IF I DO SAMPLING FROM FREQUENCIES DIRECTLY? IN THIS CASE I DON'T NEED TO COME BACK TO NORMALIZATIONS AND DO EVERYTHING SIMPLER?

rm(list=ls(all=TRUE))

pdf("../../Body/4Figures/DuplexSeqDataAnalyses.Downsampling.01.R.pdf", width = 12, height = 14)
par(mfrow=c(3,1))

##### YOUNG MICE 
MusYoung = read.table("../../Body/1Raw/DuplexSeqData/DataS3.txt", sep = '\t', header = TRUE)
names(MusYoung)
MusYoungMajorArc=MusYoung[MusYoung$BIN_NUM >= 53 & MusYoung$BIN_NUM <= 153,]
MusYoungMajorArc$MajorArcBins = seq(1:101)

Sampling = MusYoungMajorArc
## start of the utilites (identical for different datasets)
TotalAhGh = sum(Sampling$T.C_COUNT); TotalAhGh # 369
TotalChTh = sum(Sampling$G.A_COUNT); TotalChTh # 2068

##### weather there is a correlation between ChTh and AhGh denominators and bins?
cor.test(Sampling$G.A_DENOM,Sampling$MajorArcBins, method = 'spearman') # negative trend... Ch>Th we divide by more and more small numbers during normalization - this increase frequencies of late and increases slope
cor.test(Sampling$T.C_DENOM,Sampling$MajorArcBins, method = 'spearman') # negative weak trend...
t.test(Sampling$G.A_DENOM,Sampling$T.C_DENOM) # mean of x = 5355255 mean of y = 13179759; 13179759/5355255 = 2.46 bigger DENOM in Ah>Gh
### decoding into individual mutations, sampling and decoding back
TotalChThVectorOfBins = c()
for (i in 1:nrow(Sampling))
{ # i = 1
  NumberOfChThPerBin = Sampling$G.A_COUNT[i]
  BinNumber = Sampling$MajorArcBins[i]
  TotalChThVectorOfBins = c(TotalChThVectorOfBins,rep(BinNumber,NumberOfChThPerBin))
}
length(TotalChThVectorOfBins)   # == TotalChTh, good

VecOfSlopesChTh=c()
for (inter in 1:1000)
{
  DownSampledChThVectorOfBins = sample(TotalChThVectorOfBins,TotalAhGh)
  length(DownSampledChThVectorOfBins)
  DownSampledChTh = data.frame(table(DownSampledChThVectorOfBins))
  names(DownSampledChTh) = c('MajorArcBins','DownSampledChThCounts')
  
  Sampling = Sampling[,1:51]  # delete every time extra column with name DownSampledChThCounts
  dim(Sampling)
  Sampling = merge(Sampling,DownSampledChTh, all.x = TRUE)
  Sampling$DownSampledChThCounts[is.na(Sampling$DownSampledChThCounts)] <- 0
  dim(Sampling)
  Sampling$DownSampledChThFreq = Sampling$DownSampledChThCounts / Sampling$G.A_DENOM
#  Sampling$DownSampledChThFreq = Sampling$DownSampledChThCounts # comment it 
   VecOfSlopesChTh = c(VecOfSlopesChTh,summary(lm(Sampling$DownSampledChThFreq ~ Sampling$MajorArcBins))$coefficients[2])
#  VecOfSlopesChTh = c(VecOfSlopesChTh,summary(lm(Sampling$DownSampledChThFreq ~ 0 + Sampling$MajorArcBins))$coefficients[1])
#  VecOfSlopesChTh = c(VecOfSlopesChTh,summary(lm(scale(Sampling$DownSampledChThFreq) ~ 0 + scale(Sampling$MajorArcBins)))$coefficients[1])
}

Sampling$AhGhFreq = Sampling$T.C_COUNT / Sampling$T.C_DENOM
#Sampling$AhGhFreq = Sampling$T.C_COUNT  # comment it 

## end of the utilites:

  YoungMiceSlopeAhGh = summary(lm(Sampling$AhGhFreq ~ Sampling$MajorArcBins))$coefficients[2]
# YoungMiceSlopeAhGh = summary(lm(Sampling$AhGhFreq ~ 0 + Sampling$MajorArcBins))$coefficients[1]
# YoungMiceSlopeAhGh = summary(lm(scale(Sampling$AhGhFreq) ~ 0 +  scale(Sampling$MajorArcBins)))$coefficients[1]

YoungMiceVecOfSlopesChTh = VecOfSlopesChTh
hist(YoungMiceVecOfSlopesChTh,xlim = c(0,4.282667e-07))
abline(v=YoungMiceSlopeAhGh, col = 'red', lwd = 4)
abline(v=median(YoungMiceVecOfSlopesChTh), col = 'dark grey', lwd = 4)

###### OLD MICE 
MusOld = read.table("../../Body/1Raw/DuplexSeqData/DataS2.txt", sep = '\t', header = TRUE)
names(MusOld)
MusOldMajorArc=MusOld[MusOld$BIN_NUM >= 53 & MusOld$BIN_NUM <= 153,]
MusOldMajorArc$MajorArcBins = seq(1:101)

Sampling = MusOldMajorArc
## start of the utilites (identical for different datasets)
TotalAhGh = sum(Sampling$T.C_COUNT); TotalAhGh # 1440
TotalChTh = sum(Sampling$G.A_COUNT); TotalChTh # 6599

##### weather there is a correlation between ChTh and AhGh denominators and bins?
cor.test(Sampling$G.A_DENOM,Sampling$MajorArcBins, method = 'spearman') # negative trend... Ch>Th we divide by more and more small numbers during normalization - this increase frequencies of late and increases slope
cor.test(Sampling$T.C_DENOM,Sampling$MajorArcBins, method = 'spearman') # negative weak trend...
t.test(Sampling$G.A_DENOM,Sampling$T.C_DENOM) # mean of x = 5637479 mean of y = 13965400; 13965400/5637479 = 2.47 bigger DENOM in Ah>Gh

### decoding into individual mutations, sampling and decoding back
TotalChThVectorOfBins = c()
for (i in 1:nrow(Sampling))
{ # i = 1
  NumberOfChThPerBin = Sampling$G.A_COUNT[i]
  BinNumber = Sampling$MajorArcBins[i]
  TotalChThVectorOfBins = c(TotalChThVectorOfBins,rep(BinNumber,NumberOfChThPerBin))
}
length(TotalChThVectorOfBins)   # == TotalChTh, good

VecOfSlopesChTh=c()
for (inter in 1:1000)
{
  DownSampledChThVectorOfBins = sample(TotalChThVectorOfBins,TotalAhGh)
  length(DownSampledChThVectorOfBins)
  DownSampledChTh = data.frame(table(DownSampledChThVectorOfBins))
  names(DownSampledChTh) = c('MajorArcBins','DownSampledChThCounts')
  
  Sampling = Sampling[,1:51]  # delete every time extra column with name DownSampledChThCounts
  dim(Sampling)
  Sampling = merge(Sampling,DownSampledChTh, all.x = TRUE)
  Sampling$DownSampledChThCounts[is.na(Sampling$DownSampledChThCounts)] <- 0
  dim(Sampling)
  Sampling$DownSampledChThFreq = Sampling$DownSampledChThCounts / Sampling$G.A_DENOM
#  Sampling$DownSampledChThFreq = Sampling$DownSampledChThCounts # comment it 
   VecOfSlopesChTh = c(VecOfSlopesChTh,summary(lm(Sampling$DownSampledChThFreq ~ Sampling$MajorArcBins))$coefficients[2])
#  VecOfSlopesChTh = c(VecOfSlopesChTh,summary(lm(Sampling$DownSampledChThFreq ~ 0 + Sampling$MajorArcBins))$coefficients[1])
#  VecOfSlopesChTh = c(VecOfSlopesChTh,summary(lm(scale(Sampling$DownSampledChThFreq) ~ 0 + scale(Sampling$MajorArcBins)))$coefficients[1])
}

Sampling$AhGhFreq = Sampling$T.C_COUNT / Sampling$T.C_DENOM
#Sampling$AhGhFreq = Sampling$T.C_COUNT  # comment it 

## end of the utilites:

OldMiceSlopeAhGh = summary(lm(Sampling$AhGhFreq ~ Sampling$MajorArcBins))$coefficients[2]
# OldMiceSlopeAhGh = summary(lm(Sampling$AhGhFreq ~ 0 + Sampling$MajorArcBins))$coefficients[1]
# OldMiceSlopeAhGh = summary(lm(scale(Sampling$AhGhFreq) ~ 0 +  scale(Sampling$MajorArcBins)))$coefficients[1]
OldMiceVecOfSlopesChTh = VecOfSlopesChTh
hist(OldMiceVecOfSlopesChTh,xlim = c(0,4.282667e-07))
abline(v=OldMiceSlopeAhGh, col = 'red', lwd = 4)
abline(v=median(OldMiceVecOfSlopesChTh), col = 'dark grey', lwd = 4)


##### HUMAN 
Homo = read.table("../../Body/1Raw/DuplexSeqData/DataS6.txt", sep = '\t', header = TRUE)
names(Homo)
HomoMajorArc=Homo[Homo$BIN_NUM >= 28 & Homo$BIN_NUM <= 76,] # 53/2 == 26.5 > start from 28 153/2 = 76.5 > end at 76
HomoMajorArc$MajorArcBins = seq(1:nrow(HomoMajorArc))

Sampling = HomoMajorArc
## start of the utilites (identical for different datasets)
TotalAhGh = sum(Sampling$T.C_COUNT); TotalAhGh # 1553
TotalChTh = sum(Sampling$G.A_COUNT); TotalChTh # 2442

##### weather there is a correlation between ChTh and AhGh denominators and bins?
cor.test(Sampling$G.A_DENOM,Sampling$MajorArcBins, method = 'spearman') # negative  trend... Ch>Th we divide by more and more small numbers during normalization - this increase frequencies of late and increases slope
cor.test(Sampling$T.C_DENOM,Sampling$MajorArcBins, method = 'spearman') # nothing
t.test(Sampling$G.A_DENOM,Sampling$T.C_DENOM) # mean of x = 1095040   mean of y = 2301627; 2301627/1095040 = 2.1 bigger DENOM in Ah>Gh 

### decoding into individual mutations, sampling and decoding back
TotalChThVectorOfBins = c()
for (i in 1:nrow(Sampling))
{ # i = 1
  NumberOfChThPerBin = Sampling$G.A_COUNT[i]
  BinNumber = Sampling$MajorArcBins[i]
  TotalChThVectorOfBins = c(TotalChThVectorOfBins,rep(BinNumber,NumberOfChThPerBin))
}
length(TotalChThVectorOfBins)   # == TotalChTh, good

VecOfSlopesChTh=c()
for (inter in 1:1000)
{
  DownSampledChThVectorOfBins = sample(TotalChThVectorOfBins,TotalAhGh)
  length(DownSampledChThVectorOfBins)
  DownSampledChTh = data.frame(table(DownSampledChThVectorOfBins))
  names(DownSampledChTh) = c('MajorArcBins','DownSampledChThCounts')
  
  Sampling = Sampling[,1:50]  # delete every time extra column with name DownSampledChThCounts
  dim(Sampling)
  Sampling = merge(Sampling,DownSampledChTh, all.x = TRUE)
  Sampling$DownSampledChThCounts[is.na(Sampling$DownSampledChThCounts)] <- 0
  dim(Sampling)
  Sampling$DownSampledChThFreq = Sampling$DownSampledChThCounts / Sampling$G.A_DENOM
# Sampling$DownSampledChThFreq = Sampling$DownSampledChThCounts # comment it 
  VecOfSlopesChTh = c(VecOfSlopesChTh,summary(lm(Sampling$DownSampledChThFreq ~ Sampling$MajorArcBins))$coefficients[2])
# VecOfSlopesChTh = c(VecOfSlopesChTh,summary(lm(Sampling$DownSampledChThFreq ~ 0 + Sampling$MajorArcBins))$coefficients[1])
# VecOfSlopesChTh = c(VecOfSlopesChTh,summary(lm(scale(Sampling$DownSampledChThFreq) ~ 0 + scale(Sampling$MajorArcBins)))$coefficients[1])
}

Sampling$AhGhFreq = Sampling$T.C_COUNT / Sampling$T.C_DENOM
#Sampling$AhGhFreq = Sampling$T.C_COUNT  # comment it 

## end of the utilites:

  HumanSlopeAhGh = summary(lm(Sampling$AhGhFreq ~ Sampling$MajorArcBins))$coefficients[2]
# HumanSlopeAhGh = summary(lm(Sampling$AhGhFreq ~ 0 + Sampling$MajorArcBins))$coefficients[1]
# HumanSlopeAhGh = summary(lm(scale(Sampling$AhGhFreq) ~ 0 +  scale(Sampling$MajorArcBins)))$coefficients[1]
HumanVecOfSlopesChTh = VecOfSlopesChTh
summary(HumanVecOfSlopesChTh)
hist(HumanVecOfSlopesChTh,xlim = c(0,4.282667e-07)) 
abline(v=HumanSlopeAhGh, col = 'red', lwd = 4)
abline(v=median(HumanVecOfSlopesChTh), col = 'dark grey', lwd = 4)

##### COMPARISONS OF SLOPES:

HumanSlopeAhGh/YoungMiceSlopeAhGh   # 70
HumanSlopeAhGh/OldMiceSlopeAhGh     # 11 => too low, why? how was before?
OldMiceSlopeAhGh/YoungMiceSlopeAhGh # 6

mean(HumanVecOfSlopesChTh)/mean(YoungMiceVecOfSlopesChTh)   # 49
mean(HumanVecOfSlopesChTh)/mean(OldMiceVecOfSlopesChTh)     # 20
mean(OldMiceVecOfSlopesChTh)/mean(YoungMiceVecOfSlopesChTh) # 2.43

dev.off()
