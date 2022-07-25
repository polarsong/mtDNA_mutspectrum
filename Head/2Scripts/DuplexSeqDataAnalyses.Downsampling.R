##############################################################################################################
################################################### DOWNSAMPLIG OF C>T TOWARDS A>G IN YOUNG MICE
#############################################################################################################

rm(list=ls(all=TRUE))

################## START OF SAMPLING ##################
####  we assume that in each bin for each type of substitution there are black (COUNT) and white (DENOM) balls.
####  total number of balls fluctuate between types of mutations and bins => we want to sample it making similar and rerun key analyses
####################################################### 

########## YOUNG MICE 
MusYoung = read.table("../../Body/1Raw/DuplexSeqData/DataS3.txt", sep = '\t', header = TRUE)
names(MusYoung)
MusYoungMajorArc=MusYoung[MusYoung$BIN_NUM >= 53 & MusYoung$BIN_NUM <= 153,]
MusYoungMajorArc$MajorArcBins = seq(1:101)

MusYoungMajorArc = MusYoungMajorArc[names(MusYoungMajorArc) %in% c('G.A_COUNT','T.C_COUNT','G.A_DENOM','T.C_DENOM','MajorArcBins')]
min(MusYoungMajorArc$G.A_DENOM)
min(MusYoungMajorArc$T.C_DENOM)

MusYoungAhGhDownSampledSlopes = c()
MusYoungChThDownSampledSlopes = c()
for (iteration in 1:500)
{
  NumberOfRandomlySampledMolecules = 400000 # this is less than minimum, however we can change it
  VectorChThDownSampledFreq = c()
  VectorAhGhDownSampledFreq = c()
  for (i in 1:nrow(MusYoungMajorArc))
    {
    AhGhTotal = c(rep(0,MusYoungMajorArc$T.C_DENOM[i]),rep(1,MusYoungMajorArc$T.C_COUNT[i])); 
    AhGhSampling = sample(AhGhTotal,NumberOfRandomlySampledMolecules)
    VectorAhGhDownSampledFreq = c(VectorAhGhDownSampledFreq,mean(AhGhSampling))

    ChThTotal = c(rep(0,MusYoungMajorArc$G.A_DENOM[i]),rep(1,MusYoungMajorArc$G.A_COUNT[i])); 
    ChThSampling = sample(ChThTotal,NumberOfRandomlySampledMolecules)
    VectorChThDownSampledFreq = c(VectorChThDownSampledFreq,mean(ChThSampling))
    }
  MusYoungAhGhDownSampledSlopes = c(MusYoungAhGhDownSampledSlopes,summary(lm(VectorAhGhDownSampledFreq ~ MusYoungMajorArc$MajorArcBins))$coefficients[2])
  MusYoungChThDownSampledSlopes = c(MusYoungChThDownSampledSlopes,summary(lm(VectorChThDownSampledFreq ~ MusYoungMajorArc$MajorArcBins))$coefficients[2])
}
summary(MusYoungAhGhDownSampledSlopes)
summary(MusYoungChThDownSampledSlopes)


########## OLD MICE 
MusOld = read.table("../../Body/1Raw/DuplexSeqData/DataS2.txt", sep = '\t', header = TRUE)
names(MusOld)
MusOldMajorArc=MusOld[MusOld$BIN_NUM >= 53 & MusOld$BIN_NUM <= 153,]
MusOldMajorArc$MajorArcBins = seq(1:101)

MusOldMajorArc = MusOldMajorArc[names(MusOldMajorArc) %in% c('G.A_COUNT','T.C_COUNT','G.A_DENOM','T.C_DENOM','MajorArcBins')]
min(MusOldMajorArc$G.A_DENOM) 
min(MusOldMajorArc$T.C_DENOM)

MusOldAhGhDownSampledSlopes = c()
MusOldChThDownSampledSlopes = c()
for (iteration in 1:500)
{
  NumberOfRandomlySampledMolecules = 400000 # this is less than minimum, however we can change it
  VectorChThDownSampledFreq = c()
  VectorAhGhDownSampledFreq = c()
  for (i in 1:nrow(MusOldMajorArc))
  {
    AhGhTotal = c(rep(0,MusOldMajorArc$T.C_DENOM[i]),rep(1,MusOldMajorArc$T.C_COUNT[i])); 
    AhGhSampling = sample(AhGhTotal,NumberOfRandomlySampledMolecules)
    VectorAhGhDownSampledFreq = c(VectorAhGhDownSampledFreq,mean(AhGhSampling))
    
    ChThTotal = c(rep(0,MusOldMajorArc$G.A_DENOM[i]),rep(1,MusOldMajorArc$G.A_COUNT[i])); 
    ChThSampling = sample(ChThTotal,NumberOfRandomlySampledMolecules)
    VectorChThDownSampledFreq = c(VectorChThDownSampledFreq,mean(ChThSampling))
  }
  MusOldAhGhDownSampledSlopes = c(MusOldAhGhDownSampledSlopes,summary(lm(VectorAhGhDownSampledFreq ~ MusOldMajorArc$MajorArcBins))$coefficients[2])
  MusOldChThDownSampledSlopes = c(MusOldChThDownSampledSlopes,summary(lm(VectorChThDownSampledFreq ~ MusOldMajorArc$MajorArcBins))$coefficients[2])
}

summary(MusOldAhGhDownSampledSlopes)
summary(MusOldChThDownSampledSlopes)


########## HUMAN
Human = read.table("../../Body/1Raw/DuplexSeqData/DataS6.txt", sep = '\t', header = TRUE)
names(Human)
HumanMajorArc=Human[Human$BIN_NUM >= 28 & Human$BIN_NUM <= 76,] # 53/2 == 26.5 > start from 28 153/2 = 76.5 > end at 76
HumanMajorArc$MajorArcBins = seq(1:49)

HumanMajorArc = HumanMajorArc[names(HumanMajorArc) %in% c('G.A_COUNT','T.C_COUNT','G.A_DENOM','T.C_DENOM','MajorArcBins')]
min(HumanMajorArc$G.A_DENOM) # 417934
min(HumanMajorArc$T.C_DENOM) # 1536386

HumanAhGhDownSampledSlopes = c()
HumanChThDownSampledSlopes = c()
for (iteration in 1:500)
{
  NumberOfRandomlySampledMolecules = 400000 # this is less than minimum, however we can change it
  VectorChThDownSampledFreq = c()
  VectorAhGhDownSampledFreq = c()
  for (i in 1:nrow(HumanMajorArc))
  {
    AhGhTotal = c(rep(0,HumanMajorArc$T.C_DENOM[i]),rep(1,HumanMajorArc$T.C_COUNT[i])); 
    AhGhSampling = sample(AhGhTotal,NumberOfRandomlySampledMolecules)
    VectorAhGhDownSampledFreq = c(VectorAhGhDownSampledFreq,mean(AhGhSampling))
    
    ChThTotal = c(rep(0,HumanMajorArc$G.A_DENOM[i]),rep(1,HumanMajorArc$G.A_COUNT[i])); 
    ChThSampling = sample(ChThTotal,NumberOfRandomlySampledMolecules)
    VectorChThDownSampledFreq = c(VectorChThDownSampledFreq,mean(ChThSampling))
  }
  HumanAhGhDownSampledSlopes = c(HumanAhGhDownSampledSlopes,summary(lm(VectorAhGhDownSampledFreq ~ HumanMajorArc$MajorArcBins))$coefficients[2])
  HumanChThDownSampledSlopes = c(HumanChThDownSampledSlopes,summary(lm(VectorChThDownSampledFreq ~ HumanMajorArc$MajorArcBins))$coefficients[2])
}

summary(HumanAhGhDownSampledSlopes)
summary(HumanChThDownSampledSlopes)

##### Save all slopes into a data-frame and write it into the file.

AllSlopes = data.frame(MusYoungAhGhDownSampledSlopes, MusYoungChThDownSampledSlopes,MusOldAhGhDownSampledSlopes, MusOldChThDownSampledSlopes, HumanAhGhDownSampledSlopes,HumanChThDownSampledSlopes)
write.table(AllSlopes,"../../Body/3Results/DuplexSeqDataAnalyses.Downsampling.R.txt")

