rm(list=ls(all=TRUE))

AASubs = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecAminoAcidSubs.txt', header = TRUE)
AANumbers = read.table('../../Body/3Results/AminoAcidFreqsChordata.txt', header = TRUE, sep=",") 

AASubsCB =  AASubs[AASubs$Gene == "CytB",]
table(AASubsCB$Gene)
AANumbersCB = AANumbers[AANumbers$Gene == "CytB",]
table(AANumbersCB$Gene)

MergedAAData = merge(AASubsCB, AANumbersCB, by='Species')
length(unique(MergedAAData$Species))

Temp = data.frame(t(data.frame(strsplit(as.character(MergedAAData$TypesOfAASub), "_"))), row.names = NULL)
MergedAAData = cbind(MergedAAData, Temp)
head(MergedAAData)
tail(MergedAAData)

table(MergedAAData$X1)
table(MergedAAData$X2)
MergedAAData$X1 = as.character(MergedAAData$X1)
MergedAAData$X2= as.character(MergedAAData$X2)
str(MergedAAData$X1)
str(MergedAAData$X2)

MergedAAData[MergedAAData$X1 == "LeuCT",]$X1 = "Leu"
MergedAAData[MergedAAData$X1 == "LeuTT",]$X1 = "Leu"
MergedAAData[MergedAAData$X1 == "SerAG",]$X1 = "Ser"
MergedAAData[MergedAAData$X1 == "SerTC",]$X1 = "Ser"

MergedAAData[MergedAAData$X2 == "LeuCT",]$X2 = "Leu"
MergedAAData[MergedAAData$X2 == "LeuTT",]$X2 = "Leu"
MergedAAData[MergedAAData$X2 == "SerAG",]$X2 = "Ser"
MergedAAData[MergedAAData$X2 == "SerTC",]$X2 = "Ser"

table(MergedAAData$X1)
table(MergedAAData$X2)

head(MergedAAData)
names(MergedAAData)[c(length(MergedAAData)-1, length(MergedAAData))] <- c("AncestralAA", "DescendantAA")
names(MergedAAData)

M = MergedAAData
for (i in 1:nrow(M)){
  NeedColName = names(M)[M$AncestralAA[i] == names(M)]
  if(M$AncestralAA[i] == NeedColName ){
    M$NumberAncAA[i] = M[NeedColName][i,]
  }
  NeedColName2 = names(M)[M$DescendantAA[i] == names(M)]
  if(M$DescendantAA[i] == NeedColName2 ){
    M$NumberAncAABackward[i] = M[NeedColName2][i,]
  }
}

tail(M, 10)

M$TypesOfAASubBackward = M$TypesOfAASub
Temp = data.frame(t(data.frame(strsplit(as.character(MergedAAData$TypesOfAASub), "_"))), row.names = NULL)
Temp$X3 = paste(Temp$X2, "_", Temp$X1, sep="")
M$TypesOfAASubBackward = Temp$X3
M$FreqOfSubBackward

Final = c()
ListOfSpecies = as.character(unique(M$Species))
#y="Abbottina_rivularis" 
#z = "LeuCT_Phe"
for (y in ListOfSpecies){
  TempFrame = M[M$Species == y,]
  for (i in 1:nrow(TempFrame)){
    for (z in as.character(TempFrame$TypesOfAASub)){
    if(TempFrame$TypesOfAASubBackward[i] == z){
      TempFrame$FreqOfSubBackward[i] = as.character(TempFrame[TempFrame$TypesOfAASub == z,]$FreqOfSub) 
      } else {
      TempFrame$FreqOfSubBackward[i] = 0
      }
    }
  }
  Final = rbind(Final, TempFrame)
}

table(Final$FreqOfSubBackward)

names(Final)

FinalToWrite=Final[,c("Species", "TypesOfAASub", "FreqOfSub", "NumberAncAA", "TypesOfAASubBackward", "FreqOfSubBackward", "NumberAncAABackward")]
write.table(FinalToWrite, file = '../../Body/3Results/VertebratePolymorphisms.MutSpecAminoAcidSubsWithRefSeqData.txt', quote = FALSE)

