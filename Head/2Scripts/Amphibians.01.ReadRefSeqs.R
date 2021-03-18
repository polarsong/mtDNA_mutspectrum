rm(list=ls(all=TRUE))

library(seqinr)

#### read fasta file

# MutSpec = read.fasta("../../Body/1Raw/AmphibiansMtDnaRefSeqs/sequence.fasta")  # Amphibians
# MutSpec = read.fasta("../../Body/1Raw/LepidosauriaMtDnaRefSeqs/sequence.fasta")  # Lepidosauria
MutSpec = read.fasta("../../Body/1Raw/PrimatesMtDnaRefSeqs/sequence.fasta")  # Primates

Seq = getSequence(MutSpec)
Names = getName(MutSpec)

#### make simple data-frame with RefSeqId and WholeGenome
DataFrame = c()
for (i in 1:length(Names))
{ # i = 1
seq=paste(unlist(Seq[i]), collapse = '')
name = Names[i]
DataFrame = rbind(DataFrame,c(name,seq))
}
DataFrame = data.frame(DataFrame)
names(DataFrame)=c('RefSeqId','WholeGenome')

#### looking for di-nucleotide patterns
for (i in 1:nrow(DataFrame))
{ # i = 1
Temp = unlist(strsplit(DataFrame$WholeGenome[i], split = ''))
FreqOfNuc=data.frame(table(Temp))
DataFrame$NucA[i] = FreqOfNuc[FreqOfNuc$Temp == 'a',]$Freq
DataFrame$NucT[i] = FreqOfNuc[FreqOfNuc$Temp == 't',]$Freq
DataFrame$NucG[i] = FreqOfNuc[FreqOfNuc$Temp == 'g',]$Freq
DataFrame$NucC[i] = FreqOfNuc[FreqOfNuc$Temp == 'c',]$Freq

TempWithShift = c('ZERO',Temp[1:(length(Temp)-1)])
data = data.frame(cbind(Temp,TempWithShift))
data=data[-1,]
data$TwoNeighbors = paste(data$Temp,data$TempWithShift,sep='')
FreqOfDinucl = data.frame(table(data$TwoNeighbors))
DataFrame$DinucAA[i] = FreqOfDinucl[FreqOfDinucl$Var1 == 'aa',]$Freq
DataFrame$DinucTT[i] = FreqOfDinucl[FreqOfDinucl$Var1 == 'tt',]$Freq
DataFrame$DinucGG[i] = FreqOfDinucl[FreqOfDinucl$Var1 == 'gg',]$Freq
DataFrame$DinucCC[i] = FreqOfDinucl[FreqOfDinucl$Var1 == 'cc',]$Freq
}

# let's assume that ExpAA = fr(A)*fr(A)*GenomeLength
DataFrame$ObsExpAA = DataFrame$DinucAA/(((DataFrame$NucA/(DataFrame$NucA + DataFrame$NucT + DataFrame$NucG + DataFrame$NucC))^2)*(DataFrame$NucA + DataFrame$NucT + DataFrame$NucG + DataFrame$NucC))
DataFrame$ObsExpTT = DataFrame$DinucTT/(((DataFrame$NucT/(DataFrame$NucA + DataFrame$NucT + DataFrame$NucG + DataFrame$NucC))^2)*(DataFrame$NucA + DataFrame$NucT + DataFrame$NucG + DataFrame$NucC))
DataFrame$ObsExpGG = DataFrame$DinucGG/(((DataFrame$NucG/(DataFrame$NucA + DataFrame$NucT + DataFrame$NucG + DataFrame$NucC))^2)*(DataFrame$NucA + DataFrame$NucT + DataFrame$NucG + DataFrame$NucC))
DataFrame$ObsExpCC = DataFrame$DinucCC/(((DataFrame$NucC/(DataFrame$NucA + DataFrame$NucT + DataFrame$NucG + DataFrame$NucC))^2)*(DataFrame$NucA + DataFrame$NucT + DataFrame$NucG + DataFrame$NucC))

## Amphibia have excess of AA and TT  (as compared to Lepidosauria and Primates) which is the potential result of UV radiation 
summary(DataFrame$ObsExpAA) # 1.0068 (Lep); 1.0395 (Amp); 1.0076 (Primates)
summary(DataFrame$ObsExpTT) # 1.0861 (Lep); 1.1196 (Amp); 1.0169 (Primates)
summary(DataFrame$ObsExpGG) # 1.611 (Lep);  1.489  (Amp); 1.4848 (Primates)
summary(DataFrame$ObsExpCC) # 1.1058 (Lep); 1.2041 (Amp); 1.1435 (Primates)

## Can UV light go through 