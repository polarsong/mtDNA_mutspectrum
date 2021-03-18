rm(list=ls(all=TRUE))

library(seqinr)

#### read fasta file

# MutSpec = read.fasta("../../Body/1Raw/AmphibiansMtDnaRefSeqs/sequence.fasta")  # Amphibians
MutSpec = read.fasta("../../Body/1Raw/LepidosauriaMtDnaRefSeqs/sequence.fasta")  # Lepidosauria
# MutSpec = read.fasta("../../Body/1Raw/PrimatesMtDnaRefSeqs/sequence.fasta")  # Primates

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
data$TwoNeighbors = paste(data$TempWithShift,data$Temp,sep='')  # this order will make a sence for light strand (from start to end)
FreqOfDinucl = data.frame(table(data$TwoNeighbors))
DataFrame$DinucAA[i] = FreqOfDinucl[FreqOfDinucl$Var1 == 'aa',]$Freq
DataFrame$DinucAT[i] = FreqOfDinucl[FreqOfDinucl$Var1 == 'at',]$Freq
DataFrame$DinucTA[i] = FreqOfDinucl[FreqOfDinucl$Var1 == 'ta',]$Freq
DataFrame$DinucTT[i] = FreqOfDinucl[FreqOfDinucl$Var1 == 'tt',]$Freq
DataFrame$DinucTG[i] = FreqOfDinucl[FreqOfDinucl$Var1 == 'tg',]$Freq
DataFrame$DinucTC[i] = FreqOfDinucl[FreqOfDinucl$Var1 == 'tc',]$Freq
DataFrame$DinucGT[i] = FreqOfDinucl[FreqOfDinucl$Var1 == 'gt',]$Freq
DataFrame$DinucGG[i] = FreqOfDinucl[FreqOfDinucl$Var1 == 'gg',]$Freq
DataFrame$DinucGC[i] = FreqOfDinucl[FreqOfDinucl$Var1 == 'gc',]$Freq
DataFrame$DinucCG[i] = FreqOfDinucl[FreqOfDinucl$Var1 == 'cg',]$Freq
DataFrame$DinucCC[i] = FreqOfDinucl[FreqOfDinucl$Var1 == 'cc',]$Freq
}

# let's assume that ExpAA = fr(A)*fr(A)*GenomeLength etc..
DataFrame$GenomeLength = DataFrame$NucA + DataFrame$NucT + DataFrame$NucG + DataFrame$NucC
DataFrame$FrA = DataFrame$NucA/DataFrame$GenomeLength
DataFrame$FrT = DataFrame$NucT/DataFrame$GenomeLength
DataFrame$FrG = DataFrame$NucG/DataFrame$GenomeLength
DataFrame$FrC = DataFrame$NucC/DataFrame$GenomeLength

DataFrame$ObsExpAA = DataFrame$DinucAA/(DataFrame$FrA*DataFrame$FrA*DataFrame$GenomeLength)
DataFrame$ObsExpAT = DataFrame$DinucAT/(DataFrame$FrA*DataFrame$FrT*DataFrame$GenomeLength)
DataFrame$ObsExpTA = DataFrame$DinucTA/(DataFrame$FrA*DataFrame$FrT*DataFrame$GenomeLength)

DataFrame$ObsExpTT = DataFrame$DinucTT/(DataFrame$FrT*DataFrame$FrT*DataFrame$GenomeLength)
DataFrame$ObsExpTG = DataFrame$DinucTG/(DataFrame$FrT*DataFrame$FrG*DataFrame$GenomeLength)
DataFrame$ObsExpGT = DataFrame$DinucGT/(DataFrame$FrT*DataFrame$FrG*DataFrame$GenomeLength)
DataFrame$ObsExpTC = DataFrame$DinucTC/(DataFrame$FrT*DataFrame$FrC*DataFrame$GenomeLength)

DataFrame$ObsExpGG = DataFrame$DinucGG/(DataFrame$FrG*DataFrame$FrG*DataFrame$GenomeLength)
DataFrame$ObsExpGC = DataFrame$DinucGC/(DataFrame$FrG*DataFrame$FrC*DataFrame$GenomeLength)
DataFrame$ObsExpCG = DataFrame$DinucCG/(DataFrame$FrG*DataFrame$FrC*DataFrame$GenomeLength)

DataFrame$ObsExpCC = DataFrame$DinucCC/(DataFrame$FrC*DataFrame$FrC*DataFrame$GenomeLength)

## Amphibia have excess of AA and TT  (as compared to Lepidosauria and Primates) which is the potential result of UV radiation 
median(DataFrame$ObsExpAA) # 1.0068 (Lep); 1.0395 (Amp); 1.0076 (Primates)
median(DataFrame$ObsExpTT) # 1.0861 (Lep); 1.1196 (Amp); 1.0169 (Primates)
median(DataFrame$ObsExpGG) # 1.611 (Lep);  1.489  (Amp); 1.4848 (Primates)
median(DataFrame$ObsExpCC) # 1.1058 (Lep); 1.2041 (Amp); 1.1435 (Primates)

### very different !!!!! WHY ????

###### prediction 1:
### CpG > TpG on heavy chain and thus on light chain should be a deficit of CpG
median(DataFrame$ObsExpCG) # 0.6!!!!!!
median(DataFrame$ObsExpGC) # 1.067  - just a control

###### prediction 2:
### ApC > GpC on heavy chain => deficit of GpT on light chain (yes!!)
median(DataFrame$ObsExpGT) # 0.78
median(DataFrame$ObsExpTG) # 1.01 - just a control

###### prediction 3:
### CpC > CpT on heavy chain and thus on light chain should be a deficit of GpG
median(DataFrame$ObsExpGG) # there is an excess (selection? repeat it on neutral sites? D loop, synonymous sites) G is super rare now


## Can UV light go through the egg skin? 