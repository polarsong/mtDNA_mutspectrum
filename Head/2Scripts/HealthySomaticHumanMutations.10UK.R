rm(list=ls(all=TRUE))

### reading Input

Data = read.csv("../../Body/1Raw/12864_2017_4287_MOESM2_ESM.csv")
# Data = read.csv('12864_2017_4287_MOESM2_ESM.csv')
NumberOfIndividuals = length(unique(Data$sampleID)); NumberOfIndividuals

### Adding random age to each sampleID
AgeVector = runif(n=NumberOfIndividuals, min=17, max=85)
Age = data.frame(unique(Data$sampleID),AgeVector); names(Age)=c('sampleID','Age')
Data = merge(Data,Age)

### Analyses
Data$Subs = paste(Data$reference_allele,Data$alternative_allele, sep = '_')
VecOfSubs = unique(Data$Subs); length(VecOfSubs) # 12
ALL = data.frame(VecOfSubs); names(ALL) = c('Subs')

AltFreqThresholds = seq(0,0.9, by=0.1)
for (i in 1:length(AltFreqThresholds))
{ # i = 1
Input = Data[Data$alt_frequency > AltFreqThresholds[i] & Data$alt_frequency <= (AltFreqThresholds[i]+0.1),]
D1 = data.frame(table(Input[Input$Age < quantile(Input$Age,0.1),]$Subs));                                         names(D1) =  c('Subs','FreqD1')
D2 = data.frame(table(Input[Input$Age >= quantile(Input$Age,0.1) & Input$Age <  quantile(Input$Age,0.2),]$Subs)); names(D2) =  c('Subs','FreqD2')
D3 = data.frame(table(Input[Input$Age >= quantile(Input$Age,0.2) & Input$Age <  quantile(Input$Age,0.3),]$Subs)); names(D3) =  c('Subs','FreqD3')
D4 = data.frame(table(Input[Input$Age >= quantile(Input$Age,0.3) & Input$Age <  quantile(Input$Age,0.4),]$Subs)); names(D4) =  c('Subs','FreqD4')
D5 = data.frame(table(Input[Input$Age >= quantile(Input$Age,0.4) & Input$Age <  quantile(Input$Age,0.5),]$Subs)); names(D5) =  c('Subs','FreqD5')
D6 = data.frame(table(Input[Input$Age >= quantile(Input$Age,0.5) & Input$Age <  quantile(Input$Age,0.6),]$Subs)); names(D6) =  c('Subs','FreqD6')
D7 = data.frame(table(Input[Input$Age >= quantile(Input$Age,0.6) & Input$Age <  quantile(Input$Age,0.7),]$Subs)); names(D7) =  c('Subs','FreqD7')
D8 = data.frame(table(Input[Input$Age >= quantile(Input$Age,0.7) & Input$Age <  quantile(Input$Age,0.8),]$Subs)); names(D8) =  c('Subs','FreqD8')
D9 = data.frame(table(Input[Input$Age >= quantile(Input$Age,0.8) & Input$Age <  quantile(Input$Age,0.9),]$Subs)); names(D9) =  c('Subs','FreqD9')
D10 =data.frame(table(Input[Input$Age >= quantile(Input$Age,0.9),]$Subs));                                        names(D10) =  c('Subs','FreqD10')

RES = merge(ALL,D1, all.x = TRUE); RES = merge(RES,D2, all.x = TRUE); RES = merge(RES,D3, all.x = TRUE); RES = merge(RES,D4, all.x = TRUE); RES = merge(RES,D5, all.x = TRUE);
RES = merge(RES,D6, all.x = TRUE); RES = merge(RES,D7, all.x = TRUE); RES = merge(RES,D8, all.x = TRUE); RES = merge(RES,D9, all.x = TRUE); RES = merge(RES,D10, all.x = TRUE);

RES[is.na(RES)] <- 0
RES$AltFreqThreshold = AltFreqThresholds[i]
if (i == 1) {Final = RES}
if (i >  1) {Final = rbind(Final,RES)}
}
write.table(Final, 'Final.txt')

