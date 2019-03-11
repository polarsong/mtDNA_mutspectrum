rm(list=ls(all=TRUE))

### reading Input

Data = read.table("../../Body/2Derived/Final_age.1.txt")

#Alt = Data[Data$AltFreqThreshold < 0.3,] # I expect that alternative alleles are rare.
#Alt = Data[Data$AltFreqThreshold < 0.5,] # I expect that alternative alleles are rare. 
Alt = Data[Data$AltFreqThreshold < 1,] # I expect that alternatie alleles are rare. 

Agg = aggregate(list(Alt$FreqD1,Alt$FreqD2,Alt$FreqD3,Alt$FreqD4,Alt$FreqD5,Alt$FreqD6,Alt$FreqD7,Alt$FreqD8,Alt$FreqD9,Alt$FreqD10), by = list(Alt$Subs), FUN = sum)
names(Agg) = c('Subs','FreqD1','FreqD2','FreqD3','FreqD4','FreqD5','FreqD6','FreqD7','FreqD8','FreqD9','FreqD10')
Agg$Subs
AggT = data.frame(t(Agg[,c(2:11)])); names(AggT) = Agg$Subs
AggT$Age = c(1:10)
cor.test(AggT$G_A,AggT$Age, method = 'spearman') # positive trend
cor.test(AggT$T_C,AggT$Age, method = 'spearman') # positive trend
cor.test((AggT$G_A + AggT$T_C),AggT$Age, method = 'spearman') # positive trend
cor.test(AggT$A_G,AggT$Age, method = 'spearman') # negative trend
cor.test(AggT$C_T,AggT$Age, method = 'spearman') # negative trend
cor.test((AggT$C_T + AggT$A_G),AggT$Age, method = 'spearman') # negative trend

colnames(AggT)

AggT$Tv = (AggT$A_C + AggT$C_A + AggT$A_T + AggT$T_A + AggT$G_C + AggT$C_G + AggT$G_T + AggT$T_G)
AggT$Ts = (AggT$G_A + AggT$A_G + AggT$T_C + AggT$C_T)

cor.test(AggT$Ts,AggT$Age, method = 'spearman') # zero!
cor.test(AggT$Tv,AggT$Age, method = 'spearman') # a bit positive
