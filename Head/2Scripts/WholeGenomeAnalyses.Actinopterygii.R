rm(list=ls(all=TRUE))

library(ggplot2)

MUT = read.table('/home/diliushchenko/mtDNA_mutspectrum/Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt', header = TRUE)

TEMPE = read.table('/home/diliushchenko/mtDNA_mutspectrum/Body/1Raw/FishBaseTemperature.txt', header = TRUE)
class(TEMPE$Temperature)
class(TEMPE$Species)
TEMPE = aggregate(Temperature ~ ., median, data = TEMPE)
TemperMut = merge(MUT, TEMPE)

taxa = read.table("/home/diliushchenko/mtDNA_mutspectrum/Body/1Raw/TaxaFromKostya.Names.stat", sep = '\t',header = FALSE) 
taxa = as.data.frame(taxa[grepl('Actinopteri',taxa$V1),])
names(taxa) = 'taxa'

taxa$Species = gsub(";(.*)",'',taxa$taxa);
taxa$Species = gsub(" ",'_',taxa$Species);
taxa$Family = gsub(";Actinopteri(.*)",'',taxa$taxa)
taxa$Family = gsub("(.*);",'',taxa$Family)
table(taxa$Family) ### some species with freq of families > 3
length(unique(taxa$Family)) ### 43 Families 
taxa = taxa[,-1]


### merge Taxa with Temp and MutSpec

finaldf = merge(TemperMut, taxa, by = 'Species')
nrow(finaldf) # 110 sps with all inf
for (i in 1:nrow(finaldf))   {  finaldf$FamilyShort[i] = paste(unlist(strsplit(finaldf$Family[i],''))[c(1:3)],collapse = '')  } # make Short names
table(finaldf$Family)
length(unique(finaldf$Family))# 30 families

FamFreq = data.frame(table(finaldf$Family));
FrequentFamilies = FamFreq[FamFreq$Freq >= 3,]$Var1; length(FrequentFamilies) # delete families, that have less than 3 species!!!
FamFreq = FamFreq[FamFreq$Var1 %in% FrequentFamilies,]
names(FamFreq)=c('Family','NumberOfSpecies')

Fishes = finaldf[finaldf$Family %in% FrequentFamilies,]
Fishes$TsTv = (Fishes$T_C + Fishes$C_T + Fishes$G_A + Fishes$A_G) / (Fishes$T_A + Fishes$A_T + Fishes$G_C + Fishes$C_G + Fishes$G_T + Fishes$T_G + Fishes$C_A + Fishes$A_C)

### Calculate Ts/Tv and A<G/T<C

Fishes$A_G.T_C = (Fishes$T_C/Fishes$A_G)
Fishes = Fishes[Fishes$A_G.T_C > 0,]
Fishes = Fishes[Fishes$A_G.T_C < Inf,]
Fishes = Fishes[Fishes$TsTv < Inf,]

agg = aggregate(list(Fishes$TsTv,Fishes$A_G.T_C,Fishes$Temperature), by = list(Fishes$Family), FUN = median)

names(agg) = c('Family','TsTv','A_G.T_C','Temperature')

cor.test(agg$TsTv,agg$Temperature,method = 'spearman')
a =cor.test(agg$A_G.T_C,agg$Temperature,method = 'spearman')

(a$estimate)^2
ggplot(data = agg, aes(x = Temperature, y = A_G.T_C))+
  geom_point()+
  geom_smooth(method="lm", se=T, col = 'red')+
  theme_classic()+
  labs(x = 'Mean annual water temperature, °C', y = 'A>G/T>C')

ggplot(data = agg, aes(x = Temperature, y = TsTv))+
  geom_point()+
  geom_smooth(method="lm", se=T)+
  theme_classic()+
  labs(x = 'Mean annual water temperature, °C', y = 'TsTv')



#plots for Families
ggplot(data = Fishes, aes(x = Temperature, y = A_G.T_C, group = FamilyShort, fill = FamilyShort))+
  geom_boxplot(alpha=0.3)+
  theme_classic()+
  labs(x = 'Temperature, °C', y = 'A>G/T>C')
  

ggplot(data = Fishes, aes(x = Temperature, y = TsTv, group = FamilyShort, fill = FamilyShort))+
  geom_boxplot(alpha=0.3)+
  theme(legend.position="none")+
  theme_classic()+
  labs(x = 'Temperature, °C', y = 'TsTv')



##### boxplots by quartiles
boxplot(Fishes[Fishes$Temperature<=quantile(Fishes$Temperature,0.25),]$TsTv,
        Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.25) & Fishes$Temperature<=quantile(Fishes$Temperature,0.5),]$TsTv,
        Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.5) & Fishes$Temperature<=quantile(Fishes$Temperature,0.75),]$TsTv,
        Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.75),]$TsTv,
        names=c('1','2','3','4'), outline = FALSE, notch = TRUE, ylab = 'TsTv')

boxplot(Fishes[Fishes$Temperature<=quantile(Fishes$Temperature,0.25),]$A_G.T_C,
        Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.25) & Fishes$Temperature<=quantile(Fishes$Temperature,0.5),]$A_G.T_C,
        Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.5) & Fishes$Temperature<=quantile(Fishes$Temperature,0.75),]$A_G.T_C,
        Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.75),]$A_G.T_C,
        names=c('1','2','3','4'), outline = FALSE, notch = TRUE, ylab = 'A>G/T>C')



wilcox.test(
  Fishes[Fishes$Temperature<=quantile(Fishes$Temperature,0.25),]$TsTv,
  Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.25) & Fishes$Temperature<=quantile(Fishes$Temperature,0.5),]$TsTv)
wilcox.test(
  Fishes[Fishes$Temperature<=quantile(Fishes$Temperature,0.25),]$TsTv,
  Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.5) & Fishes$Temperature<=quantile(Fishes$Temperature,0.75),]$TsTv)
wilcox.test(
  Fishes[Fishes$Temperature<=quantile(Fishes$Temperature,0.25),]$TsTv,
  Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.75),]$TsTv)
wilcox.test(
  Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.25) & Fishes$Temperature<=quantile(Fishes$Temperature,0.5),]$TsTv,
  Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.5) & Fishes$Temperature<=quantile(Fishes$Temperature,0.75),]$TsTv)
wilcox.test(
  Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.25) & Fishes$Temperature<=quantile(Fishes$Temperature,0.5),]$TsTv,
  Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.75),]$TsTv)


##### boxplots by median
boxplot(Fishes[Fishes$Temperature<=quantile(Fishes$Temperature,0.5),]$TsTv,
        Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.5),]$TsTv,
        names=c('colder fishes','warmer fishes'), outline = FALSE, notch = TRUE, cex = 3, ylab = 'TsTv')

boxplot(Fishes[Fishes$Temperature<=quantile(Fishes$Temperature,0.5),]$A_G.T_C,
        Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.5),]$A_G.T_C,
        names=c('colder fishes','warmer fishes'), outline = FALSE, notch = TRUE, cex = 3, ylab = 'A>G/T>C')


quantile(Fishes$Temperature,0.5) # 15 days

wilcox.test(Fishes[Fishes$Temperature<=quantile(Fishes$Temperature,0.5),]$TsTv,Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.5),]$TsTv) # p-value = 0,0755
wilcox.test(Fishes[Fishes$Temperature<=quantile(Fishes$Temperature,0.5),]$A_G.T_C,Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.5),]$A_G.T_C) # 1.687e-06


### reading whole genomes database

unzip("/home/diliushchenko/mtDNA_mutspectrum/Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip")
SynNuc = read.table("/home/diliushchenko/mtDNA_mutspectrum/Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")) {file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")}

SynNuc = SynNuc[SynNuc$Gene != 'ND6',]

####### obtaining neutral nucleotide fractions
SynNuc = aggregate(list(SynNuc$NeutralA,SynNuc$NeutralT,SynNuc$NeutralG,SynNuc$NeutralC), by = list(SynNuc$Species), FUN = sum)
names(SynNuc) = c('Species','NeutralA','NeutralT','NeutralG','NeutralC')
SynNuc$FrA = SynNuc$NeutralA / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC)
SynNuc$FrT = SynNuc$NeutralT / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 
SynNuc$FrG = SynNuc$NeutralG / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 
SynNuc$FrC = SynNuc$NeutralC / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 

### merge whole genomes with temperature
TEMPE = read.table('/home/diliushchenko/mtDNA_mutspectrum/Body/1Raw/FishBaseTemperature.txt', header = TRUE)
TEMPE = aggregate(Temperature ~ ., median, data = TEMPE); summary(TEMPE$Temperature)
SynNuc = merge(TEMPE,SynNuc, by = 'Species'); summary(SynNuc$Temperature)

###### merge whole genomes and temperature with time of maturation
MATUTM = read.table('/home/diliushchenko/mtDNA_mutspectrum/Body/1Raw/FishBaseMaturity_Tm.txt',  header = TRUE)
MATUTM = aggregate(Tm ~ ., median, data = MATUTM) 
SynNuc = merge(MATUTM,SynNuc, by = 'Species'); nrow(SynNuc)

SynNuc$CtoTSkew = (SynNuc$FrC-SynNuc$FrT)/(SynNuc$FrC+SynNuc$FrT); summary(SynNuc$CtoTSkew) 
SynNuc$GtoASkew = (SynNuc$FrG-SynNuc$FrA)/(SynNuc$FrG+SynNuc$FrA); summary(SynNuc$GtoASkew)
SynNuc$AtoGSkew = (SynNuc$FrA-SynNuc$FrG)/(SynNuc$FrA+SynNuc$FrG); summary(SynNuc$AtoGSkew)
SynNuc$TtoCSkew = (SynNuc$FrT-SynNuc$FrC)/(SynNuc$FrT+SynNuc$FrC); summary(SynNuc$TtoCSkew)
SynNuc$TG = SynNuc$FrT+SynNuc$FrG
SynNuc$AC = SynNuc$FrA+SynNuc$FrC
SynNuc$TG_ACSkew = (SynNuc$TG-SynNuc$AC)/(SynNuc$TG+SynNuc$AC); summary(SynNuc$TG_ACSkew)
SynNuc$AC_TGSkew = -(SynNuc$TG-SynNuc$AC)/(SynNuc$TG+SynNuc$AC); summary(SynNuc$AC_TGSkew)

medmatur = median(SynNuc$Tm.x)

#### function to separate by Maturated
mmm = function(x)
{
  if (x < medmatur){return(0)}
  else{return(1)}
}

SynNuc$Maturated = apply(as.matrix(SynNuc$Tm.x),1,mmm)

ShortMaturated = SynNuc[SynNuc$Maturated == 0,]

LongMaturated = SynNuc[SynNuc$Maturated == 1,]


boxplot(ShortMaturated[ShortMaturated$Temperature.x<=quantile(ShortMaturated$Temperature.x,0.5),]$AC_TGSkew,
        ShortMaturated[ShortMaturated$Temperature.x>quantile(ShortMaturated$Temperature.x,0.5),]$AC_TGSkew,
        names=c('colder fishes','warmer fishes'), outline = FALSE, notch = TRUE, cex = 3, ylab = 'AC-TG skew', main = 'ShortMaturated')

boxplot(LongMaturated[LongMaturated$Temperature.x<=quantile(LongMaturated$Temperature.x,0.5),]$AC_TGSkew,
        LongMaturated[LongMaturated$Temperature.x>quantile(LongMaturated$Temperature.x,0.5),]$AC_TGSkew,
        names=c('colder fishes','warmer fishes'), outline = FALSE, notch = TRUE, cex = 3, ylab = 'AC-TG skew', main = 'LongMaturated')


