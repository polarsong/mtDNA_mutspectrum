rm(list=ls(all=TRUE))

library(ggplot2)

MUT = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt', header = TRUE)

TEMPE = read.table('../../Body/1Raw/FishBaseTemperature.txt', header = TRUE)
class(TEMPE$Temperature)
class(TEMPE$Species)
TEMPE = aggregate(Temperature ~ ., median, data = TEMPE)
TemperMut = merge(MUT, TEMPE)

taxa = read.table("../../Body/1Raw/TaxaFromKostya.Names.stat", sep = '\t',header = FALSE) 
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
for (i in 1:nrow(agg))   {  agg$FamilyShort[i] = paste(unlist(strsplit(agg$Family[i],''))[c(1:3)],collapse = '')  }

cor.test(agg$A_G.T_C, agg$Temperature,method = 'spearman') ### positive and good p - value = 0.0227

svg("../../Body/4Figures/.svg")
ggplot(data = agg, aes(x = Temperature, y = A_G.T_C, label = FamilyShort, color = Family))+
  geom_point()+
  geom_smooth(method="lm", se=F, col = 'red')+
  theme_classic()+
  labs(x = 'Median annual water temperature, °C', y = 'Median A>G/T>C')+
  geom_text(aes(label=FamilyShort),hjust=-0.15, vjust=-0.5, show.legend = F)
  

#plots for Families
ggplot(data = Fishes, aes(x = Temperature, y = A_G.T_C, group = FamilyShort, fill = FamilyShort))+
  geom_boxplot(alpha=0.3)+
  theme_classic()+
  labs(x = 'Temperature, Â°C', y = 'A>G/T>C')
  

ggplot(data = Fishes, aes(x = Temperature, y = TsTv, group = FamilyShort, fill = FamilyShort))+
  geom_boxplot(alpha=0.3)+
  theme(legend.position="none")+
  theme_classic()+
  labs(x = 'Temperature, Â°C', y = 'TsTv')



##### boxplots by quartiles

boxplot(Fishes[Fishes$Temperature<=quantile(Fishes$Temperature,0.25),]$A_G.T_C,
        Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.25) & Fishes$Temperature<=quantile(Fishes$Temperature,0.5),]$A_G.T_C,
        Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.5) & Fishes$Temperature<=quantile(Fishes$Temperature,0.75),]$A_G.T_C,
        Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.75),]$A_G.T_C,
        names=c('1','2','3','4'), outline = FALSE, notch = FALSE, ylab = 'A>G/T>C')



wilcox.test(
  Fishes[Fishes$Temperature<=quantile(Fishes$Temperature,0.25),]$A_G.T_C,
  Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.25) & Fishes$Temperature<=quantile(Fishes$Temperature,0.5),]$A_G.T_C) # p -value 0.025

wilcox.test(
  Fishes[Fishes$Temperature<=quantile(Fishes$Temperature,0.25),]$A_G.T_C,
  Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.5) & Fishes$Temperature<=quantile(Fishes$Temperature,0.75),]$A_G.T_C) # p-value = 0.006

wilcox.test(
  Fishes[Fishes$Temperature<=quantile(Fishes$Temperature,0.25),]$A_G.T_C,
  Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.75),]$A_G.T_C) # p-value = almost 0.000

wilcox.test(
  Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.25) & Fishes$Temperature<=quantile(Fishes$Temperature,0.5),]$A_G.T_C,
  Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.5) & Fishes$Temperature<=quantile(Fishes$Temperature,0.75),]$A_G.T_C) # p-value = 0.93

wilcox.test(
  Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.25) & Fishes$Temperature<=quantile(Fishes$Temperature,0.5),]$A_G.T_C,
  Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.75),]$A_G.T_C) # p-value = 0.016


##### boxplots by median

boxplot(Fishes[Fishes$Temperature<=quantile(Fishes$Temperature,0.5),]$A_G.T_C,
        Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.5),]$A_G.T_C,
        names=c('colder fishes','warmer fishes'), outline = FALSE, notch = TRUE, cex = 3, ylab = 'A>G/T>C')


quantile(Fishes$Temperature,0.5) # 15 days

wilcox.test(Fishes[Fishes$Temperature<=quantile(Fishes$Temperature,0.5),]$A_G.T_C,
            Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.5),]$A_G.T_C) # p-value = 0.00026


### reading whole genomes database

unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
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
TEMPE = read.table('../../Body/1Raw/FishBaseTemperature.txt', header = TRUE)
TEMPE = aggregate(Temperature ~ ., median, data = TEMPE); summary(TEMPE$Temperature)
SynNuc = merge(TEMPE,SynNuc, by = 'Species'); summary(SynNuc$Temperature)

###### merge whole genomes and temperature with time of maturation
MATUTM = read.table('../../Body/1Raw/FishBaseMaturity_Tm.txt',  header = TRUE)
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

medmatur = median(SynNuc$Tm)

#### function to separate by Maturated
mmm = function(x)
{
  if (x < medmatur){return('Short Maturated')}
  else{return('Long Maturated')}
}

SynNuc$Maturated = apply(as.matrix(SynNuc$Tm),1,mmm)

for ( i in 1:nrow(SynNuc)){
  if (SynNuc$Temperature.x[i]<=median(SynNuc$Temperature.x)){SynNuc$Temp[i] = 'Colder Fishes'}
  else {SynNuc$Temp[i] = 'Warmer Fishes'}
}

ggplot(data = SynNuc, aes(x = Temp, y = AC_TGSkew))+
  geom_boxplot()+
  facet_wrap(~Maturated)+
  theme_test()+
  labs(y = 'ACí-TGí skew')
  

wilcox.test(SynNuc[SynNuc$Maturated =='Long Maturated' & SynNuc$Temp == 'Colder Fishes',]$AC_TGSkew,
            SynNuc[SynNuc$Maturated =='Long Maturated' & SynNuc$Temp == 'Warmer Fishes',]$AC_TGSkew) ### inside maturated p - value = 0.0057

wilcox.test(SynNuc[SynNuc$Maturated =='Short Maturated' & SynNuc$Temp == 'Colder Fishes',]$AC_TGSkew,
            SynNuc[SynNuc$Maturated =='Short Maturated' & SynNuc$Temp == 'Warmer Fishes',]$AC_TGSkew) ### inside maturated p - value = 0.0039

wilcox.test(SynNuc[SynNuc$Maturated =='Short Maturated' & SynNuc$Temp == 'Colder Fishes',]$AC_TGSkew,
            SynNuc[SynNuc$Maturated =='Long Maturated' & SynNuc$Temp == 'Colder Fishes',]$AC_TGSkew) ### between maturated just cold p - value = 0.16 

wilcox.test(SynNuc[SynNuc$Maturated =='Short Maturated' & SynNuc$Temp == 'Warmer Fishes',]$AC_TGSkew,
            SynNuc[SynNuc$Maturated =='Long Maturated' & SynNuc$Temp == 'Warmer Fishes',]$AC_TGSkew) ### between maturated just warm p - value = 0.4545

wilcox.test(SynNuc[SynNuc$Maturated =='Short Maturated' & SynNuc$Temp == 'Colder Fishes',]$AC_TGSkew,
            SynNuc[SynNuc$Maturated =='Long Maturated' & SynNuc$Temp == 'Warmer Fishes',]$AC_TGSkew) ### between maturated warm and cold p - value = 0.0003

wilcox.test(SynNuc[SynNuc$Maturated =='Short Maturated' & SynNuc$Temp == 'Warmer Fishes',]$AC_TGSkew,
            SynNuc[SynNuc$Maturated =='Long Maturated' & SynNuc$Temp == 'Colder Fishes',]$AC_TGSkew) ### between maturated warm and cold p - value = 0.04



