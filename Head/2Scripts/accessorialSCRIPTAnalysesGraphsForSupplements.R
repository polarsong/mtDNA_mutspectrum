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

t1 = wilcox.test(
  Fishes[Fishes$Temperature<=quantile(Fishes$Temperature,0.25),]$A_G.T_C,
  Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.25) & Fishes$Temperature<=quantile(Fishes$Temperature,0.5),]$A_G.T_C) 

t2 = wilcox.test(
  Fishes[Fishes$Temperature<=quantile(Fishes$Temperature,0.25),]$A_G.T_C,
  Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.5) & Fishes$Temperature<=quantile(Fishes$Temperature,0.75),]$A_G.T_C) 

t3 = wilcox.test(
  Fishes[Fishes$Temperature<=quantile(Fishes$Temperature,0.25),]$A_G.T_C,
  Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.75),]$A_G.T_C) 

t4 = wilcox.test(
  Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.25) & Fishes$Temperature<=quantile(Fishes$Temperature,0.5),]$A_G.T_C,
  Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.5) & Fishes$Temperature<=quantile(Fishes$Temperature,0.75),]$A_G.T_C) 

t5 = wilcox.test(
  Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.25) & Fishes$Temperature<=quantile(Fishes$Temperature,0.5),]$A_G.T_C,
  Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.75),]$A_G.T_C) 

t6 = wilcox.test(
  Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.5) & Fishes$Temperature<=quantile(Fishes$Temperature,0.75),]$A_G.T_C,
  Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.75),]$A_G.T_C) 

round(p.adjust(c(t1$p.value,t2$p.value,t3$p.value,t4$p.value,t5$p.value,t6$p.value), method = "bonferroni"), 3)

##### boxplots by median

boxplot(Fishes[Fishes$Temperature<=quantile(Fishes$Temperature,0.5),]$A_G.T_C,
        Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.5),]$A_G.T_C,
        names=c('colder fishes','warmer fishes'), outline = FALSE, notch = TRUE, cex = 3, ylab = 'A>G/T>C')


quantile(Fishes$Temperature,0.5) # 15 days

wilcox.test(Fishes[Fishes$Temperature<=quantile(Fishes$Temperature,0.5),]$A_G.T_C,
            Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.5),]$A_G.T_C) 


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

m_temp_s = median(SynNuc[SynNuc$Maturated == 'Short Maturated',]$Temperature)
m_temp_l = median(SynNuc[SynNuc$Maturated == 'Long Maturated',]$Temperature)

for (i in 1:nrow(SynNuc))
{
  if (SynNuc$Maturated[i] == 'Short Maturated')
  {
    if (SynNuc$Temperature[i] <= m_temp_s) {SynNuc$Temp[i] = 'Colder Fishes'}
    else {SynNuc$Temp[i] = 'Warmer Fishes'}
  }
  else if (SynNuc$Maturated[i] == 'Long Maturated') 
  {
    if (SynNuc$Temperature[i] <= m_temp_l) {SynNuc$Temp[i] = 'Colder Fishes'}
    else {SynNuc$Temp[i] = 'Warmer Fishes'}
  }
}
ggplot(data = SynNuc, aes(x = Temp, y = AC_TGSkew))+
  geom_boxplot()+
  facet_wrap(~Maturated)+
  theme_test()+
  labs(y = 'ACí-TGí skew')
  

t1 = wilcox.test(SynNuc[SynNuc$Maturated =='Long Maturated' & SynNuc$Temp == 'Colder Fishes',]$AC_TGSkew,
            SynNuc[SynNuc$Maturated =='Long Maturated' & SynNuc$Temp == 'Warmer Fishes',]$AC_TGSkew) ### inside maturated p - value = 0.0057

t2 = wilcox.test(SynNuc[SynNuc$Maturated =='Short Maturated' & SynNuc$Temp == 'Colder Fishes',]$AC_TGSkew,
            SynNuc[SynNuc$Maturated =='Short Maturated' & SynNuc$Temp == 'Warmer Fishes',]$AC_TGSkew) ### inside maturated p - value = 0.0039

t3 = wilcox.test(SynNuc[SynNuc$Maturated =='Short Maturated' & SynNuc$Temp == 'Colder Fishes',]$AC_TGSkew,
            SynNuc[SynNuc$Maturated =='Long Maturated' & SynNuc$Temp == 'Colder Fishes',]$AC_TGSkew) ### between maturated just cold p - value = 0.16 

t4 = wilcox.test(SynNuc[SynNuc$Maturated =='Short Maturated' & SynNuc$Temp == 'Warmer Fishes',]$AC_TGSkew,
            SynNuc[SynNuc$Maturated =='Long Maturated' & SynNuc$Temp == 'Warmer Fishes',]$AC_TGSkew) ### between maturated just warm p - value = 0.4545

t5 = wilcox.test(SynNuc[SynNuc$Maturated =='Short Maturated' & SynNuc$Temp == 'Colder Fishes',]$AC_TGSkew,
            SynNuc[SynNuc$Maturated =='Long Maturated' & SynNuc$Temp == 'Warmer Fishes',]$AC_TGSkew) ### between maturated warm and cold p - value = 0.0003

t6 = wilcox.test(SynNuc[SynNuc$Maturated =='Short Maturated' & SynNuc$Temp == 'Warmer Fishes',]$AC_TGSkew,
            SynNuc[SynNuc$Maturated =='Long Maturated' & SynNuc$Temp == 'Colder Fishes',]$AC_TGSkew) ### between maturated warm and cold p - value = 0.04

p.adjust(c(t1$p.value,t2$p.value,t3$p.value,t4$p.value,t5$p.value,t6$p.value), method = "bonferroni")

#### reading df for drawing mammals plots

unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip")
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")) {file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")}

kuptsovtable = read.table("../../Body/2Derived/EcologyMammalianTable01_KuptsovA_ver2_Full.txt", sep='\t', header=TRUE)

SynNuc = SynNuc[SynNuc$Gene != 'ND6',]
SynNuc = aggregate(list(SynNuc$NeutralA,SynNuc$NeutralT,SynNuc$NeutralG,SynNuc$NeutralC), by = list(SynNuc$Species), FUN = sum)
names(SynNuc) = c('Species','NeutralA','NeutralT','NeutralG','NeutralC')
SynNuc$FrA = SynNuc$NeutralA / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC)
SynNuc$FrT = SynNuc$NeutralT / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 
SynNuc$FrG = SynNuc$NeutralG / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 
SynNuc$FrC = SynNuc$NeutralC / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 

kuptsovtable$FrT= NULL
allparameters = merge(kuptsovtable, SynNuc, by="Species")
allparameters$Temper = as.numeric(gsub(",", ".", allparameters$Temperature.C._White2003.2006.other.close.species))
allparameters$GenerationLength_d = as.numeric(gsub(",", ".", allparameters$GenerationLength_d))
allparameters$TG = allparameters$FrT+allparameters$FrG
allparameters$AC = allparameters$FrA+allparameters$FrC
allparameters$TG_ACSkew = (allparameters$TG-allparameters$AC)/(allparameters$TG+allparameters$AC); summary(allparameters$TG_ACSkew)
allparameters$TtoCSkew = (allparameters$FrT-allparameters$FrC)/(allparameters$FrT+allparameters$FrC); summary(allparameters$TtoCSkew)
allparameters$AC_TGSkew = -(allparameters$TG-allparameters$AC)/(allparameters$TG+allparameters$AC); summary(allparameters$AC_TGSkew)

allparameters$MarsMono = allparameters$Mars + allparameters$Mono
table(allparameters$MarsMono)

formediantemperature = allparameters[!is.na(allparameters$Temper),]$Temper
coldspeciesnames = allparameters[allparameters$Temper <= mean(formediantemperature) & !is.na(allparameters$Temper),]$Species
allparameters$colddummy = 0
allparameters[allparameters$Species %in% coldspeciesnames,]$colddummy = 1

allparameters$allcolddummy = allparameters$Hib.unconfirmedHib + allparameters$Daily.unconfirmedDaily + allparameters$MarsMono + allparameters$colddummy
table(allparameters$allcolddummy)
allparameters[allparameters$allcolddummy > 0,]$allcolddummy = 1

for (i in 1:nrow(allparameters)){if(allparameters$allcolddummy[i] == 0) {allparameters$allcolddummy[i] = 'Colder mammals'} else {allparameters$allcolddummy[i] = 'Warmer mammals'}}

medGL = median(allparameters$GenerationLength_d)

allparameters$GLgroups = "Long-lived mammals"
allparameters[allparameters$GenerationLength_d < medGL,]$GLgroups = "Short-lived mammals"
allparameters$gl_f = factor(allparameters$GLgroups, levels = c('Short-lived mammals','Long-lived mammals'))

#pdf('../../Body/4Figures/Supplements_Fig4.pdf')

ggplot(data = allparameters, aes(x = allcolddummy, y = AC_TGSkew))+
  geom_boxplot()+
  facet_wrap(~gl_f)+
  theme_test()+
  labs(x = '',y = 'TGh_ACH skew')

#dev.off()

t1 = wilcox.test(allparameters[allparameters$GLgroups =='Long-lived mammals' & allparameters$allcolddummy == 'Colder mammals',]$AC_TGSkew,
                 allparameters[allparameters$GLgroups =='Long-lived mammals' & allparameters$allcolddummy == 'Warmer mammals',]$AC_TGSkew) ### inside lifespan 

t2 = wilcox.test(allparameters[allparameters$GLgroups =='Short-lived mammals' & allparameters$allcolddummy == 'Colder mammals',]$AC_TGSkew,
                 allparameters[allparameters$GLgroups =='Short-lived mammals' & allparameters$allcolddummy == 'Warmer mammals',]$AC_TGSkew) ### inside lifespan 

t3 = wilcox.test(allparameters[allparameters$GLgroups =='Short-lived mammals' & allparameters$allcolddummy == 'Colder mammals',]$AC_TGSkew,
                 allparameters[allparameters$GLgroups =='Long-lived mammals' & allparameters$allcolddummy == 'Colder mammals',]$AC_TGSkew) ### between lifespan just cold 

t4 = wilcox.test(allparameters[allparameters$GLgroups =='Short-lived mammals' & allparameters$allcolddummy == 'Warmer mammals',]$AC_TGSkew,
                 allparameters[allparameters$GLgroups =='Long-lived mammals' & allparameters$allcolddummy == 'Warmer mammals',]$AC_TGSkew) ### between lifespan just warm 

t5 = wilcox.test(allparameters[allparameters$GLgroups =='Short-lived mammals' & allparameters$allcolddummy == 'Colder mammals',]$AC_TGSkew,
                 allparameters[allparameters$GLgroups =='Long-lived mammals' & allparameters$allcolddummy == 'Warmer mammals',]$AC_TGSkew) ### between lifespan warm and cold 

t6 = wilcox.test(allparameters[allparameters$GLgroups =='Short-lived mammals' & allparameters$allcolddummy == 'Warmer mammals',]$AC_TGSkew,
                 allparameters[allparameters$GLgroups =='Long-lived mammals' & allparameters$allcolddummy == 'Colder mammals',]$AC_TGSkew) ### between lifespan warm and cold 

round(p.adjust(c(t1$p.value,t2$p.value,t3$p.value,t4$p.value,t5$p.value,t6$p.value), method = "bonferroni"), 10)


