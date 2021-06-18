rm(list=ls(all=TRUE))

library(ggplot2)

MUT = read.table('C:/Users/MitoClubHouse/Documents/data/DIliushchenko/mtDNA_mutspectrum/Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt', header = TRUE)

AA = read.table("C:/Users/MitoClubHouse/Documents/data/DIliushchenko/mtDNA_mutspectrum/Body/1Raw/anage_data.txt", header = TRUE, sep = '\t')
AA$Species = paste(AA$Genus,AA$Species,sep = '_')
df = merge(MUT,AA)
df = df[df$Class == 'Actinopterygii',]

taxacl = read.table('C:/Users/MitoClubHouse/Documents/data/DIliushchenko/mtDNA_mutspectrum/Body/2Derived/TaxaWithClasses.txt')
taxamut = merge(MUT,taxacl)
taxamut = taxamut[taxamut$Class == 'Actinopterygii',]
taxamut =  merge(AA, taxamut, all.y = TRUE)


TEMPE = read.table('C:/Users/MitoClubHouse/Documents/data/DIliushchenko/mtDNA_mutspectrum/Body/1Raw/FishBaseTemperature.txt', header = TRUE)
class(TEMPE$Temperature)
class(TEMPE$Species)
TEMPE = aggregate(Temperature ~ ., median, data = TEMPE)
TemperMut = merge(MUT, TEMPE)

taxa = read.table("DIliushchenko/mtDNA_mutspectrum/Body/1Raw/TaxaFromKostya.Names.stat", sep = '\t',header = FALSE) 
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
cor.test(agg$A_G.T_C,agg$Temperature,method = 'spearman')


ggplot(data = agg, aes(x = Temperature, y = A_G.T_C))+
  geom_point()+
  geom_smooth(method="lm", se=T)

ggplot(data = agg, aes(x = Temperature, y = TsTv))+
  geom_point()+
  geom_smooth(method="lm", se=T)


#plots for Families
ggplot(data = Fishes, aes(x = Temperature, y = A_G.T_C, group = FamilyShort, fill = FamilyShort))+
  geom_boxplot(alpha=0.3)+
  theme(legend.position="none")
  

ggplot(data = Fishes, aes(x = Temperature, y = TsTv, group = FamilyShort, fill = FamilyShort))+
  geom_boxplot(alpha=0.3)+
  theme(legend.position="none")



##### boxplots by quartiles
boxplot(Fishes[Fishes$Temperature<=quantile(Fishes$Temperature,0.25),]$TsTv,
        Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.25) & Fishes$Temperature<=quantile(Fishes$Temperature,0.5),]$TsTv,
        Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.5) & Fishes$Temperature<=quantile(Fishes$Temperature,0.75),]$TsTv,
        Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.75),]$TsTv,
        names=c('1','2','3','4'), outline = FALSE, notch = TRUE)

boxplot(Fishes[Fishes$Temperature<=quantile(Fishes$Temperature,0.25),]$A_G.T_C,
        Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.25) & Fishes$Temperature<=quantile(Fishes$Temperature,0.5),]$A_G.T_C,
        Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.5) & Fishes$Temperature<=quantile(Fishes$Temperature,0.75),]$A_G.T_C,
        Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.75),]$A_G.T_C,
        names=c('1','2','3','4'), outline = FALSE, notch = TRUE)



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
        names=c('cold','warm'), outline = FALSE, notch = TRUE, cex = 3)

boxplot(Fishes[Fishes$Temperature<=quantile(Fishes$Temperature,0.5),]$A_G.T_C,
        Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.5),]$A_G.T_C,
        names=c('cold','warm'), outline = FALSE, notch = TRUE, cex = 3)


quantile(Fishes$Temperature,0.5) # 15 days

wilcox.test(Fishes[Fishes$Temperature<=quantile(Fishes$Temperature,0.5),]$TsTv,Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.5),]$TsTv) # p-value = 0,0755
wilcox.test(Fishes[Fishes$Temperature<=quantile(Fishes$Temperature,0.5),]$A_G.T_C,Fishes[Fishes$Temperature>quantile(Fishes$Temperature,0.5),]$A_G.T_C) # 1.687e-06

