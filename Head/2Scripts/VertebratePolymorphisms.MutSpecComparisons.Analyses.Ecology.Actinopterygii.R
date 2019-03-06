rm(list=ls(all=TRUE))


MUT = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt', header = TRUE)

### longevity in fishes

AnAge = read.table('../../Body/1Raw/anage_data.txt', sep = '\t', header = TRUE)
AnAge$Species = paste(AnAge$Genus,AnAge$Species,sep = '_')
AnAge = AnAge[AnAge$Class == 'Actinopterygii',]
AnAgeMut = merge(MUT,AnAge) # 110

cor.test(AnAgeMut$A_T,AnAgeMut$Maximum.longevity..yrs., method = 'spearman')
cor.test(AnAgeMut$A_G,AnAgeMut$Maximum.longevity..yrs., method = 'spearman')
cor.test(AnAgeMut$A_C,AnAgeMut$Maximum.longevity..yrs., method = 'spearman')
cor.test(AnAgeMut$T_A,AnAgeMut$Maximum.longevity..yrs., method = 'spearman')
cor.test(AnAgeMut$T_G,AnAgeMut$Maximum.longevity..yrs., method = 'spearman')
cor.test(AnAgeMut$T_C,AnAgeMut$Maximum.longevity..yrs., method = 'spearman')
cor.test(AnAgeMut$G_A,AnAgeMut$Maximum.longevity..yrs., method = 'spearman')
cor.test(AnAgeMut$G_T,AnAgeMut$Maximum.longevity..yrs., method = 'spearman')
cor.test(AnAgeMut$G_C,AnAgeMut$Maximum.longevity..yrs., method = 'spearman')
cor.test(AnAgeMut$C_A,AnAgeMut$Maximum.longevity..yrs., method = 'spearman')
cor.test(AnAgeMut$C_T,AnAgeMut$Maximum.longevity..yrs., method = 'spearman')
cor.test(AnAgeMut$C_G,AnAgeMut$Maximum.longevity..yrs., method = 'spearman')

### temperature in fishes: T_C is positive!?

Temper = read.table('../../Body/1Raw/TemperatureAllChordataDataSet.txt', sep = '\t', header = TRUE)
Temper = Temper[Temper$Tax == 'fishes',]
Temper$Species = gsub(' ','_',Temper$Species)
TemperMut = merge(MUT,Temper) # 41 just

cor.test(TemperMut$T_C,TemperMut$T..oC., method = 'spearman') # positive!!!! Why????
plot(TemperMut$T_C,TemperMut$T..oC.)
cor.test(TemperMut$T_A,TemperMut$T..oC., method = 'spearman')
cor.test(TemperMut$T_G,TemperMut$T..oC., method = 'spearman')

cor.test(TemperMut$A_T,TemperMut$T..oC., method = 'spearman')
cor.test(TemperMut$A_G,TemperMut$T..oC., method = 'spearman') # a bit negative
cor.test(TemperMut$A_C,TemperMut$T..oC., method = 'spearman')

cor.test(TemperMut$G_A,TemperMut$T..oC., method = 'spearman')
cor.test(TemperMut$G_T,TemperMut$T..oC., method = 'spearman')
cor.test(TemperMut$G_C,TemperMut$T..oC., method = 'spearman')

cor.test(TemperMut$C_A,TemperMut$T..oC., method = 'spearman')
cor.test(TemperMut$C_T,TemperMut$T..oC., method = 'spearman')
cor.test(TemperMut$C_G,TemperMut$T..oC., method = 'spearman')

### temperature in fishes: from cold and tropical waters:
ColdFish = read.table('../../Body/1Raw/cold_water_fishes.txt', sep = '\t', header = FALSE); names(ColdFish)=c('Species')
ColdFish$Species = gsub(' ','_',ColdFish$Species)

WarmFish = read.table('../../Body/1Raw/tropical_water_fishes.txt', sep = '\t', header = FALSE); names(WarmFish)=c('Species')
WarmFish$Species = gsub(' ','_',WarmFish$Species)

boxplot(MUT[MUT$Species %in% ColdFish$Species,]$G_A,MUT[MUT$Species %in% WarmFish$Species,]$G_A)








