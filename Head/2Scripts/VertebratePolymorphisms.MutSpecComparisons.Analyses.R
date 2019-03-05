rm(list=ls(all=TRUE))

MUT = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.txt')

## filter set of species with more than minimal number of mutations (minimum is 15) for different settings
length(unique(MUT$Species))                                    # 2118  
length(unique(MUT[MUT$NumOfSynMut >= 15,]$Species))            # 1560
length(unique(MUT[MUT$NumOfFourFoldMut >= 15,]$Species))       # 1631
length(unique(MUT[MUT$NumOfSynInCytB >= 15,]$Species))         # 1122
length(unique(MUT[MUT$NumOfFourFoldMutInCytB >= 15,]$Species)) # 1172

Mut = MUT[MUT$MutType == 'FourFold',]
  
Mut = Mut[(Mut$NumOfFourFoldMut) >= 15,] 

#### generate MutSpec for each species

## get ancestral nucleotide 
Mut$Subs = as.character(Mut$Subs)
TRIM<-function(x)	 {unlist(strsplit(x,'_'))[1]}
Mut$AncestralNuc = apply(as.matrix(Mut$Subs), 1 , FUN = TRIM)

## normalize by freq of ancestral nucleotides
A = Mut[Mut$AncestralNuc == 'A',]; A$Number = A$A/(A$A+A$T+A$G+A$C)
T = Mut[Mut$AncestralNuc == 'T',]; T$Number = T$T/(T$A+T$T+T$G+T$C)
G = Mut[Mut$AncestralNuc == 'G',]; G$Number = G$G/(G$A+G$T+G$G+G$C)
C = Mut[Mut$AncestralNuc == 'C',]; C$Number = C$C/(C$A+C$T+C$G+C$C)

Mut = rbind(A,T); Mut = rbind(Mut,G); Mut = rbind(Mut,C);
agg = aggregate(Mut$Number, by = list(Mut$Species,Mut$Subs), FUN = sum)
names(agg)=c('Species','Subs','Freq')

## make vector of 12 Subs for each species
Template = data.frame(unique(Mut$Subs)); names(Template) = c('Subs'); Template$Freq = 0;
VecOfSpecies = unique(agg$Species)
for (i in 1:length(VecOfSpecies))
{ # i = 2
  Temp = agg[agg$Species == VecOfSpecies[i],]
  Template$Species = VecOfSpecies[i]
  Temp = merge(Temp,Template, by = c('Subs','Species'), all = TRUE)
  Temp[is.na(Temp)] <- 0
  Temp$Freq = Temp$Freq.x + Temp$Freq.y
  ALL = sum(Temp$Freq)
  
  Line = data.frame(VecOfSpecies[i], Temp[Temp$Subs == 'A_T',]$Freq/ALL, Temp[Temp$Subs == 'A_G',]$Freq/ALL, Temp[Temp$Subs == 'A_C',]$Freq/ALL, Temp[Temp$Subs == 'T_A',]$Freq/ALL, Temp[Temp$Subs == 'T_G',]$Freq/ALL, Temp[Temp$Subs == 'T_C',]$Freq/ALL, Temp[Temp$Subs == 'G_A',]$Freq/ALL, Temp[Temp$Subs == 'G_T',]$Freq/ALL, Temp[Temp$Subs == 'G_C',]$Freq/ALL, Temp[Temp$Subs == 'C_A',]$Freq/ALL, Temp[Temp$Subs == 'C_T',]$Freq/ALL, Temp[Temp$Subs == 'C_G',]$Freq/ALL)
  names(Line)=c('Species','A_T','A_G','A_C','T_A','T_G','T_C','G_A','G_T','G_C','C_A','C_T','C_G')
  if (i == 1) {Final = Line}
  if (i >  1) {Final = rbind(Final,Line)}
}

write.table(Final, '../../Body/3Results/VertebratePolymorphisms.MutSpecData.FourFoldSynMutationsAllGenes.txt', quote = FALSE, row.names = FALSE)
