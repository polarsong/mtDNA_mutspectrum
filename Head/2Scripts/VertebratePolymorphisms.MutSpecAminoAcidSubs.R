#####  COUNT NUCLEOTIDE CONTENT CAREFULLY (BODY/2Derived/polarizedbr_data => external + More Shallow => codons, 4fold nucl, FrA,T,G,C) => barplot?
#####  normalized average MutSpec (pie charts or 12 boxplots for each class)?

rm(list=ls(all=TRUE))

### neutral ATGC
NeutralATGC = read.table('../../Body/3Results/VertebratePolymorphisms.Normalization.NeutralATGC.txt', header = TRUE)

####### READ
MUT = read.table("../../Body/3Results/Mutational_spectra_in_Chordata_ML.txt", header = TRUE)
length(unique(MUT$Species)) # 2404  SOME SPECIES HAVE THREE WORDS => CUT THE LAST AND MERGE WITH TAXA, OR EVEN CUT TWO LAST AND LEAVE JUST GENUS

##### FILTER 1: to take only normal substitutions and filter out species with too high fraction (> 5%) of unnormal substitutions
VecOfNormalSubstitutions <- c('A_C','C_A',
                              'A_G','G_A',
                              'C_G','G_C',
                              'C_T','T_C',
                              'G_T','T_G',
                              'T_A','A_T')
nrow(MUT)
table(MUT$Subs)   # MANY CRAPPY SUBSTITUTIONS!!!!!!!!!!!!!!!!!! WHY?????????????????????
SP = data.frame(table(MUT$Species)); names(SP) = c('Species','NumberOfAllSubst')
SPN = data.frame(table(MUT[MUT$Subs %in% VecOfNormalSubstitutions,]$Species)); names(SPN) = c('Species','NumberOfNormalSubst')
SP = merge(SP,SPN); SP$FractionOfNormal = SP$NumberOfNormalSubst/SP$NumberOfAllSubst
hist(SP$FractionOfNormal)
summary(SP$FractionOfNormal) # how many to delete? ask to have more than 95% of substitutions as normal
SpeciesToDelete = SP[SP$FractionOfNormal <=0.95,]$Species; length(SpeciesToDelete) # 279 - delete
MUT = MUT[!MUT$Species %in% SpeciesToDelete,]
MUT = MUT[MUT$Subs %in% VecOfNormalSubstitutions,]



AACart = read.table("../../Body/1Raw/AminoAcidTranslation.txt")

##### FILTER 2: Type of Substitutions - NonSyn 
MUT$MutTypeAA = "S"
MUT[MUT$AncestralAA != MUT$DescendantAA,]$MutTypeAA = "N"
MUT = MUT[MUT$AncestralAA != "*",]
MUT = MUT[MUT$DescendantAA != "*",]
NonSMut = MUT[MUT$MutTypeAA == "N",]
NonSMut$AncesralAATriple = "None"
NonSMut$DescendantAATriple = "Nonå"


######## translating single-letter aminoacid code into triple-letter
for (i in 1:nrow(NonSMut)){
  AncestralAA = as.character(NonSMut$AncestralAA[i])
  DescendantAA = as.character(NonSMut$DescendantAA[i])
  for (y in 1:nrow(AACart)){
    Triple = as.character(AACart$V1[y])
    Single = as.character(AACart$V2[y])
    if (AncestralAA == Single){
      NonSMut$AncesralAATriple[i]= Triple
    }
    if (DescendantAA == Single){
      NonSMut$DescendantAATriple[i]= Triple
    }
  }
}

length(table(NonSMut$DescendantAATriple))

gooddata = NonSMut

###### Making Leu & Ser aminoacids into four different types depending on codon usage
for (i in 1:nrow(NonSMut)){ # i = 1
  if (NonSMut$AncesralAATriple[i] == 'Leu')
  {
    if (tolower(NonSMut$AncestorCodon[i]) == 'ctt' | tolower(NonSMut$AncestorCodon[i]) == 'ctc' | tolower(NonSMut$AncestorCodon[i]) == 'cta' | tolower(NonSMut$AncestorCodon[i]) == 'ctg') {NonSMut$AncesralAATriple[i] = 'LeuCT'}
    if (tolower(NonSMut$AncestorCodon[i]) == 'tta' | tolower(NonSMut$AncestorCodon[i]) == 'ttg') {NonSMut$AncesralAATriple[i] = 'LeuTT'}
  }
  if (NonSMut$AncesralAATriple[i] == 'Ser')
  {
    if (tolower(NonSMut$AncestorCodon[i]) == 'tct' | tolower(NonSMut$AncestorCodon[i]) == 'tcc' | tolower(NonSMut$AncestorCodon[i]) == 'tca' | tolower(NonSMut$AncestorCodon[i]) == 'tcg') {NonSMut$AncesralAATriple[i] = 'SerTC'}
    if (tolower(NonSMut$AncestorCodon[i]) == 'agt' | tolower(NonSMut$AncestorCodon[i]) == 'agc') {NonSMut$AncesralAATriple[i] = 'SerAG'}
  }
}
table(NonSMut$AncesralAATriple)

for (i in 1:nrow(NonSMut)){ # i = 1
  if (NonSMut$DescendantAATriple[i] == 'Leu')
  {
    if (tolower(NonSMut$DescendantCodon[i]) == 'ctt' | tolower(NonSMut$DescendantCodon[i]) == 'ctc' | tolower(NonSMut$DescendantCodon[i]) == 'cta' | tolower(NonSMut$DescendantCodon[i]) == 'ctg') {NonSMut$DescendantAATriple[i] = 'LeuCT'}
    if (tolower(NonSMut$DescendantCodon[i]) == 'tta' | tolower(NonSMut$DescendantCodon[i]) == 'ttg') {NonSMut$DescendantAATriple[i] = 'LeuTT'}
  }
  if (NonSMut$DescendantAATriple[i] == 'Ser')
  {
    if (tolower(NonSMut$DescendantCodon[i]) == 'tct' | tolower(NonSMut$DescendantCodon[i]) == 'tcc' | tolower(NonSMut$DescendantCodon[i]) == 'tca' | tolower(NonSMut$DescendantCodon[i]) == 'tcg') {NonSMut$DescendantAATriple[i] = 'SerTC'}
    if (tolower(NonSMut$DescendantCodon[i]) == 'agt' | tolower(NonSMut$DescendantCodon[i]) == 'agc') {NonSMut$DescendantAATriple[i] = 'SerAG'}
  }
}

length(table(NonSMut$DescendantAATriple))



### counting numbers of AA substitutions in each species
InterData = data.frame(NonSMut$Species, NonSMut$Gene, NonSMut$AncesralAATriple, NonSMut$DescendantAATriple,
                         Number = 1)
InterData$AASub = paste(InterData$NonSMut.AncesralAATriple, InterData$NonSMut.DescendantAATriple, sep = "_")
TypesOfSub = data.frame(table(InterData$AASub))
TypesOfSub = TypesOfSub$Var1
InterData = InterData[InterData$NonSMut.Gene != "ND6",]

Final = aggregate(InterData$Number, by = list (InterData$NonSMut.Species, InterData$NonSMut.Gene, InterData$AASub), FUN = sum)
names(Final) = c("Species", "Gene", "TypesOfAASub", "FreqOfSub")

Taxa =  read.table('../../Body/3Results/TaxaWithClasses.txt', header = TRUE)

Final = merge(Final, Taxa, all.x = TRUE)
table(Final$Class)


write.table(Final, file = '../../Body/3Results/VertebratePolymorphisms.MutSpecAminoAcidSubs.txt', quote = FALSE)


