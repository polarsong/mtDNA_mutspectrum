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

##### FILTER 2: Synonymous Substitutions
nrow(MUT) # 461215
MUT = MUT[MUT$AncestralAA == MUT$DescendantAA,]; nrow(MUT) # 395157
table(MUT$AncestralAA)
table(MUT$DescendantAA)

######## add several columns: Syn, 4fold;  A T G C neutral fractions, number of syn, number of 4fold per species,  number of syn, number of 4fold per species per CYTB 

##### fourfold degenerate sites
VecOfSynFourFoldDegenerateSites <- c('CTT', 'CTC', 'CTA', 'CTG', 
                                     'GTT', 'GTC', 'GTA', 'GTG', 
                                     'TCT', 'TCC', 'TCA', 'TCG', 
                                     'CCT', 'CCC', 'CCA', 'CCG', 
                                     'ACT', 'ACC', 'ACA', 'ACG', 
                                     'GCT', 'GCC', 'GCA', 'GCG', 
                                     'CGT', 'CGC', 'CGA', 'CGG', 
                                     'GGT', 'GGC', 'GGA', 'GGG')
length(unique(VecOfSynFourFoldDegenerateSites)) # 32

VecOfAllOthersSites <-              c('TTT', 'TTC', 'TTA', 'TTG', 
                                     'ATT', 'ATC', 'ATA', 'ATG', 
                                     'TAT', 'TAC', 'TAA', 'TAG', 
                                     'CAT', 'CAC', 'CAA', 'CAG', 
                                     'AAT', 'AAC', 'AAA', 'AAG', 
                                     'GAT', 'GAC', 'GAA', 'GAG', 
                                     'TGT', 'TGC', 'TGA', 'TGG', 
                                     'AGT', 'AGC', 'AGA', 'AGG')
length(unique(VecOfAllOthersSites)) # 32
length(unique(c(VecOfAllOthersSites,VecOfSynFourFoldDegenerateSites))) # 64

nrow(MUT) # 395157 = 209120 + 186037
MUT4fold = MUT[MUT$AncestorCodon %in% VecOfSynFourFoldDegenerateSites & MUT$DescendantCodon %in% VecOfSynFourFoldDegenerateSites,]; nrow(MUT4fold) # 209120
MUTOthers = MUT[MUT$AncestorCodon %in% VecOfAllOthersSites | MUT$DescendantCodon %in% VecOfAllOthersSites,]; nrow(MUTOthers) # 186037
MUT4fold$MutType = 'FourFold'; MUTOthers$MutType = 'Syn'
MUT = rbind(MUT4fold,MUTOthers)
MUT$Number = 1

# Count number of different mutations
AGG1 = aggregate(MUT[MUT$MutType == 'FourFold',]$Number, by = list(MUT[MUT$MutType == 'FourFold',]$Species), FUN = sum); names(AGG1) = c('Species','NumOfFourFoldMut')
AGG2 = aggregate(MUT[MUT$MutType == 'Syn',]$Number, by = list(MUT[MUT$MutType == 'Syn',]$Species), FUN = sum); names(AGG2) = c('Species','NumOfSynMut')
AGG3 = aggregate(MUT[MUT$MutType == 'FourFold' & MUT$Gene == 'CytB',]$Number, by = list(MUT[MUT$MutType == 'FourFold' & MUT$Gene == 'CytB',]$Species), FUN = sum); names(AGG3) = c('Species','NumOfFourFoldMutInCytB')
AGG4 = aggregate(MUT[MUT$MutType == 'Syn' & MUT$Gene == 'CytB',]$Number, by = list(MUT[MUT$MutType == 'Syn' & MUT$Gene == 'CytB',]$Species), FUN = sum); names(AGG4) = c('Species','NumOfSynInCytB')
AGG = merge(AGG1,AGG2, all = TRUE); AGG = merge(AGG,AGG3, all = TRUE); AGG = merge(AGG,AGG4, all = TRUE);
AGG[is.na(AGG)] <- 0

# merge with Neutral A T G C
nrow(MUT) # 395157
MUT = merge(MUT,NeutralATGC[,c(1:6)], by = c("Species", "Gene"))
nrow(MUT) # 395269 ?????

MUT = merge(MUT,AGG, by = 'Species')
write.table(MUT, file = '../../Body/3Results/VertebratePolymorphisms.MutSpecData.txt', quote = FALSE)

