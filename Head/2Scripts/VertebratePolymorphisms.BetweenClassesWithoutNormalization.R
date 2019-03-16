###################################
###################################

rm(list=ls(all=TRUE))

pdf("../../Body/4Figures/VertebratePolymorphisms.BetweenClassesWithoutNormalization.R.01.pdf", width = 30, height = 20)

ColG = rgb(0.1,0.1,0.1,0.5)
ColT = rgb(0.1,0.1,1,0.5)
ColC = rgb(0.1,1,0.1,0.5)
ColA = rgb(1,0.1,0.1,0.5)

MUT = read.table("../../Body/3Results/Mutational_spectra_in_Chordata_ML.txt", header = TRUE)
length(unique(MUT$Species)) # 2404  SOME SPECIES HAVE THREE WORDS => CUT THE LAST AND MERGE WITH TAXA, OR EVEN CUT TWO LAST AND LEAVE JUST GENUS

####### associate species name with Class
### Taxa 1, Cut out the third world!!!!!!!!!!!!!!!!!
Taxa = read.table("../../Body/1Raw/TaxaFromKostya.Names.stat", sep = '\t',header = FALSE) 
Taxa$Species = gsub(";.*",'',Taxa$V1); 
for (i in (1:nrow(Taxa)))  {Taxa$Species[i] = paste(unlist(strsplit(as.character(Taxa$Species[i]),split = ' '))[1],unlist(strsplit(as.character(Taxa$Species[i]),split = ' '))[2], sep = '_')}
Taxa$Class = gsub(";Chordata;.*",'',Taxa$V1); Taxa$Class = gsub(".*;",'',Taxa$Class); table(Taxa$Class)
Taxa$Class = gsub('Actinopteri','Actinopterygii',Taxa$Class)
Taxa$Class = gsub("Testudines|Squamata|Crocodylia|Sphenodontia",'Reptilia',Taxa$Class)
length(unique(Taxa$Species)) # 1708
table(Taxa$Class)
Taxa = Taxa[,-1]

### Taxa 2, Cut out the third world!!!!!!!!!!!!!!!!!
TaxaMore = read.table("../../Body/1Raw/TaxaFromKostya.2NeedTaxa.tax.txt", sep = '\t',header = FALSE) 
TaxaMore$Species = ''
for (i in (1:nrow(TaxaMore)))  
{TaxaMore$Species[i] = paste(unlist(strsplit(as.character(TaxaMore$V1[i]),split = ' '))[1],unlist(strsplit(as.character(TaxaMore$V1[i]),split = ' '))[2], sep = '_')}
TaxaMore$Class = gsub("; Chordata;.*",'',TaxaMore$V2); 
TaxaMore$Class = gsub(".*; ",'',TaxaMore$Class); 
TaxaMore$Class = gsub('Actinopteri','Actinopterygii',TaxaMore$Class)
TaxaMore$Class = gsub("Testudines|Squamata|Crocodylia",'Reptilia',TaxaMore$Class)
table(TaxaMore$Class)
TaxaMore = TaxaMore[,-c(1,2)]

Taxa = rbind(Taxa,TaxaMore); Taxa = unique(Taxa)
table(TaxaMore$Class)

SpeciesWithoutTaxonomy = setdiff(unique(MUT$Species),unique(Taxa$Species)); length(SpeciesWithoutTaxonomy) # 63
write.table(SpeciesWithoutTaxonomy, "../../Body/2Derived/NeedTaxa.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)

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

##### FILTER 3: fourfold degenerate sites:
VecOfSynFourFoldDegenerateSites <- c('CTT', 'CTC', 'CTA', 'CTG', 
                                     'GTT', 'GTC', 'GTA', 'GTG', 
                                     'TCT', 'TCC', 'TCA', 'TCG', 
                                     'CCT', 'CCC', 'CCA', 'CCG', 
                                     'ACT', 'ACC', 'ACA', 'ACG', 
                                     'GCT', 'GCC', 'GCA', 'GCG', 
                                     'CGT', 'CGC', 'CGA', 'CGG', 
                                     'GGT', 'GGC', 'GGA', 'GGG')
length(unique(VecOfSynFourFoldDegenerateSites)) # 32
nrow(MUT) # 395157
MUT = MUT[MUT$AncestorCodon %in% VecOfSynFourFoldDegenerateSites & MUT$DescendantCodon %in% VecOfSynFourFoldDegenerateSites,]; nrow(MUT) # 209120

##### FILTER: Gene
nrow(MUT) # 209120
table(MUT$Gene)
# ATP6   ATP8   COX1   COX2   COX3   CytB    ND1    ND2    ND3    ND4   ND4L 
# 9752    625  13233   4264   4873 109049  14314  33809   2635  14511   2055  # there are quite many genes - do I need to analyse all of them focus on CYTB? Better - all and from time to time to check that results based on CYTB only are robust!!!
# fraction of CytB is 109049 / 209120  = 0.52 

MUT = merge(MUT,Taxa, all.x = TRUE)  ##### NOT ALL SPECIES HAVE TAXONOMY!!!!
nrow(MUT[MUT$Class == 'Mammalia',]) # 70053
table(MUT[MUT$Class == 'Mammalia',]$Gene)
# ATP6  ATP8  COX1  COX2  COX3  CytB   ND1   ND2   ND3   ND4  ND4L 
# 1592   184  4481  2044  1809 39112  3111  3624  1048  4934   897 
# fraction of CytB: 39112 / 70053 - 0.558 PAPER

##### derive observed number of mutations for each class (no normalization)

PieChartTable = c()
Equil = c()
par(mfrow = c(1,5))
VecOfClasses = c('Actinopterygii','Amphibia','Reptilia','Mammalia','Aves')
for (i in 1:length(VecOfClasses))
{ # i = 1
  Temp = MUT[MUT$Class == VecOfClasses[i],]
  NumberOfSpecies = length(unique(Temp$Species))
  title = paste(VecOfClasses[i],', N = ',NumberOfSpecies, sep = '')
  Temp$number = 1
  Agg = aggregate(Temp$number, by = list(Temp$Subs), FUN = sum); names(Agg) = c('Subs','Number')
  Agg$Number = Agg$Number/sum(Agg$Number)
  sum(Agg$Number) # 1 - 100%
  pie(Agg$Number, labels = Agg$Subs, main = title, col=rainbow(length(Agg$Subs)))
  
  Agg$Class = VecOfClasses[i];
  PieChartTable = rbind(PieChartTable,Agg)
  
  ToGFromG = sum(Agg[Agg$Subs %in% c('A_G','T_G','C_G'),]$Number) /  sum(Agg[Agg$Subs %in% c('G_A','G_T','G_C'),]$Number)
  ToTFromT = sum(Agg[Agg$Subs %in% c('A_T','G_T','C_T'),]$Number) /  sum(Agg[Agg$Subs %in% c('T_A','T_G','T_C'),]$Number)
  ToAFromA = sum(Agg[Agg$Subs %in% c('T_A','G_A','C_A'),]$Number) /  sum(Agg[Agg$Subs %in% c('A_T','A_G','A_C'),]$Number)
  ToCFromC =  sum(Agg[Agg$Subs %in% c('T_C','G_C','A_C'),]$Number) / sum(Agg[Agg$Subs %in% c('C_T','C_G','C_A'),]$Number)
  
  Equil = rbind(Equil,c(i,VecOfClasses[i],ToGFromG,ToTFromT,ToAFromA,ToCFromC))
}

PieChartTable = data.frame(PieChartTable)
write.table(PieChartTable,"../../Body/3Results/VertebratePolymorphisms.BetweenClassesWithoutNormalization.PieChartTable.txt")
Equil = data.frame(Equil); names(Equil) = c('i','Classes','G','T','C','A')

# i        Classes                 G                 T                C                 A
# 1 Actinopterygii 0.594741697416974  0.88695652173913 1.64319645994664  1.04326864396436
# 2       Amphibia   0.5125284738041  1.04096989966555 1.76627712854758 0.917777777777778
# 3       Reptilia 0.532318741450068 0.869721115537849 1.69920494699647   1.0426957687153
# 4       Mammalia 0.538739462063428 0.899331180811808 1.77489626556017  1.02578796561605
# 5           Aves 0.569868001820665 0.766223612197029  1.8578811369509  1.08779093507554

### equilibrium for each species:
Equilibrium = c()
MUT$number = 1
AGG = aggregate(MUT$number, by = list(MUT$Subs,MUT$Species,MUT$Class), FUN = sum); names(AGG) = c('Subs','Species','Class','Number')
VecOfSpecies = unique(MUT$Species)
for (i in 1:length(VecOfSpecies))
  { # i = 1
  Agg = AGG[AGG$Species == VecOfSpecies[i],]
  if (sum(Agg$Number) >= 15) ### do I need to add the same requirement for piechart, probably now, the first figure is expected to be completely general - all together
    {
    ToGFromG = log2(sum(Agg[Agg$Subs %in% c('A_G','T_G','C_G'),]$Number) /  sum(Agg[Agg$Subs %in% c('G_A','G_T','G_C'),]$Number))
    ToTFromT = log2(sum(Agg[Agg$Subs %in% c('A_T','G_T','C_T'),]$Number) /  sum(Agg[Agg$Subs %in% c('T_A','T_G','T_C'),]$Number))
    ToAFromA = log2(sum(Agg[Agg$Subs %in% c('T_A','G_A','C_A'),]$Number) /  sum(Agg[Agg$Subs %in% c('A_T','A_G','A_C'),]$Number))
    ToCFromC =  log2(sum(Agg[Agg$Subs %in% c('T_C','G_C','A_C'),]$Number) / sum(Agg[Agg$Subs %in% c('C_T','C_G','C_A'),]$Number))
    
    Equilibrium = rbind(Equilibrium,c(as.character(VecOfSpecies[i]),Agg$Class[1],ToGFromG,ToTFromT,ToAFromA,ToCFromC))
    }
  }
Equilibrium = data.frame(Equilibrium); names(Equilibrium)=c('Species','Class','G','T','A','C')
Equilibrium[,3:6] <- sapply(Equilibrium[,3:6], function(x) as.numeric(as.character(x)))
write.table(Equilibrium,"../../Body/3Results/VertebratePolymorphisms.BetweenClassesWithoutNormalization.EquilibriumLog2ToFrom.txt", quote = FALSE, row.names = FALSE)

par(mfrow=c(1,1))
Equilibrium = Equilibrium[Equilibrium$Class %in% VecOfClasses,]
boxplot(
Equilibrium[Equilibrium$Class == 'Actinopterygii',]$G,Equilibrium[Equilibrium$Class == 'Actinopterygii',]$T, Equilibrium[Equilibrium$Class == 'Actinopterygii',]$C, Equilibrium[Equilibrium$Class == 'Actinopterygii',]$A, 
Equilibrium[Equilibrium$Class == 'Amphibia',]$G,Equilibrium[Equilibrium$Class == 'Amphibia',]$T, Equilibrium[Equilibrium$Class == 'Amphibia',]$C, Equilibrium[Equilibrium$Class == 'Amphibia',]$A, 
Equilibrium[Equilibrium$Class == 'Reptilia',]$G,Equilibrium[Equilibrium$Class == 'Reptilia',]$T, Equilibrium[Equilibrium$Class == 'Reptilia',]$C, Equilibrium[Equilibrium$Class == 'Reptilia',]$A, 
Equilibrium[Equilibrium$Class == 'Mammalia',]$G,Equilibrium[Equilibrium$Class == 'Mammalia',]$T, Equilibrium[Equilibrium$Class == 'Mammalia',]$C, Equilibrium[Equilibrium$Class == 'Mammalia',]$A, 
Equilibrium[Equilibrium$Class == 'Aves',]$G,Equilibrium[Equilibrium$Class == 'Aves',]$T, Equilibrium[Equilibrium$Class == 'Aves',]$C, Equilibrium[Equilibrium$Class == 'Aves',]$A, 
outline = FALSE, notch = TRUE, col = c(ColG,ColT,ColC,ColA), names = rep(c('G','T','C','A'),5))
abline(h = 0, col = 'red')

dev.off()

############# DELETE BELOW: 
### little experiment: equilibrium and generation time - only log2(ToG/FromG) is increasing with GT -> short lived are more fard from equilibrium... strange...
GT = read.table("../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt", sep = '\t', header = TRUE)
GT$Species = gsub(' ','_',GT$Scientific_name)
GT = merge(Equilibrium,GT)
cor.test(GT$G,GT$GenerationLength_d, method = 'spearman') # positive a bit - small are running away faster...
cor.test(GT$C,GT$GenerationLength_d, method = 'spearman') # nothing
cor.test(GT$T,GT$GenerationLength_d, method = 'spearman') # nothing
cor.test(GT$A,GT$GenerationLength_d, method = 'spearman') # nothing


