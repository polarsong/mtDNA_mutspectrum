#statistics for birds with pheno, ratio and AnAge
rm(list=ls(all=TRUE))
aves_pheno = read.table('../../Body/1Raw/DataFromValya/ALL_PHENOTYPES.txt')
anage = read.table('../../Body/1Raw/anage_data.txt', header = TRUE, sep = '\t')
anage$Species = paste(anage$Genus, anage$Species, sep = '_')
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip") 
codon_usage = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t') 
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")) {file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")} 
codon_usage[codon_usage$Gene == 'ND6',] = NA
codon_usage = na.omit(codon_usage)
unique_names = unique(codon_usage$Species)
counter = 1
RatioA = c()
RatioG = c()
RatioC = c()
RatioT = c()
for (i in unique_names)
{
  df = codon_usage[codon_usage$Species == i,]
  Ad = (sum(df$NeutralA))/(sum(df$NeutralA)+sum(df$NeutralG)+sum(df$NeutralC)+sum(df$NeutralT))
  Gu = (sum(df$NeutralG))/(sum(df$NeutralA)+sum(df$NeutralG)+sum(df$NeutralC)+sum(df$NeutralT))
  Ci = (sum(df$NeutralC))/(sum(df$NeutralA)+sum(df$NeutralG)+sum(df$NeutralC)+sum(df$NeutralT))
  Ti = (sum(df$NeutralT))/(sum(df$NeutralA)+sum(df$NeutralG)+sum(df$NeutralC)+sum(df$NeutralT))
  RatioA[counter] = Ad
  RatioG[counter] = Gu
  RatioC[counter] = Ci
  RatioT[counter] = Ti
  counter = counter + 1
}
all_ratio = data.frame(unique_names)
names(all_ratio) = 'Species'
all_ratio$RatioA = RatioA
all_ratio$RatioC = RatioC
all_ratio$RatioG = RatioG
all_ratio$RatioT = RatioT
names(aves_pheno) = c('Species', 'Phenotype')
aves_ratio_and_pheno = merge(aves_pheno, all_ratio, by = 'Species' )
aves_ratio_and_pheno_and_anage = merge(aves_ratio_and_pheno, anage, by = 'Species')

cor.test(aves_ratio_and_pheno_and_anage$RatioG, aves_ratio_and_pheno_and_anage$Maximum.longevity..yrs.., method = 'spearman')
#cannot be counted, because some birds have several phenotypes




#statistics for birds with pheno, mutspec and AnAge
mutspec = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt', header = TRUE)
birds_mutspec_and_pheno = merge(mutspec, aves_pheno, by = 'Species')
birds_mutspec_and_pheno_and_anage = merge(birds_mutspec_and_pheno, anage, by = 'Species')
cor.test(birds_mutspec_and_pheno_and_anage$G_A, birds_mutspec_and_pheno_and_anage$Maximum.longevity..yrs., method = 'spearman')
cor.test(birds_mutspec_and_pheno_and_anage$T_C, birds_mutspec_and_pheno_and_anage$Maximum.longevity..yrs., method = 'spearman')
#i can do more analysis, but what parameters are better to choose?



#statistics for birds with ratio, mutspec and anage
aves_anage = anage[anage$Class == 'Aves',]
aves_anage_and_ratio = merge(aves_anage, all_ratio, by = 'Species')
aves_anage_and_ratio_and_mutspec = merge(aves_anage_and_ratio, mutspec, by = 'Species')
cor.test(aves_anage_and_ratio_and_mutspec$G_A, aves_anage_and_ratio_and_mutspec$Maximum.longevity..yrs., method = 'spearman')
cor.test(aves_anage_and_ratio_and_mutspec$T_C, aves_anage_and_ratio_and_mutspec$Maximum.longevity..yrs., method = 'spearman')
cor.test(aves_anage_and_ratio_and_mutspec$RatioA, aves_anage_and_ratio_and_mutspec$Maximum.longevity..yrs., method = 'spearman')
cor.test(aves_anage_and_ratio_and_mutspec$RatioC, aves_anage_and_ratio_and_mutspec$Maximum.longevity..yrs., method = 'spearman')
cor.test(aves_anage_and_ratio_and_mutspec$RatioG, aves_anage_and_ratio_and_mutspec$Maximum.longevity..yrs., method = 'spearman')
cor.test(aves_anage_and_ratio_and_mutspec$RatioT, aves_anage_and_ratio_and_mutspec$Maximum.longevity..yrs., method = 'spearman')


cor.test(aves_anage_and_ratio_and_mutspec$G_A, aves_anage_and_ratio_and_mutspec$Metabolic.rate..W., method = 'spearman')
cor.test(aves_anage_and_ratio_and_mutspec$T_C, aves_anage_and_ratio_and_mutspec$Metabolic.rate..W., method = 'spearman')
cor.test(aves_anage_and_ratio_and_mutspec$RatioA, aves_anage_and_ratio_and_mutspec$Metabolic.rate..W., method = 'spearman')
cor.test(aves_anage_and_ratio_and_mutspec$RatioC, aves_anage_and_ratio_and_mutspec$Metabolic.rate..W., method = 'spearman')
cor.test(aves_anage_and_ratio_and_mutspec$RatioG, aves_anage_and_ratio_and_mutspec$Metabolic.rate..W., method = 'spearman')
cor.test(aves_anage_and_ratio_and_mutspec$RatioT, aves_anage_and_ratio_and_mutspec$Metabolic.rate..W., method = 'spearman')


cor.test(aves_anage_and_ratio_and_mutspec$G_A, aves_anage_and_ratio_and_mutspec$Body.mass..g., method = 'spearman')
cor.test(aves_anage_and_ratio_and_mutspec$T_C, aves_anage_and_ratio_and_mutspec$Body.mass..g., method = 'spearman')
cor.test(aves_anage_and_ratio_and_mutspec$RatioA, aves_anage_and_ratio_and_mutspec$Body.mass..g., method = 'spearman')
cor.test(aves_anage_and_ratio_and_mutspec$RatioC, aves_anage_and_ratio_and_mutspec$Body.mass..g., method = 'spearman')
cor.test(aves_anage_and_ratio_and_mutspec$RatioG, aves_anage_and_ratio_and_mutspec$Body.mass..g., method = 'spearman')
cor.test(aves_anage_and_ratio_and_mutspec$RatioT, aves_anage_and_ratio_and_mutspec$Body.mass..g., method = 'spearman')
#i can do more analysis, but what parameters are better to choose?
