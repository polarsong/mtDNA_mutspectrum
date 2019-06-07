rm(list=ls(all=TRUE))

# PURPOSE: to understand which nucleotide content correlate better with Generation Length:  WholeGenome (A T G C) or UnderSelection (A T G C) or Neutral (A T G C)

########### 1 read a table with 649 unique mammalian species with A T G C per whole mit genome, A T G C in neutral sites (with removed overlaps between genes) and generation length.

AtgcGl = read.table("../../Body/3Results/WholeGenomeAnalysesCodonUsage.MammalsGenTimeATGCNeutralOrUnderSelection.R.txt")

########### 2 derive fractions (12 types):
AtgcGl$NeutralAfr = AtgcGl$NeutralA / (AtgcGl$NeutralA + AtgcGl$NeutralT + AtgcGl$NeutralG + AtgcGl$NeutralC)
AtgcGl$NeutralTfr = AtgcGl$NeutralT / (AtgcGl$NeutralA + AtgcGl$NeutralT + AtgcGl$NeutralG + AtgcGl$NeutralC)
AtgcGl$NeutralGfr = AtgcGl$NeutralG / (AtgcGl$NeutralA + AtgcGl$NeutralT + AtgcGl$NeutralG + AtgcGl$NeutralC)
AtgcGl$NeutralCfr = AtgcGl$NeutralC / (AtgcGl$NeutralA + AtgcGl$NeutralT + AtgcGl$NeutralG + AtgcGl$NeutralC)

AtgcGl$UnderSelectionAfr = AtgcGl$UnderSelectionA / (AtgcGl$UnderSelectionA + AtgcGl$UnderSelectionT + AtgcGl$UnderSelectionG + AtgcGl$UnderSelectionC)
AtgcGl$UnderSelectionTfr = AtgcGl$UnderSelectionT / (AtgcGl$UnderSelectionA + AtgcGl$UnderSelectionT + AtgcGl$UnderSelectionG + AtgcGl$UnderSelectionC)
AtgcGl$UnderSelectionGfr = AtgcGl$UnderSelectionG / (AtgcGl$UnderSelectionA + AtgcGl$UnderSelectionT + AtgcGl$UnderSelectionG + AtgcGl$UnderSelectionC)
AtgcGl$UnderSelectionCfr = AtgcGl$UnderSelectionC / (AtgcGl$UnderSelectionA + AtgcGl$UnderSelectionT + AtgcGl$UnderSelectionG + AtgcGl$UnderSelectionC)

AtgcGl$WholeGenomeAfr = AtgcGl$WholeGenomeA / (AtgcGl$WholeGenomeA + AtgcGl$WholeGenomeT + AtgcGl$WholeGenomeG + AtgcGl$WholeGenomeC)
AtgcGl$WholeGenomeTfr = AtgcGl$WholeGenomeT / (AtgcGl$WholeGenomeA + AtgcGl$WholeGenomeT + AtgcGl$WholeGenomeG + AtgcGl$WholeGenomeC)
AtgcGl$WholeGenomeGfr = AtgcGl$WholeGenomeG / (AtgcGl$WholeGenomeA + AtgcGl$WholeGenomeT + AtgcGl$WholeGenomeG + AtgcGl$WholeGenomeC)
AtgcGl$WholeGenomeCfr = AtgcGl$WholeGenomeC / (AtgcGl$WholeGenomeA + AtgcGl$WholeGenomeT + AtgcGl$WholeGenomeG + AtgcGl$WholeGenomeC)

########### 3 pairwise cor analyses 
# A T are decreasing, G C are increasing
# UnderSelectionT is decreasing with GL faster than NeutralT! why!!!??? PICs!?

cor.test(AtgcGl$GenerationLength_d, AtgcGl$WholeGenomeA, method = 'spearman')
cor.test(AtgcGl$GenerationLength_d, AtgcGl$UnderSelectionA, method = 'spearman')
cor.test(AtgcGl$GenerationLength_d, AtgcGl$NeutralA, method = 'spearman')
cor.test(AtgcGl$GenerationLength_d, AtgcGl$WholeGenomeAfr, method = 'spearman')
cor.test(AtgcGl$GenerationLength_d, AtgcGl$UnderSelectionAfr, method = 'spearman')
cor.test(AtgcGl$GenerationLength_d, AtgcGl$NeutralAfr, method = 'spearman')

cor.test(AtgcGl$GenerationLength_d, AtgcGl$WholeGenomeT, method = 'spearman')
cor.test(AtgcGl$GenerationLength_d, AtgcGl$UnderSelectionT, method = 'spearman')
cor.test(AtgcGl$GenerationLength_d, AtgcGl$NeutralT, method = 'spearman')
cor.test(AtgcGl$GenerationLength_d, AtgcGl$WholeGenomeTfr, method = 'spearman')
cor.test(AtgcGl$GenerationLength_d, AtgcGl$UnderSelectionTfr, method = 'spearman')
cor.test(AtgcGl$GenerationLength_d, AtgcGl$NeutralTfr, method = 'spearman')

cor.test(AtgcGl$GenerationLength_d, AtgcGl$WholeGenomeG, method = 'spearman')
cor.test(AtgcGl$GenerationLength_d, AtgcGl$UnderSelectionG, method = 'spearman')
cor.test(AtgcGl$GenerationLength_d, AtgcGl$NeutralG, method = 'spearman')
cor.test(AtgcGl$GenerationLength_d, AtgcGl$WholeGenomeGfr, method = 'spearman')
cor.test(AtgcGl$GenerationLength_d, AtgcGl$UnderSelectionGfr, method = 'spearman')
cor.test(AtgcGl$GenerationLength_d, AtgcGl$NeutralGfr, method = 'spearman')

cor.test(AtgcGl$GenerationLength_d, AtgcGl$WholeGenomeC, method = 'spearman')
cor.test(AtgcGl$GenerationLength_d, AtgcGl$UnderSelectionC, method = 'spearman')
cor.test(AtgcGl$GenerationLength_d, AtgcGl$NeutralC, method = 'spearman')
cor.test(AtgcGl$GenerationLength_d, AtgcGl$WholeGenomeCfr, method = 'spearman')
cor.test(AtgcGl$GenerationLength_d, AtgcGl$UnderSelectionCfr, method = 'spearman')
cor.test(AtgcGl$GenerationLength_d, AtgcGl$NeutralCfr, method = 'spearman')

########### 4 LM analyses: Neutral vs UnderSelection 

a <- lm(log2(AtgcGl$GenerationLength_d) ~ AtgcGl$NeutralA + AtgcGl$UnderSelectionA); summary(a)
a <- lm(log2(AtgcGl$GenerationLength_d) ~ AtgcGl$NeutralT + AtgcGl$UnderSelectionT); summary(a)
a <- lm(log2(AtgcGl$GenerationLength_d) ~ AtgcGl$NeutralG + AtgcGl$UnderSelectionG); summary(a)
a <- lm(log2(AtgcGl$GenerationLength_d) ~ AtgcGl$NeutralC + AtgcGl$UnderSelectionC); summary(a)

########### 5 which pairs correlate better with each other
cor.test(AtgcGl$NeutralA, AtgcGl$NeutralT, method = 'spearman')
cor.test(AtgcGl$NeutralA, AtgcGl$NeutralG, method = 'spearman')
cor.test(AtgcGl$NeutralA, AtgcGl$NeutralC, method = 'spearman') # -0.53
cor.test(AtgcGl$NeutralT, AtgcGl$NeutralG, method = 'spearman')
cor.test(AtgcGl$NeutralT, AtgcGl$NeutralC, method = 'spearman') # -0.52

########### 6 which pairs correlate better with each other
a <- lm(log2(AtgcGl$GenerationLength_d) ~ AtgcGl$TotalNeutralA + AtgcGl$TotalNeutralT + AtgcGl$TotalNeutralG + AtgcGl$TotalNeutralC + AtgcGl$A + AtgcGl$T + AtgcGl$G + AtgcGl$C); summary(a)
a <- lm(log2(AtgcGl$GenerationLength_d) ~ scale(AtgcGl$TotalNeutralA) + scale(AtgcGl$TotalNeutralT) + scale(AtgcGl$TotalNeutralG) + scale(AtgcGl$TotalNeutralC) + scale(AtgcGl$A) + scale(AtgcGl$T) + scale(AtgcGl$G) + scale(AtgcGl$C)); summary(a)
a <- lm(log2(AtgcGl$GenerationLength_d) ~ scale(AtgcGl$TotalNeutralA) + scale(AtgcGl$TotalNeutralT) + scale(AtgcGl$TotalNeutralG) + scale(AtgcGl$A) + scale(AtgcGl$T) + scale(AtgcGl$G) + scale(AtgcGl$C)); summary(a)
a <- lm(log2(AtgcGl$GenerationLength_d) ~ scale(AtgcGl$TotalNeutralA) + scale(AtgcGl$TotalNeutralG) + scale(AtgcGl$A) + scale(AtgcGl$T) + scale(AtgcGl$G) + scale(AtgcGl$C)); summary(a)

