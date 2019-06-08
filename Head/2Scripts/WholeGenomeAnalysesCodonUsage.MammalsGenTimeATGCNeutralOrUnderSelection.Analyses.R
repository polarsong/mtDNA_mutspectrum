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


##############################################################################################
############## PICs 

library(ape) # install.packages('ape') 

tree <- read.tree("../../Body/1Raw/mtalign.aln.treefile.rooted")

data <- AtgcGl[which(as.character(AtgcGl$species) %in% tree$tip.label),]
row.names(data) <- data$species
data$GenerationLength_d = log2(data$GenerationLength_d)

df_vec <- as.character(data$species)
tree_vec <- tree$tip.label

a <- setdiff(df_vec, tree_vec)
b <- setdiff(tree_vec, df_vec)

tree2 <- drop.tip(tree, b)

library(pacman)
p_load(tibble, dplyr, magrittr, purrr)
contrasts <- data %>% 
  select(-c(species, taxonomy)) %>% 
#   mutate_if(is.numeric, log2) %>% 
  map(pic, tree2)






cor.test(contrasts$GenerationLength_d, contrasts$WholeGenomeA, method = 'spearman')
cor.test(contrasts$GenerationLength_d, contrasts$UnderSelectionA, method = 'spearman')
cor.test(contrasts$GenerationLength_d, contrasts$NeutralA, method = 'spearman')
cor.test(contrasts$GenerationLength_d, contrasts$WholeGenomeAfr, method = 'spearman')
cor.test(contrasts$GenerationLength_d, contrasts$UnderSelectionAfr, method = 'spearman')
cor.test(contrasts$GenerationLength_d, contrasts$NeutralAfr, method = 'spearman')

cor.test(contrasts$GenerationLength_d, contrasts$WholeGenomeT, method = 'spearman')
cor.test(contrasts$GenerationLength_d, contrasts$UnderSelectionT, method = 'spearman')
cor.test(contrasts$GenerationLength_d, contrasts$NeutralT, method = 'spearman')
cor.test(contrasts$GenerationLength_d, contrasts$WholeGenomeTfr, method = 'spearman')
cor.test(contrasts$GenerationLength_d, contrasts$UnderSelectionTfr, method = 'spearman')
cor.test(contrasts$GenerationLength_d, contrasts$NeutralTfr, method = 'spearman')

cor.test(contrasts$GenerationLength_d, contrasts$WholeGenomeG, method = 'spearman')
cor.test(contrasts$GenerationLength_d, contrasts$UnderSelectionG, method = 'spearman')
cor.test(contrasts$GenerationLength_d, contrasts$NeutralG, method = 'spearman')
cor.test(contrasts$GenerationLength_d, contrasts$WholeGenomeGfr, method = 'spearman')
cor.test(contrasts$GenerationLength_d, contrasts$UnderSelectionGfr, method = 'spearman')
cor.test(contrasts$GenerationLength_d, contrasts$NeutralGfr, method = 'spearman')

cor.test(contrasts$GenerationLength_d, contrasts$WholeGenomeC, method = 'spearman')
cor.test(contrasts$GenerationLength_d, contrasts$UnderSelectionC, method = 'spearman')
cor.test(contrasts$GenerationLength_d, contrasts$NeutralC, method = 'spearman')
cor.test(contrasts$GenerationLength_d, contrasts$WholeGenomeCfr, method = 'spearman')
cor.test(contrasts$GenerationLength_d, contrasts$UnderSelectionCfr, method = 'spearman')
cor.test(contrasts$GenerationLength_d, contrasts$NeutralCfr, method = 'spearman')

########### 4 LM analyses: Neutral vs UnderSelection 

a <- lm((contrasts$GenerationLength_d) ~ contrasts$NeutralA + contrasts$UnderSelectionA); summary(a)
a <- lm((contrasts$GenerationLength_d) ~ contrasts$NeutralT + contrasts$UnderSelectionT); summary(a)
a <- lm((contrasts$GenerationLength_d) ~ contrasts$NeutralG + contrasts$UnderSelectionG); summary(a)
a <- lm((contrasts$GenerationLength_d) ~ contrasts$NeutralC + contrasts$UnderSelectionC); summary(a)

########### 5 which pairs correlate better with each other
cor.test(contrasts$NeutralA, contrasts$NeutralT, method = 'spearman')
cor.test(contrasts$NeutralA, contrasts$NeutralG, method = 'spearman')
cor.test(contrasts$NeutralA, contrasts$NeutralC, method = 'spearman') # -0.53
cor.test(contrasts$NeutralT, contrasts$NeutralG, method = 'spearman')
cor.test(contrasts$NeutralT, contrasts$NeutralC, method = 'spearman') # -0.52

########### 6 which pairs correlate better with each other
a <- lm((contrasts$GenerationLength_d) ~ contrasts$TotalNeutralA + contrasts$TotalNeutralT + contrasts$TotalNeutralG + contrasts$TotalNeutralC + contrasts$A + contrasts$T + contrasts$G + contrasts$C); summary(a)
a <- lm((contrasts$GenerationLength_d) ~ scale(contrasts$TotalNeutralA) + scale(contrasts$TotalNeutralT) + scale(contrasts$TotalNeutralG) + scale(contrasts$TotalNeutralC) + scale(contrasts$A) + scale(contrasts$T) + scale(contrasts$G) + scale(contrasts$C)); summary(a)
a <- lm((contrasts$GenerationLength_d) ~ scale(contrasts$TotalNeutralA) + scale(contrasts$TotalNeutralT) + scale(contrasts$TotalNeutralG) + scale(contrasts$A) + scale(contrasts$T) + scale(contrasts$G) + scale(contrasts$C)); summary(a)
a <- lm((contrasts$GenerationLength_d) ~ scale(contrasts$TotalNeutralA) + scale(contrasts$TotalNeutralG) + scale(contrasts$A) + scale(contrasts$T) + scale(contrasts$G) + scale(contrasts$C)); summary(a)

