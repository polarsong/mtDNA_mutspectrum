rm(list=ls(all=TRUE))
library('ggplot2')
#merging some phenotypes and life history traits to find some good correlations
#reading phenotypes
pheno = read.table('../../Body/1Raw/DataFromValya/ALL_PHENOTYPES.txt')
#reading anage
anage = read.table('../../Body/1Raw/anage_data.txt',header = TRUE, sep = '\t')
names(pheno) = c('Species', 'Phenotype')
#getting only tropical birds
only_tropical = pheno[pheno$Phenotype == ",\"TROPIC\"",]
#getting only high_altitude birds
only_altitude = pheno[pheno$Phenotype == ",\"HI\"",]
#reading codon usage
unzip("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt.zip") 
codus = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t') 
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")) {file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")} 
codus = codus[codus$Gene != 'ND6',]
#### aggregate (summ up) neutral nucleotides for each species
SynNuc = aggregate(list(codus$NeutralA,codus$NeutralT,codus$NeutralG,codus$NeutralC), by = list(codus$Species,codus$Class,codus$Taxonomy), FUN = sum)
names(SynNuc)=c('Species','Class','Taxonomy','NeutralA','NeutralT','NeutralG','NeutralC')
table(SynNuc$Class)

#### estimate fraction of each nucleotide 
SynNuc$FrA = SynNuc$NeutralA / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 
SynNuc$FrT = SynNuc$NeutralT / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC)
SynNuc$FrG = SynNuc$NeutralG / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC)
SynNuc$FrC = SynNuc$NeutralC / (SynNuc$NeutralA + SynNuc$NeutralT + SynNuc$NeutralG + SynNuc$NeutralC) 
#merging high altitude birds and tropical birds with SynNuc
tropical_and_fr = merge(SynNuc, only_tropical, by = 'Species')
altitude_and_fr = merge(SynNuc, only_altitude, by = 'Species')
random_rows = sample(1:81, 18, replace = FALSE)
random_birds_from_tropical = data.frame()
for (i in random_rows)
{
  random_birds_from_tropical = rbind(random_birds_from_tropical, tropical_and_fr[i,])
}
#graphics
boxplot(altitude_and_fr$FrA, random_birds_from_tropical$FrA)
boxplot(altitude_and_fr$FrG, random_birds_from_tropical$FrG)
boxplot(altitude_and_fr$FrC, random_birds_from_tropical$FrC)
boxplot(altitude_and_fr$FrT, random_birds_from_tropical$FrT)

boxplot(altitude_and_fr$FrA, tropical_and_fr$FrA)
boxplot(altitude_and_fr$FrG, tropical_and_fr$FrG)
boxplot(altitude_and_fr$FrC, tropical_and_fr$FrC)
boxplot(altitude_and_fr$FrT, tropical_and_fr$FrT)


#getting some extra groups
pheno1 = pheno[pheno$Phenotype != ",\"HI\"",]
pheno1 = pheno1[pheno1$Phenotype != ",\"WI\"",]
pheno1 = pheno1[pheno1$Phenotype != ",\"AI\"",]
pheno1 = pheno1[pheno1$Phenotype != ",\"DI\"",]
pheno1 = pheno1[pheno1$Phenotype != ",\"FM\"",]
nf_and_tropic_ratio = merge(pheno1, SynNuc, by = 'Species')
pheno2 = pheno[pheno$Phenotype != ",\"DI\"",]
pheno2 = pheno2[pheno2$Phenotype != ",\"FM\"",]
pheno2 = pheno2[pheno2$Phenotype != ",\"TROPIC\"",]
pheno2 = pheno2[pheno2$Phenotype != ",\"NF\"",]
pheno2 = pheno2[pheno2$Phenotype != ",\"AI\"",]
hi_and_wintering_ratio = merge(pheno2, SynNuc, by = 'Species')
#boxploting
boxplot(nf_and_tropic_ratio$FrA, hi_and_wintering_ratio$FrA)
boxplot(nf_and_tropic_ratio$FrG, hi_and_wintering_ratio$FrG)
boxplot(nf_and_tropic_ratio$FrC, hi_and_wintering_ratio$FrC)
boxplot(nf_and_tropic_ratio$FrT, hi_and_wintering_ratio$FrT)


#looking at birds from lipid paper
first_species = SynNuc[SynNuc$Species == 'Tachycineta_bicolor',]
second_species = SynNuc[SynNuc$Species == 'Tachycineta_albilinea',]
boxplot(first_species$FrA, second_species$FrA)
boxplot(first_species$FrG, second_species$FrG)
boxplot(first_species$FrC, second_species$FrC)
boxplot(first_species$FrT, second_species$FrT)

unique(pheno$Phenotype)
#lots of drawing
pheno_and_SynNuc = merge(SynNuc, pheno, by = 'Species')
pdf("../../Body/3Results/birds boxplots no notch.pdf")
par(mfrow)=c(1,4) #didn't work
boxplot(pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"TROPIC\"",]$FrA, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"WI\"",]$FrA, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"HI\"",]$FrA,
        pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"FM\"",]$FrA, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"NF\"",]$FrA, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"DI\"",]$FrA, 
        pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"AI\"",]$FrA, names=c('Trop94','Wint13','HiAlt18', 'LD64', 'NF20', 'DI10', 'AI45'), ylab = 'FrA')
boxplot(pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"TROPIC\"",]$FrG, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"WI\"",]$FrG, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"HI\"",]$FrG,
        pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"FM\"",]$FrG, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"NF\"",]$FrG, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"DI\"",]$FrG, 
        pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"AI\"",]$FrG, names=c('Trop94','Wint13','HiAlt18', 'LD64', 'NF20', 'DI10', 'AI45'), ylab = 'FrG')
boxplot(pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"TROPIC\"",]$FrC, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"WI\"",]$FrC, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"HI\"",]$FrC,
        pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"FM\"",]$FrC, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"NF\"",]$FrC, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"DI\"",]$FrC, 
        pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"AI\"",]$FrC, names=c('Trop94','Wint13','HiAlt18', 'LD64', 'NF20', 'DI10', 'AI45'), ylab = 'FrC')
boxplot(pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"TROPIC\"",]$FrT, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"WI\"",]$FrT, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"HI\"",]$FrT,
        pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"FM\"",]$FrT, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"NF\"",]$FrT, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"DI\"",]$FrT, 
        pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"AI\"",]$FrT, names=c('Trop94','Wint13','HiAlt18', 'LD64', 'NF20', 'DI10', 'AI45'), ylab = 'FrT')
dev.off()





pdf("../../Body/3Results/birds boxplots with notch.pdf")
par(mfrow)=c(1,4) #didn't work
boxplot(pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"TROPIC\"",]$FrA, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"WI\"",]$FrA, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"HI\"",]$FrA,
        pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"FM\"",]$FrA, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"NF\"",]$FrA, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"DI\"",]$FrA, 
        pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"AI\"",]$FrA, names=c('Trop94','Wint13','HiAlt18', 'LD64', 'NF20', 'DI10', 'AI45'), ylab = 'FrA',notch = TRUE, outline = FALSE)
boxplot(pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"TROPIC\"",]$FrG, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"WI\"",]$FrG, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"HI\"",]$FrG,
        pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"FM\"",]$FrG, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"NF\"",]$FrG, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"DI\"",]$FrG, 
        pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"AI\"",]$FrG, names=c('Trop94','Wint13','HiAlt18', 'LD64', 'NF20', 'DI10', 'AI45'), ylab = 'FrG',notch = TRUE, outline = FALSE)
boxplot(pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"TROPIC\"",]$FrC, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"WI\"",]$FrC, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"HI\"",]$FrC,
        pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"FM\"",]$FrC, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"NF\"",]$FrC, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"DI\"",]$FrC, 
        pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"AI\"",]$FrC, names=c('Trop94','Wint13','HiAlt18', 'LD64', 'NF20', 'DI10', 'AI45'), ylab = 'FrC',notch = TRUE, outline = FALSE)
boxplot(pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"TROPIC\"",]$FrT, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"WI\"",]$FrT, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"HI\"",]$FrT,
        pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"FM\"",]$FrT, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"NF\"",]$FrT, pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"DI\"",]$FrT, 
        pheno_and_SynNuc[pheno_and_SynNuc$Phenotype == ",\"AI\"",]$FrT, names=c('Trop94','Wint13','HiAlt18', 'LD64', 'NF20', 'DI10', 'AI45'), ylab = 'FrT',notch = TRUE, outline = FALSE)
dev.off()


#working with mutspec
#transitions
MutSpec = read.table("../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt", header = TRUE) 
mutspec_and_pheno = merge(MutSpec, pheno, by = 'Species')
pdf("../../Body/3Results/birds boxplots mutspec no notch correct.pdf")
boxplot(mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"TROPIC\"",]$A_G, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"WI\"",]$A_G,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"HI\"",]$A_G, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"FM\"",]$A_G, 
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"NF\"",]$A_G, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"DI\"",]$A_G,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"AI\"",]$A_G, names=c('Trop14','Wint7','HiAlt3', 'LD10', 'NF1', 'DI1', 'AI5'), ylab = 'A->G')
boxplot(mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"TROPIC\"",]$G_A, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"WI\"",]$G_A,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"HI\"",]$G_A, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"FM\"",]$G_A, 
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"NF\"",]$G_A, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"DI\"",]$G_A,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"AI\"",]$G_A, names=c('Trop14','Wint7','HiAlt3', 'LD10', 'NF1', 'DI1', 'AI5'), ylab = 'G->A')
boxplot(mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"TROPIC\"",]$T_C, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"WI\"",]$T_C,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"HI\"",]$T_C, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"FM\"",]$T_C, 
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"NF\"",]$T_C, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"DI\"",]$T_C,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"AI\"",]$T_C, names=c('Trop14','Wint7','HiAlt3', 'LD10', 'NF1', 'DI1', 'AI5'), ylab = 'T->C')
boxplot(mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"TROPIC\"",]$C_T, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"WI\"",]$C_T,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"HI\"",]$C_T, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"FM\"",]$C_T, 
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"NF\"",]$C_T, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"DI\"",]$C_T,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"AI\"",]$C_T, names=c('Trop14','Wint7','HiAlt3', 'LD10', 'NF1', 'DI1', 'AI5'), ylab = 'C->T')
dev.off()

pdf("../../Body/3Results/birds boxplots mutspec with notch correct.pdf")
boxplot(mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"TROPIC\"",]$A_G, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"WI\"",]$A_G,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"HI\"",]$A_G, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"FM\"",]$A_G, 
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"NF\"",]$A_G, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"DI\"",]$A_G,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"AI\"",]$A_G, names=c('Trop14','Wint7','HiAlt3', 'LD10', 'NF1', 'DI1', 'AI5'), ylab = 'A->G', notch = TRUE, outline = FALSE)
boxplot(mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"TROPIC\"",]$G_A, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"WI\"",]$G_A,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"HI\"",]$G_A, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"FM\"",]$G_A, 
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"NF\"",]$G_A, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"DI\"",]$G_A,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"AI\"",]$G_A, names=c('Trop14','Wint7','HiAlt3', 'LD10', 'NF1', 'DI1', 'AI5'), ylab = 'G->A', notch = TRUE, outline = FALSE)
boxplot(mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"TROPIC\"",]$T_C, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"WI\"",]$T_C,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"HI\"",]$T_C, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"FM\"",]$T_C, 
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"NF\"",]$T_C, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"DI\"",]$T_C,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"AI\"",]$T_C, names=c('Trop14','Wint7','HiAlt3', 'LD10', 'NF1', 'DI1', 'AI5'), ylab = 'T->C', notch = TRUE, outline = FALSE)
boxplot(mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"TROPIC\"",]$C_T, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"WI\"",]$C_T,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"HI\"",]$C_T, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"FM\"",]$C_T, 
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"NF\"",]$C_T, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"DI\"",]$C_T,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"AI\"",]$C_T, names=c('Trop14','Wint7','HiAlt3', 'LD10', 'NF1', 'DI1', 'AI5'), ylab = 'C->T', notch = TRUE, outline = FALSE)
dev.off()

#transvertions
pdf("../../Body/3Results/birds boxplots mutspec transvertions no notch correct.pdf")
boxplot(mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"TROPIC\"",]$A_T, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"WI\"",]$A_T,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"HI\"",]$A_T, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"FM\"",]$A_T, 
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"NF\"",]$A_T, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"DI\"",]$A_T,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"AI\"",]$A_T, names=c('Trop14','Wint7','HiAlt3', 'LD10', 'NF1', 'DI1', 'AI5'), ylab = 'A->T')
boxplot(mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"TROPIC\"",]$A_C, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"WI\"",]$A_C,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"HI\"",]$A_C, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"FM\"",]$A_C, 
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"NF\"",]$A_C, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"DI\"",]$A_C,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"AI\"",]$A_C, names=c('Trop14','Wint7','HiAlt3', 'LD10', 'NF1', 'DI1', 'AI5'), ylab = 'A->C')
boxplot(mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"TROPIC\"",]$G_C, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"WI\"",]$G_C,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"HI\"",]$G_C, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"FM\"",]$G_C, 
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"NF\"",]$G_C, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"DI\"",]$G_C,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"AI\"",]$G_C, names=c('Trop14','Wint7','HiAlt3', 'LD10', 'NF1', 'DI1', 'AI5'), ylab = 'G->C')
boxplot(mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"TROPIC\"",]$G_T, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"WI\"",]$G_T,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"HI\"",]$G_T, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"FM\"",]$G_T, 
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"NF\"",]$G_T, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"DI\"",]$G_T,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"AI\"",]$G_T, names=c('Trop14','Wint7','HiAlt3', 'LD10', 'NF1', 'DI1', 'AI5'), ylab = 'G->T')
boxplot(mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"TROPIC\"",]$T_G, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"WI\"",]$T_G,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"HI\"",]$T_G, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"FM\"",]$T_G, 
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"NF\"",]$T_G, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"DI\"",]$T_G,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"AI\"",]$T_G, names=c('Trop14','Wint7','HiAlt3', 'LD10', 'NF1', 'DI1', 'AI5'), ylab = 'T->G')
boxplot(mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"TROPIC\"",]$T_A, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"WI\"",]$T_A,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"HI\"",]$T_A, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"FM\"",]$T_A, 
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"NF\"",]$T_A, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"DI\"",]$T_A,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"AI\"",]$T_A, names=c('Trop14','Wint7','HiAlt3', 'LD10', 'NF1', 'DI1', 'AI5'), ylab = 'T->A')
boxplot(mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"TROPIC\"",]$C_A, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"WI\"",]$C_A,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"HI\"",]$C_A, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"FM\"",]$C_A, 
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"NF\"",]$C_A, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"DI\"",]$C_A,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"AI\"",]$C_A, names=c('Trop14','Wint7','HiAlt3', 'LD10', 'NF1', 'DI1', 'AI5'), ylab = 'C->A')
boxplot(mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"TROPIC\"",]$C_G, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"WI\"",]$C_G,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"HI\"",]$C_G, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"FM\"",]$C_G, 
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"NF\"",]$C_G, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"DI\"",]$C_G,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"AI\"",]$C_G, names=c('Trop14','Wint7','HiAlt3', 'LD10', 'NF1', 'DI1', 'AI5'), ylab = 'C->G')
dev.off()

pdf("../../Body/3Results/birds boxplots mutspec transvertions with notch correct.pdf")
boxplot(mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"TROPIC\"",]$A_T, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"WI\"",]$A_T,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"HI\"",]$A_T, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"FM\"",]$A_T, 
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"NF\"",]$A_T, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"DI\"",]$A_T,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"AI\"",]$A_T, names=c('Trop14','Wint7','HiAlt3', 'LD10', 'NF1', 'DI1', 'AI5'), ylab = 'A->T', notch = TRUE, outline = FALSE)
boxplot(mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"TROPIC\"",]$A_C, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"WI\"",]$A_C,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"HI\"",]$A_C, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"FM\"",]$A_C, 
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"NF\"",]$A_C, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"DI\"",]$A_C,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"AI\"",]$A_C, names=c('Trop14','Wint7','HiAlt3', 'LD10', 'NF1', 'DI1', 'AI5'), ylab = 'A->C', notch = TRUE, outline = FALSE)
boxplot(mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"TROPIC\"",]$G_C, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"WI\"",]$G_C,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"HI\"",]$G_C, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"FM\"",]$G_C, 
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"NF\"",]$G_C, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"DI\"",]$G_C,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"AI\"",]$G_C, names=c('Trop14','Wint7','HiAlt3', 'LD10', 'NF1', 'DI1', 'AI5'), ylab = 'G->C', notch = TRUE, outline = FALSE)
boxplot(mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"TROPIC\"",]$G_T, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"WI\"",]$G_T,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"HI\"",]$G_T, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"FM\"",]$G_T, 
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"NF\"",]$G_T, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"DI\"",]$G_T,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"AI\"",]$G_T, names=c('Trop14','Wint7','HiAlt3', 'LD10', 'NF1', 'DI1', 'AI5'), ylab = 'G->T', notch = TRUE, outline = FALSE)
boxplot(mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"TROPIC\"",]$T_G, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"WI\"",]$T_G,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"HI\"",]$T_G, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"FM\"",]$T_G, 
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"NF\"",]$T_G, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"DI\"",]$T_G,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"AI\"",]$T_G, names=c('Trop14','Wint7','HiAlt3', 'LD10', 'NF1', 'DI1', 'AI5'), ylab = 'T->G', notch = TRUE, outline = FALSE)
boxplot(mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"TROPIC\"",]$T_A, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"WI\"",]$T_A,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"HI\"",]$T_A, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"FM\"",]$T_A, 
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"NF\"",]$T_A, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"DI\"",]$T_A,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"AI\"",]$T_A, names=c('Trop14','Wint7','HiAlt3', 'LD10', 'NF1', 'DI1', 'AI5'), ylab = 'T->A', notch = TRUE, outline = FALSE)
boxplot(mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"TROPIC\"",]$C_A, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"WI\"",]$C_A,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"HI\"",]$C_A, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"FM\"",]$C_A, 
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"NF\"",]$C_A, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"DI\"",]$C_A,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"AI\"",]$C_A, names=c('Trop14','Wint7','HiAlt3', 'LD10', 'NF1', 'DI1', 'AI5'), ylab = 'C->A', notch = TRUE, outline = FALSE)
boxplot(mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"TROPIC\"",]$C_G, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"WI\"",]$C_G,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"HI\"",]$C_G, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"FM\"",]$C_G, 
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"NF\"",]$C_G, mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"DI\"",]$C_G,
        mutspec_and_pheno[mutspec_and_pheno$Phenotype == ",\"AI\"",]$C_G, names=c('Trop14','Wint7','HiAlt3', 'LD10', 'NF1', 'DI1', 'AI5'), ylab = 'C->G', notch = TRUE, outline = FALSE)
dev.off()

