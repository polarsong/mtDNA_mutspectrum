rm(list=ls(all=TRUE))

############ Syn mut
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
names(SynNuc)
MUT = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt', header = TRUE)

Codons = aggregate(list(CTT = SynNuc$CTT, CTC = SynNuc$CTC, CTA = SynNuc$CTA, CTG = SynNuc$CTG, GTT = SynNuc$GTT, GTC = SynNuc$GTC, GTA = SynNuc$GTA, GTG = SynNuc$GTG, TCT = SynNuc$TCT, TCC = SynNuc$TCC, TCA = SynNuc$TCA, TCG = SynNuc$TCG,
                        CCT = SynNuc$CCT, CCC = SynNuc$CCC, CCA = SynNuc$CCA, CCG = SynNuc$CCG, ACT = SynNuc$ACT, ACC = SynNuc$ACC, ACA = SynNuc$ACA, ACG = SynNuc$ACG, GCT = SynNuc$GCT, GCC = SynNuc$GCC, GCA = SynNuc$GCA,
                        GCG = SynNuc$GCG, CGT = SynNuc$CGT, CGC = SynNuc$CGC, CGA = SynNuc$CGA, CGG = SynNuc$CGG, GGT = SynNuc$GGT, GGC = SynNuc$GGC, GGA = SynNuc$GGA, GGG = SynNuc$GGG), by=list(Species = SynNuc$Species), FUN = sum)
CodonsMut = merge(Codons, MUT)
names(CodonsMut)

summary(lm(T_C ~ CTC + GTC + TCC + CCC + ACC + GCC + CGC + GGC, data = CodonsMut))
cor.test(CodonsMut$T_C, CodonsMut$CTC)
cor.test(CodonsMut$T_C, CodonsMut$GTC)
cor.test(CodonsMut$T_C, CodonsMut$TCC)
cor.test(CodonsMut$T_C, CodonsMut$CCC)
cor.test(CodonsMut$T_C, CodonsMut$ACC)
cor.test(CodonsMut$T_C, CodonsMut$GCC)
cor.test(CodonsMut$T_C, CodonsMut$CGC)
cor.test(CodonsMut$T_C, CodonsMut$GGC)
