################################
################## 
################################

rm(list=ls(all=TRUE))

############ Syn mut
Mut = read.table("../../Body/3Results/fulltreeCodons.csv", header = TRUE, sep = ';'); # "../../Body/3Results/fulltreeCodons.csv", sep=";"
table(Mut$note)
# gaps non-coding     normal 
# 787     261486     315003 
# Why gaps? they should be in coding or non-coding? strange names

Pc= Mut[Mut$note == 'normal',] # filter out everything except protein-coding mutations:
AncAa = data.frame(table(Pc$ancestral_aa)) #
AncAa = AncAa[order(-AncAa$Freq),]
AncAa
VectorOfAa=AncAa$Var1; VectorOfAa

table(Pc$derived_aa) # why  Ambiguous, Asn/Asp, Gln/Glu, Leu/Ile? => because of ambiguous derived AA. Whe there is no ambigios in ancestral? 
PcQuestion= Pc[Pc$derived_aa == 'Ambiguous',] 
head(PcQuestion,10)
PcQuestion= Pc[Pc$derived_aa == 'Asn/Asp' | Pc$derived_aa == 'Gln/Glu' | Pc$derived_aa == 'Leu/Ile',] 
head(PcQuestion,10)

## why problematic only in the descendant???? each branch has ancestor and descendant? 

## filter out all problems:
nrow(Pc)
PcGold = Pc[Pc$ancestral_aa %in% VectorOfAa & Pc$derived_aa %in% VectorOfAa,]
nrow(PcGold)

## syn / nons => too many nons 
table(PcGold$synonymous) 

## COX1, ND4, ND4L looks good; - others - not really 
ToStop = PcGold[PcGold$ancestral_aa != 'Stop' & PcGold$derived_aa == 'Stop',]; ToStop$synonymous = 'ToStop'
FromStop = PcGold[PcGold$ancestral_aa == 'Stop' & PcGold$derived_aa != 'Stop',]; FromStop$synonymous = 'FromStop'
NoStop = PcGold[PcGold$ancestral_aa != 'Stop' & PcGold$derived_aa != 'Stop',];
PcGold = rbind(NoStop,ToStop,FromStop)
T1=data.frame(table(PcGold$synonymous, by = PcGold$gene_info))
names(T1)=c('MutType','Gene','Freq')
T1=T1[grep('mRNA',T1$Gene),]; nrow(T1) # 60

## only ND4L and ND3 look good, 
Syn = T1[T1$MutType == 'synonymous',]; Syn$Syn = Syn$Freq; Syn = Syn[,grep("Gene|Syn", names(Syn))]
Nons = T1[T1$MutType == 'non-synonymous',]; Nons$Nons = Nons$Freq; Nons = Nons[,grep("Gene|Nons", names(Nons))]
ToStop = T1[T1$MutType == 'ToStop',]; ToStop$ToStop = ToStop$Freq; ToStop = ToStop[,grep("Gene|ToStop", names(ToStop))]
FromStop = T1[T1$MutType == 'FromStop',]; FromStop$FromStop = FromStop$Freq; FromStop = FromStop[,grep("Gene|FromStop", names(FromStop))]
MutTypes = merge(Syn,Nons, by = 'Gene', all = TRUE);
MutTypes = merge(MutTypes,ToStop, by = 'Gene', all = TRUE);
MutTypes = merge(MutTypes,FromStop, by = 'Gene', all = TRUE);
MutTypes



