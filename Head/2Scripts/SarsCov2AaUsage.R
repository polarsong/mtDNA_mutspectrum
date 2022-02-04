
rm(list=ls(all=TRUE))

AA = read.table("../../Body/1Raw/orf1ab_aa_usage.csv", sep = ';', header = TRUE)
dim(AA)
for (i in 1:nrow(AA))
{
  AA$Total[i] = sum(AA[i,4:23])  
}

names(AA)
AA$AlaFr = AA$Ala/AA$Total
AA$ArgFr = AA$Arg/AA$Total
AA$AsnFr = AA$Asn/AA$Total
AA$AspFr = AA$Asp/AA$Total
AA$CysFr = AA$Cys/AA$Total
AA$GlnFr = AA$Gln/AA$Total
AA$GluFr = AA$Glu/AA$Total
AA$GlyFr = AA$Gly/AA$Total
AA$HisFr = AA$His/AA$Total
AA$IleFr = AA$Ile/AA$Total
AA$LeuFr = AA$Leu/AA$Total
AA$LysFr = AA$Lys/AA$Total
AA$MetFr = AA$Met/AA$Total
AA$PheFr = AA$Phe/AA$Total
AA$ProFr = AA$Pro/AA$Total
AA$SerFr = AA$Ser/AA$Total
AA$ThrFr = AA$Thr/AA$Total
AA$TrpFr = AA$Trp/AA$Total
AA$TyrFr = AA$Tyr/AA$Total
AA$ValFr = AA$Val/AA$Total

AA$GainersFr = AA$Phe + AA$Tyr + AA$Ile + AA$Met
AA$LosersFr  = AA$Pro + AA$Arg + AA$Ala + AA$Gly + AA$Thr + AA$His

plot(AA$X,AA$GainersFr, col=c('red','black',rep('grey',210)))
plot(AA$X,AA$LosersFr, col=c('Green','black',rep('grey',210)))

#####################


library("seqinr")
seq1<- read.fasta(file = "../../Body/1Raw/CovidTemp/SarsCov2.Orf1ab.fasta")
seq2<- read.fasta(file = "../../Body/1Raw/CovidTemp/Ratg13.Orf1ab.fasta")
#seq1<- read.fasta(file = "../../Body/1Raw/CovidTemp/SarsCov2.S.fasta")
#seq2<- read.fasta(file = "../../Body/1Raw/CovidTemp/Ratg13.S.fasta")
seq1string <- toupper(c2s(seq1[[1]]))
seq2string <- toupper(c2s(seq2[[1]]))
library(Biostrings)
globalAlign<- pairwiseAlignment(seq1string, seq2string)
pid(globalAlign, type = "PID3")

summary(globalAlign) # pattern is Sars Cov 2, subject is a bat
data = str(globalAlign)
MisMatchTable = mismatchTable(globalAlign)
MisMatchTable$FromTo = paste(MisMatchTable$SubjectSubstring,'>',MisMatchTable$PatternSubstring,sep='')
AllSubstitutions = data.frame(table(MisMatchTable$FromTo))

Expected = c('A>S','A>V','V>L','P>L','P>S','L>F','R>L','R>C','G>V','G>W','G>C','H>Y','D>Y','S>I','R>I','R>M','T>I','T>M')
Unexpected = c('S>A','V>A','L>V','L>P','S>P','F>L','L>R','C>R','V>G','W>G','C>G','Y>H','Y>D','I>S','I>R','M>R','I>T','M>T')

sum(AllSubstitutions[AllSubstitutions$Var1 %in% Expected,]$Freq)   # 4 S, 19 orf
sum(AllSubstitutions[AllSubstitutions$Var1 %in% Unexpected,]$Freq) # 1 S, 21 orf



