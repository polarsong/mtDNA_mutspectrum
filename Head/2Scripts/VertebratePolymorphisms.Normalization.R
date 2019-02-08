#####  COUNT NUCLEOTIDE CONTENT CAREFULLY (BODY/2Derived/polarizedbr_data => external + More Shallow => codons, 4fold nucl, FrA,T,G,C) => barplot?
#####  normalized average MutSpec (pie charts or 12 boxplots for each class)?

rm(list=ls(all=TRUE))

VecOfSynFourFoldDegenerateSites <- c('CTT', 'CTC', 'CTA', 'CTG', 
                                     'GTT', 'GTC', 'GTA', 'GTG', 
                                     'TCT', 'TCC', 'TCA', 'TCG', 
                                     'CCT', 'CCC', 'CCA', 'CCG', 
                                     'ACT', 'ACC', 'ACA', 'ACG', 
                                     'GCT', 'GCC', 'GCA', 'GCG', 
                                     'CGT', 'CGC', 'CGA', 'CGG', 
                                     'GGT', 'GGC', 'GGA', 'GGG')
length(unique(VecOfSynFourFoldDegenerateSites)) # 32

unzip("../../Body/2Derived/POLARIZEDBR_DATA.zip", exdir= "../../Body/2Derived/")
List = list.files("../../Body/2Derived/POLARIZEDBR_DATA/")

for (i in 1:length(List))
{ # i = 1303
  infile = paste("../../Body/2Derived/POLARIZEDBR_DATA/",as.character(List[i]),sep='') 
  if (length(grep('POLARISED',infile)) > 0)
  {
  Species = gsub('\\..*','',as.character(List[i]))
  Gene = gsub(Species,'',as.character(List[i])); Gene = gsub('.POLARISED.txt','',Gene); Gene = gsub('\\.POLARISED.txt','',Gene); Gene = gsub('\\.','',Gene); 
  GeneSpecies = read.table(infile, header = TRUE)
  GeneSpecies = GeneSpecies[GeneSpecies$BranchPosition == 'External',]
  ExternalSeqsTogether = paste(GeneSpecies$MoreShallowNodeSeq,collapse = '')
  ExternalSeqsTogether = unlist(strsplit(ExternalSeqsTogether,'')) # 5700/3
  CodonsVec = c(); StartNuc = 1
  if (length(ExternalSeqsTogether)/3 == round(length(ExternalSeqsTogether)/3))  # if divide by 3 without the rest
  {
  for (j in 1:(length(ExternalSeqsTogether)/3))
    {
    CodonsVec = c(CodonsVec,paste(ExternalSeqsTogether[StartNuc : (StartNuc+2)],collapse = ''))
    StartNuc = StartNuc+3
    }
  AllCodons = length(CodonsVec)        # 1021
  CodonsVecNeutral = CodonsVec[CodonsVec %in% VecOfSynFourFoldDegenerateSites]
  NeutralCodons = length(CodonsVecNeutral) # 1900
  data.frame(table(CodonsVecNeutral))
  
  CodonsVecNeutral = gsub("CTA|GTA|TCA|CCA|ACA|GCA|CGA|GGA",'A',CodonsVecNeutral)
  CodonsVecNeutral = gsub("CTT|GTT|TCT|CCT|ACT|GCT|CGT|GGT",'T',CodonsVecNeutral)
  CodonsVecNeutral = gsub("CTG|GTG|TCG|CCG|ACG|GCG|CGG|GGG",'G',CodonsVecNeutral)
  CodonsVecNeutral = gsub("CTC|GTC|TCC|CCC|ACC|GCC|CGC|GGC",'C',CodonsVecNeutral)
  
  Line=c(Species,Gene,length(CodonsVecNeutral[CodonsVecNeutral == 'A']),length(CodonsVecNeutral[CodonsVecNeutral == 'T']),length(CodonsVecNeutral[CodonsVecNeutral == 'G']),length(CodonsVecNeutral[CodonsVecNeutral == 'C']), AllCodons, NeutralCodons)
  if (i == 1) {Final = Line}
  if (i >  1) {Final = rbind(Final,Line)}
  }
  }
}

Final = as.data.frame(Final); names(Final)=c('Species','Gene','A','T','G','C',"NumberOfFourFoldDegenCodons",'NumberOfAllCodons')
write.table(Final, "../../Body/3Results/VertebratePolymorphisms.Normalization.NeutralATGC.txt", quote = FALSE, row.names = FALSE)

## delete all unziped files
files <- list.files("../../Body/2Derived/POLARIZEDBR_DATA/")
for (i in 1:length(files)) 
{ # i = 1
  file = paste('../../Body/2Derived/POLARIZEDBR_DATA/',files[i],sep='')
  if (file.exists(file)) file.remove(file)
}

################### merge with classes (from Taxa & MoreTaxa), average A T G C for each species, average it for classes and draw it.
#### associate species name with Class
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

###

Final = merge(Final,Taxa)
str(Final)
Final[,3:6] <- sapply(Final[,3:6], function(x) as.numeric(as.character(x)))
Agg1 = aggregate(list(Final$A,Final$T,Final$G,Final$C), by = list(Final$Species, Final$Class), FUN = sum); names(Agg1) = c('Species','Class','A','T','G','C')

Agg1$FrA = Agg1$A/(Agg1$A+Agg1$T+Agg1$G+Agg1$C)
Agg1$FrT = Agg1$T/(Agg1$A+Agg1$T+Agg1$G+Agg1$C)
Agg1$FrG = Agg1$G/(Agg1$A+Agg1$T+Agg1$G+Agg1$C)
Agg1$FrC = Agg1$C/(Agg1$A+Agg1$T+Agg1$G+Agg1$C)

Agg2 = aggregate(list(Agg1$FrA,Agg1$FrT,Agg1$FrG,Agg1$FrC), by = list(Agg1$Class), FUN = mean); names(Agg2) = c('Class','FrA','FrT','FrG','FrC')
Agg2$Sum = Agg2$FrA + Agg2$FrT + Agg2$FrG + Agg2$FrC # 1 
Agg2 = Agg2[Agg2$Class %in% c('Actinopterygii','Amphibia','Reptilia','Mammalia','Aves'),]

row.names(Agg2)=Agg2$Class
BARPLOT = t(as.matrix(Agg2[,c(2:5)]))
BARPLOT = BARPLOT[,c(1,5,2,4,3)] # column order
BARPLOT = BARPLOT[c(3,2,4,1),] # row order

colnames(data)=c("A","B","C","D","E")
rownames(data)=c("var1","var2","var3")

# Get the stacked barplot
ColG = rgb(0.1,0.1,0.1,0.5)
ColT = rgb(0.1,0.1,1,0.5)
ColC = rgb(0.1,1,0.1,0.5)
ColA = rgb(1,0.1,0.1,0.5)

pdf("../../Body/4Figures/VertebratePolymorphisms.Normalization.R.01.pdf", width = 30, height = 20)

barplot(BARPLOT, col = c(ColG,ColT,ColC,ColA), border="white", space=0.04, font.axis=2, xlab="")

PieChartTable = read.table("../../Body/3Results/VertebratePolymorphisms.BetweenClassesWithoutNormalization.PieChartTable.txt")
PieChartTable$AncestralNuc = NA
for (i in 1:nrow(PieChartTable)) {PieChartTable$AncestralNuc[i] = unlist(strsplit(as.character(PieChartTable$Subs[i]),split = '_'))[1]}

PieChartTable = merge(PieChartTable,Agg2[,c(1:5)])
PieChartTable$Normalised1Number = 0
for (i in 1:nrow(PieChartTable))
{
  if (PieChartTable$AncestralNuc[i] == 'A') {PieChartTable$Normalised1Number[i] = PieChartTable$Number[i] / PieChartTable$FrA[i]}
  if (PieChartTable$AncestralNuc[i] == 'T') {PieChartTable$Normalised1Number[i] = PieChartTable$Number[i] / PieChartTable$FrT[i]}
  if (PieChartTable$AncestralNuc[i] == 'G') {PieChartTable$Normalised1Number[i] = PieChartTable$Number[i] / PieChartTable$FrG[i]}
  if (PieChartTable$AncestralNuc[i] == 'C') {PieChartTable$Normalised1Number[i] = PieChartTable$Number[i] / PieChartTable$FrC[i]}
}
 
Agg = aggregate(PieChartTable$Normalised1Number, by = list(PieChartTable$Class), FUN = sum); names(Agg) = c('Class','Total')
PieChartTable = merge(PieChartTable,Agg)
PieChartTable$Normalised2Number = PieChartTable$Normalised1Number/PieChartTable$Total

write.table(PieChartTable,"../../Body/3Results/VertebratePolymorphisms.Normalization.Normalized12Fractions.txt", quote = FALSE, row.names = FALSE)

par(mfrow=c(1,5))
VecOfClasses = c('Actinopterygii','Amphibia','Reptilia','Mammalia','Aves')
for (i in 1:length(VecOfClasses))
{
pie(PieChartTable[PieChartTable$Class == VecOfClasses[i],]$Normalised2Number, labels = PieChartTable[PieChartTable$Class == VecOfClasses[i],]$Subs, main = VecOfClasses[i], col=rainbow(12))
}

dev.off()

