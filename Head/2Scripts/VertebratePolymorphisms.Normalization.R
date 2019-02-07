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