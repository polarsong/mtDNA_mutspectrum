rm(list=ls(all=TRUE))

AA = read.csv('../../Body/3Results/AminoAcidFreqsChordata.csv')
names(AA)
table(AA$Class)
AA = AA[AA$Class == 'Mammalia',]

GenTime = read.table('../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt', sep = '\t', header = TRUE)
GenTime$Species = gsub(' ','_',GenTime$Scientific_name)
GenTime = GenTime[colnames(GenTime) %in% c('Species','Genus','GenerationLength_d')]

dim(AA)
AA = merge(AA,GenTime)
dim(AA)

names(AA)
for (i in 1:nrow(AA))
{
AA$CysFr[i] = AA$Cys[i]/sum(as.numeric(AA[i,3:22]))  
}
summary(AA$CysFr)

boxplot(AA$CysFr ~ AA$Gene)  # ND6 is maximal
table(AA$Gene)

# genes are sorted according to TBSS (from small to high) within Major arc:
cor.test(AA[AA$Gene == 'COX1',]$CysFr,AA[AA$Gene == 'COX1',]$GenerationLength_d) 
cor.test(AA[AA$Gene == 'COX2',]$CysFr,AA[AA$Gene == 'COX2',]$GenerationLength_d)
cor.test(AA[AA$Gene == 'ATP8',]$CysFr,AA[AA$Gene == 'ATP8',]$GenerationLength_d)
cor.test(AA[AA$Gene == 'ATP6',]$CysFr,AA[AA$Gene == 'ATP6',]$GenerationLength_d)
cor.test(AA[AA$Gene == 'COX3',]$CysFr,AA[AA$Gene == 'COX3',]$GenerationLength_d)
cor.test(AA[AA$Gene == 'ND3',]$CysFr,AA[AA$Gene == 'ND3',]$GenerationLength_d)
cor.test(AA[AA$Gene == 'ND4L',]$CysFr,AA[AA$Gene == 'ND4L',]$GenerationLength_d) 
cor.test(AA[AA$Gene == 'ND4',]$CysFr,AA[AA$Gene == 'ND4',]$GenerationLength_d)
cor.test(AA[AA$Gene == 'ND5',]$CysFr,AA[AA$Gene == 'ND5',]$GenerationLength_d)
cor.test(AA[AA$Gene == 'ND6',]$CysFr,AA[AA$Gene == 'ND6',]$GenerationLength_d)  # opposite chain!! did we use reverse complement?
cor.test(AA[AA$Gene == 'CytB',]$CysFr,AA[AA$Gene == 'CytB',]$GenerationLength_d)

# genes are sorted according to TBSS (from small to high) within Minor arc:
cor.test(AA[AA$Gene == 'ND1',]$CysFr,AA[AA$Gene == 'ND1',]$GenerationLength_d)
cor.test(AA[AA$Gene == 'ND2',]$CysFr,AA[AA$Gene == 'ND2',]$GenerationLength_d)




















