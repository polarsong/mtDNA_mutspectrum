#####  COUNT NUCLEOTIDE CONTENT CAREFULLY (BODY/2Derived/polarizedbr_data => external + More Shallow => codons, 4fold nucl, FrA,T,G,C) => barplot?
#####  normalized average MutSpec (pie charts or 12 boxplots for each class)?

rm(list=ls(all=TRUE))

MS = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.AllSynMutSepForEachGene.txt', header = TRUE)
MS$Species = gsub("\\.(.*)",'',MS$SpeciesGene)
MS$Gene = gsub("(.*)\\.",'',MS$SpeciesGene)

GT = read.table('../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt', header = TRUE, sep = '\t')
GT$Species = gsub(" ",'_',GT$Scientific_name)
  
MS = merge(MS,GT) 

MS = MS[MS$Gene == 'COX1' | MS$Gene == 'CytB',]
