###################################

rm(list=ls(all=TRUE))
  
ALL = read.table("../../Body/1Raw/mtDNA_snv_Oct2016.txt", head = TRUE, sep = '\t')  
length(unique(ALL$sample)) # 2177

Pat = read.table("../../Body/1Raw/TCGA_MTmut_Age_Coverage.txt", head = TRUE, sep = '\t')  

Pat = Pat[,grepl("sample|Consensus_age|Normal_coverage|Tumour_coverage|mt_coverage|mt_copies", names(Pat))]  # don't need tissue and Tier2
nrow(Pat)
length(unique(Pat$sample)) # 2176
setdiff(unique(ALL$sample),unique(Pat$sample)) # PCSI_0357 -> for one sample ID we have no age / mtDNA copies info

ALL = merge(ALL,Pat, by = 'sample', all.x = TRUE)
names(ALL)

write.table(ALL,"../../Body/2Derived/mtDNA_snv_Oct2016.PatientInfo.txt", sep = '\t', row.names = FALSE, quote = FALSE) 