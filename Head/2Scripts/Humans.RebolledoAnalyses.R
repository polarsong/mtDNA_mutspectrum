rm(list=ls(all=TRUE))

Mut = read.table("../../Body/1Raw/Rebolledo2014/S3.txt", header = TRUE)
Age = read.table("../../Body/1Raw/Rebolledo2014/S1.txt", sep = '\t', header = TRUE)
Age$Sample = gsub("-bl",'',Age$BloodId)
Mut = merge(Mut,Age, by = 'Sample')

