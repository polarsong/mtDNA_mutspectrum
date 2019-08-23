rm(list=ls(all=TRUE))

Mut = read.table("../../Body/1Raw/Rebolledo2014/S3.txt", header = TRUE)
Mut$FromTo = paste(Mut$MajorAllele,"_",Mut$MinorAllele,sep = '')
Mut$PosFromTo = paste(Mut$Position,Mut$MajorAllele,">",Mut$MinorAllele,sep = '')
write.table(Mut$PosFromTo,"../../Body/2Derived/Humans.RebolledoAnalyses.ForAnnotation.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
### write the file in order to annotate manually: go to https://mseqdr.org/mvtool.php, paste my variants and save file: "../../Body/1Raw/Rebolledo2014/mseqdr_mvTool_annotation_vep_2019_2_5_15_25.csv"

Ann = read.table("../../Body/1Raw/Rebolledo2014/mseqdr_mvTool_annotation_vep_2019_2_5_15_25.csv", sep = ',', header = TRUE)
table(Ann$consequence_terms)
#    -                   missense_variant non_coding_transcript_exon_variant                 synonymous_variant              upstream_gene_variant 
#   47                                 27                                 15                                 19                                108 
Ann = unique(Ann[,grep("Input|consequence_terms",colnames(Ann))]); names(Ann)=c('PosFromTo','Annotation')
Mut = merge(Mut,Ann, all.x = TRUE)

Age = read.table("../../Body/1Raw/Rebolledo2014/S1.txt", sep = '\t', header = TRUE)
Age$Sample = gsub("-bl",'',Age$BloodId)
Age$length = 0;
for (i in 1:nrow(Age))
{
  Age$length[i] = length(unlist(strsplit(Age$Sample[i],split='')))
}
AgeOfMothers = Age[Age$length <= 4,];  AgeOfMothers$AgeOfMothers = AgeOfMothers$AgeAtCollection; AgeOfMothers=AgeOfMothers[,c(4,6)]
AgeOfKids = Age[Age$length > 4,]; AgeOfKids$Sample = gsub("C1$|C2$|C3$|C4$|C5$",'',AgeOfKids$Sample)
AgeOfKids = AgeOfKids[AgeOfKids$Sample != 'M502G',]
AgeOfKids$AgeOfKids = AgeOfKids$AgeAtCollection; AgeOfKids=AgeOfKids[,c(4,6)]

Mut = merge(Mut,AgeOfMothers, by = 'Sample', all.x = TRUE)
Mut = merge(Mut,AgeOfKids, by = 'Sample', all.x = TRUE)
Mut$MotherMinusKid = Mut$AgeOfMothers - Mut$AgeOfKids
summary(Mut$MotherMinusKid)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 21.00   26.70   32.40   31.13   34.00   39.00       3 

table(Mut$FromTo)
# A_C A_G A_T C_T G_A T_C 
#  1  34   1  14  19  31
table(Mut[Mut$FromTo == 'T_C',]$Category)

Child = Mut[Mut$Category == 'child',]
# Child = Child[Child$Annotation == 'synonymous_variant',] # !!!!!
# Child = Child[Child$Mitomap == 'annotated',] # !!!!!
table(Child$FromTo)
# A_G A_T C_T G_A T_C 
# 6   1   1   4   4 

## compare reproduction age for different mutation types: We can just report this result, no more!!!! Do it. 
## Cite the main conclusion of the paper (reread it and say that T_C might be the strongest driver, because the maximal mean age)
## report the trend, not statistical result - it is difficult to get significant result with so low numbers.

mean(Child[Child$FromTo == 'T_C',]$MotherMinusKid); nrow(Child[Child$FromTo == 'T_C',]) # 34.075
mean(Child[Child$FromTo != 'T_C',]$MotherMinusKid); nrow(Child[Child$FromTo != 'T_C',])  # 32.066
mean(Child[Child$FromTo == 'A_G',]$MotherMinusKid) # 33.7
mean(Child[Child$FromTo == 'G_A',]$MotherMinusKid); nrow(Child[Child$FromTo == 'G_A',]) # 29.9
mean(Child[Child$FromTo == 'A_T',]$MotherMinusKid) # 37.9 - just one
mean(Child[Child$FromTo == 'C_T',]$MotherMinusKid) # 25.1 - just one

wilcox.test(Child[Child$FromTo == 'T_C',]$MotherMinusKid,Child[Child$FromTo == 'G_A',]$MotherMinusKid, alternative = 'greater')  # p = 0.0956 
wilcox.test(Child[Child$FromTo == 'T_C',]$MotherMinusKid, Child[Child$FromTo != 'T_C',]$MotherMinusKid, alternative = 'greater') # p = 0.19
wilcox.test(Child[Child$FromTo == 'T_C' | Child$FromTo == 'A_G',]$MotherMinusKid, Child[Child$FromTo != 'T_C' || Child$FromTo != 'A_G',]$MotherMinusKid, alternative = 'greater') # p = 0.19

TCGA = Child[Child$FromTo == 'T_C' | Child$FromTo == 'G_A',]
summary(Child$MotherMinusKid) # 33.55
nrow(Child[Child$FromTo == 'T_C' & Child$MotherMinusKid < median(Child$MotherMinusKid),])/nrow(Child[Child$FromTo == 'G_A' & Child$MotherMinusKid < median(Child$MotherMinusKid),])   # 1 to 3
nrow(Child[Child$FromTo == 'T_C' & Child$MotherMinusKid >= median(Child$MotherMinusKid),])/nrow(Child[Child$FromTo == 'G_A' & Child$MotherMinusKid >= median(Child$MotherMinusKid),]) # 3 vs 1

## Take a subset of kids with several de novo mtDNA mutations and compare them according to their VAFs:


## somatic mutations in mothers versus their age - no difference:
SomGain = Mut[Mut$Category == 'somatic-gain' & !is.na(Mut$AgeOfMothers),]
mean(SomGain[SomGain$FromTo == 'T_C',]$AgeOfMothers); nrow(SomGain[SomGain$FromTo == 'T_C',]);
mean(SomGain[SomGain$FromTo != 'T_C',]$AgeOfMothers); nrow(SomGain[SomGain$FromTo != 'T_C',]);

