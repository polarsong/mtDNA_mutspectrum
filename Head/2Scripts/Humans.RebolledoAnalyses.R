rm(list=ls(all=TRUE))

Mut = read.table("../../Body/1Raw/Rebolledo2014/S3.txt", header = TRUE)
Mut$FromTo = paste(Mut$MajorAllele,"_",Mut$MinorAllele,sep = '')
Mut$PosFromTo = paste(Mut$Position,Mut$MajorAllele,">",Mut$MinorAllele,sep = '')
write.table(Mut$PosFromTo,"../../Body/2Derived/Humans.RebolledoAnalyses.ForAnnotation.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

###### go to https://mseqdr.org/mvtool.php, paste my variants and save file: "../../Body/1Raw/Rebolledo2014/mseqdr_mvTool_annotation_vep_2019_2_5_15_25.csv"
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

pdf("../../Body/4Figures/Humans.RebolledoAnalyses.R.01.pdf", width = 14, height = 14)
par(mfrow=c(2,3)) # dev.off()

################## somatic-gain => mutations which originated once - either in mother's one tissue or in child's one tissues - we have to split them:
table(Mut[Mut$Category == 'somatic-gain',]$FromTo)
# A_G C_T G_A T_C 
#  7   1   2   5 
SomGain = Mut[Mut$Category == 'somatic-gain',]
SomGain[,c(6:9)] = sapply(SomGain[,(6:9)], as.character)
SomGain[,c(6:9)] = sapply(SomGain[,(6:9)], as.numeric)
SomGainMothers = SomGain[(SomGain$MotherCheck + SomGain$MotherBlood) > (SomGain$ChildCheck + SomGain$ChildBlood),]
SomGainKids = SomGain[(SomGain$MotherCheck + SomGain$MotherBlood) < (SomGain$ChildCheck + SomGain$ChildBlood),]
table(SomGainMothers$FromTo)
# A_G G_A T_C 
#  4   1   4 
boxplot(SomGainMothers[SomGainMothers$FromTo == 'T_C',]$MotherMinusKid, SomGainMothers[SomGainMothers$FromTo != 'T_C',]$MotherMinusKid, names = c('T>C','TheRest'), ylab = 'AgeOfMotherMinusKid', main = '"somatic-gain" in mothers')
wilcox.test(SomGainMothers[SomGainMothers$FromTo == 'T_C',]$MotherMinusKid, SomGainMothers[SomGainMothers$FromTo != 'T_C',]$MotherMinusKid) # 0.87

# plot.new() # dev.off()

############# child - de novo origin in kids
Child = Mut[Mut$Category == 'child',]
table(Child$FromTo)
# A_G A_T C_T G_A T_C 
# 6   1   1   4   4 
boxplot(Child[Child$FromTo == 'T_C',]$MotherMinusKid, Child[Child$FromTo != 'T_C',]$MotherMinusKid, names = c('T>C','TheRest'), ylab = 'AgeOfMotherMinusKid', main = 'de novo germline in "child"')
wilcox.test(Child[Child$FromTo == 'T_C',]$MotherMinusKid, Child[Child$FromTo != 'T_C',]$MotherMinusKid, alternative = 'greater') #

YoungTc = nrow(Child[Child$MotherMinusKid < median(Child$MotherMinusKid) & Child$FromTo == 'T_C',]); YoungTc       # 1
YoungNotTc = nrow(Child[Child$MotherMinusKid < median(Child$MotherMinusKid) & Child$FromTo != 'T_C',]); YoungNotTc # 7
OldTc = nrow(Child[Child$MotherMinusKid >= median(Child$MotherMinusKid) & Child$FromTo == 'T_C',]); OldTc            # 3       
OldNotTc = nrow(Child[Child$MotherMinusKid >= median(Child$MotherMinusKid) & Child$FromTo != 'T_C',]); OldNotTc      # 5 
X = cbind(c(OldTc,OldNotTc),c(YoungTc,YoungNotTc))
fisher.test(X)  # 
mosaicplot(X) # dev.off()

YoungTc = nrow(Child[Child$MotherMinusKid < median(Child$MotherMinusKid) & Child$FromTo == 'T_C',]); YoungTc       # 1
YoungGa = nrow(Child[Child$MotherMinusKid < median(Child$MotherMinusKid) & Child$FromTo == 'G_A',]); YoungGa       # 3
OldTc = nrow(Child[Child$MotherMinusKid >= median(Child$MotherMinusKid) & Child$FromTo == 'T_C',]); OldTc          # 3       
OldGa = nrow(Child[Child$MotherMinusKid >= median(Child$MotherMinusKid) & Child$FromTo == 'G_A',]); OldGa          # 1 
X = cbind(c(OldTc,OldGa),c(YoungTc,YoungGa))
fisher.test(X)  # 
mosaicplot(X) # dev.off()

############ merge together "child" + "somatic gain in mothers" 
Merge = Child
# Merge = rbind(SomGainMothers,Child)
Merge = Merge[!is.na(Merge$MotherMinusKid),] # there is one unknown age for kid (SC14)!!!! And this is T>C!!!! somatic gain in mother!!! may be ask Makova!!
table(Merge$FromTo)
# A_G A_T C_T G_A T_C 
# 10   1   1   5   8 

# remove missense variants
# Merge = Merge[Merge$Annotation != 'missense_variant',] # 25 => 14

boxplot(Merge[Merge$FromTo == 'T_C',]$MotherMinusKid,Merge[Merge$FromTo != 'T_C',]$MotherMinusKid, notch = TRUE, names = c('T>C','TheRest'), ylab = 'AgeOfMotherMinusKid', main = 'child + somatic gain in mothers') # dev.off()
wilcox.test(Merge[Merge$FromTo == 'T_C',]$MotherMinusKid,Merge[Merge$FromTo != 'T_C',]$MotherMinusKid, alternative = 'greater') # 0.24

boxplot(Merge[Merge$FromTo == 'T_C',]$MotherMinusKid,Merge[Merge$FromTo == 'A_G',]$MotherMinusKid,Merge[Merge$FromTo == 'G_A',]$MotherMinusKid, names = c('T>C','A>G','G>A'), ylab = 'AgeOfMotherMinusKid', main = 'child + somatic gain in mothers')
wilcox.test(Merge[Merge$FromTo == 'T_C',]$MotherMinusKid,Merge[Merge$FromTo == 'A_G',]$MotherMinusKid)
wilcox.test(Merge[Merge$FromTo == 'T_C',]$MotherMinusKid,Merge[Merge$FromTo == 'G_A',]$MotherMinusKid, alternative = 'greater') # 0.095!!!!

boxplot(Merge[Merge$FromTo == 'T_C' | Merge$FromTo == 'A_G',]$MotherMinusKid,Merge[Merge$FromTo != 'T_C' & Merge$FromTo != 'A_G',]$MotherMinusKid, names = c('T>C & A>G','TheRest '), ylab = 'AgeOfMotherMinusKid', main = 'child + somatic gain in mothers')
wilcox.test(Merge[Merge$FromTo == 'T_C' | Merge$FromTo == 'A_G',]$MotherMinusKid,Merge[Merge$FromTo != 'T_C' & Merge$FromTo != 'A_G',]$MotherMinusKid, alternative = 'greater')

##### try logistic regression - nothing
Merge$DummyTC = 0;
for (i in 1:nrow(Merge)) {if (Merge$FromTo[i] == 'T_C') {Merge$DummyTC[i] = 1;}}
A<-glm(Merge$DummyTC ~  scale(Merge$MotherMinusKid), family = binomial()); summary(A)

#### what else - delete deleterious and work only with synonymous? but I have no power finally...
# I can delete deleterious variants, discussed in the paper (table 2):
VectorOfLocationsOfDeleteriousMutations = c(195, 1391, 1555, 2352, 3242, 3243, 12634, 13708) # table 2, paper
# Merge = Merge[!Merge$Position %in% VectorOfLocationsOfDeleteriousMutations,] # 3 - out doesn't affect significantly...

#### what else - permute dataset? 
ObservedMean = mean(Merge[Merge$FromTo == 'A_G',]$MotherMinusKid); ObservedMean;  # 32.75
ObservedMean = mean(Merge[Merge$FromTo == 'G_A',]$MotherMinusKid); ObservedMean;  # 30.4
ObservedMean = mean(Merge[Merge$FromTo == 'T_C',]$MotherMinusKid); ObservedMean;  # 32.78
PermMeanFinal = c()
for (i in 1:10000)
{
  Merge$FromToPerm = sample(Merge$FromTo)
  PermMean = mean(Merge[Merge$FromToPerm == 'T_C',]$MotherMinusKid); 
  PermMeanFinal = c(PermMeanFinal,PermMean)
}
summary(PermMeanFinal)
hist(PermMeanFinal, breaks = 100)
abline(v=ObservedMean, col = 'red', lwd = 4)
length(PermMeanFinal[PermMeanFinal>ObservedMean])/length(PermMeanFinal) # 0.223

#### permute only T>C and G>A
Merge = Merge[Merge$FromTo == 'T_C' | Merge$FromTo == 'G_A',]
ObservedMean = mean(Merge[Merge$FromTo == 'A_G',]$MotherMinusKid); ObservedMean;  # NA
ObservedMean = mean(Merge[Merge$FromTo == 'G_A',]$MotherMinusKid); ObservedMean;  # 29.9
ObservedMean = mean(Merge[Merge$FromTo == 'T_C',]$MotherMinusKid); ObservedMean;  # 34.075
PermMeanFinal = c()
for (i in 1:10000)
{
  Merge$FromToPerm = sample(Merge$FromTo)
  PermMean = mean(Merge[Merge$FromToPerm == 'T_C',]$MotherMinusKid); 
  PermMeanFinal = c(PermMeanFinal,PermMean)
}
summary(PermMeanFinal)
hist(PermMeanFinal, breaks = 100)
abline(v=ObservedMean, col = 'red', lwd = 4)
length(PermMeanFinal[PermMeanFinal>ObservedMean])/length(PermMeanFinal) # 0.0905

dev.off()
