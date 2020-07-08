rm(list=ls(all=TRUE))

if (!require(ggpubr)) install.packages("ggpubr")

library("ggpubr")

###########Taxonomy###################################################################
Taxa = read.table("../../Body/1Raw/TaxaFromKostya.Names.stat", sep = '\t',header = FALSE) 
Taxa$Species = gsub(";.*",'',Taxa$V1); 
for (i in (1:nrow(Taxa)))  {Taxa$Species[i] = paste(unlist(strsplit(as.character(Taxa$Species[i]),split = ' '))[1],unlist(strsplit(as.character(Taxa$Species[i]),split = ' '))[2], sep = '_')}
Taxa$Class = gsub(";Chordata;.*",'',Taxa$V1); Taxa$Class = gsub(".*;",'',Taxa$Class); table(Taxa$Class)
Taxa$Class = gsub('Actinopteri','Actinopterygii',Taxa$Class)
Taxa$Class = gsub("Testudines|Squamata|Crocodylia|Sphenodontia",'Reptilia',Taxa$Class)
length(unique(Taxa$Species)) # 1708
table(Taxa$Class)
Taxa = Taxa[,-1]

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
#########################################################################################

###########################MUT spectrum in Actinoptery##################################
MUT = read.table('../../Body/3Results/VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.txt', header = TRUE)
MUTACTINOPTERITAXAFROMK = merge(MUT, Taxa)
MUTACTINOPTERITAXAFROMK = MUTACTINOPTERITAXAFROMK[MUTACTINOPTERITAXAFROMK$Class == "Actinopterygii",]
#names(MUTACTINOPTERITAXAFROMK) = c("T_G", "T_C", "T_A", "G_T", "G_C", "G_A", "C_T", "C_G", "C_A", "A_T", "A_G", "A_C")

MUTACTINOPTERITAXAFROMK=MUTACTINOPTERITAXAFROMK[,-1]
MUTACTINOPTERITAXAFROMK=MUTACTINOPTERITAXAFROMK[,-13]
temp = data.frame(t(MUTACTINOPTERITAXAFROMK))
a=c()

for(i in 1:ncol(temp)){
  b = data.frame(temp[i], rownames(temp))
  colnames(b)=c("Freq", "Sub")
  a = rbind(a, b)
}

#c("TH>GH", "TH>CH", "TH>AH", "GH>TH", "GH>CH", "GH>AH", "CH>TH", "CH>GH", "CH>AH", "AH>TH", "AH>GH", "AH>CH")
pdf("../../Body/4Figures/VertebratePolymorphisms.MutSpecComparisons.Analyses.Ecology.Actinopterygii.FishBaseData.FIGURE1A.pdf")
ggbarplot(a, x = "Sub", y = "Freq", fill = "Sub", color = "Sub",
          palette = c("#bdbdbd", "#73514f", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#055088", "#9c3d37", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#036a5b", "#bdbdbd"), 
          xlab="Substitution types", ylab="Normalised frequencies", add = "mean_se", legend = "none")  
dev.off()
##########################################################################################

svg("../../Body/4Figures/VertebratePolymorphisms.MutSpecComparisons.Analyses.Ecology.Actinopterygii.FishBaseData.FIGURE1A.svg")
ggbarplot(a, x = "Sub", y = "Freq", fill = "Sub", color = "Sub",
          palette = c("#bdbdbd", "#73514f", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#055088", "#9c3d37", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#036a5b", "#bdbdbd"), 
          xlab="Substitution types", ylab="Normalised frequencies", add = "mean_se")
dev.off()

