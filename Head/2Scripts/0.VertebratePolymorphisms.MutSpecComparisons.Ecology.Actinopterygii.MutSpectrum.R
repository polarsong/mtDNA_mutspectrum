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
MUTSPEC= MUTACTINOPTERITAXAFROMK
MUTACTINOPTERITAXAFROMK = MUTACTINOPTERITAXAFROMK[MUTACTINOPTERITAXAFROMK$Class == "Actinopterygii",]
#names(MUTACTINOPTERITAXAFROMK) = c("T_G", "T_C", "T_A", "G_T", "G_C", "G_A", "C_T", "C_G", "C_A", "A_T", "A_G", "A_C")



MUTACTINOPTERITAXAFROMK=MUTACTINOPTERITAXAFROMK[,-c(1,14)]
temp = data.frame(t(MUTACTINOPTERITAXAFROMK))
a=c()

for(i in 1:ncol(temp)){
  b = data.frame(temp[i], rownames(temp))
  colnames(b)=c("Freq", "Sub")
  a = rbind(a, b)
}

n=names(MUTACTINOPTERITAXAFROMK)
m=c()
i="A_T"
for(i in n){
  t=c(i, gsub(":", "", strsplit(summary(MUTACTINOPTERITAXAFROMK[i])[4], split=" ")[[1]][4]))
  m=rbind(m, t)
}
m=data.frame(m, row.names = NULL)
summary(MUTACTINOPTERITAXAFROMK)

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



#####Absolute numbers
ABSOLUTEMUT = read.table('../../Body/2Derived/VertebratePolymorphisms.MutSpecData.txt', header = TRUE)
AGG = ABSOLUTEMUT[ABSOLUTEMUT$MutType == "FourFold",]; AGG$N = 1
AGG = aggregate(AGG$N, by=list(AGG$Species, AGG$Subs), FUN=sum)

nam=unique(AGG$Group.1)

final=c()
for (i in nam){
  TEMP=AGG[AGG$Group.1 == i,]
  TEMP=t(TEMP)
  TEMP2=data.frame()
  TEMP2=rbind(TEMP2, TEMP[-c(1,2),])
  names(TEMP2)=TEMP[2,]
  TEMP2$Species=i
  final=dplyr::bind_rows(final, TEMP2)
  print(i)
}
final[is.na(final)] <- 0
str(final)

MUTSPEC = merge(MUTSPEC, final, by="Species")
str(MUTSPEC)
summary(MUTSPEC)

#######check for zeroes
for (i in 1:nrow(MUTSPEC)){
  if (MUTSPEC$A_T.x[i] == 0){
    MUTSPEC$Check[i] = print(as.numeric(MUTSPEC$A_T.x[i]) == as.numeric(MUTSPEC$A_T.y[i]))
  }
  if (MUTSPEC$A_G.x[i] == 0){
    MUTSPEC$Check[i] = print(as.numeric(MUTSPEC$A_G.x[i]) == as.numeric(MUTSPEC$A_G.y[i]))
  }
  if (MUTSPEC$A_C.x[i] == 0){
    MUTSPEC$Check[i] = print(as.numeric(MUTSPEC$A_C.x[i]) == as.numeric(MUTSPEC$A_C.y[i]))
  }
  if (MUTSPEC$C_A.x[i] == 0){
    MUTSPEC$Check[i] = print(as.numeric(MUTSPEC$C_A.x[i]) == as.numeric(MUTSPEC$C_A.y[i]))
  }
  if (MUTSPEC$C_G.x[i] == 0){
    MUTSPEC$Check[i] = print(as.numeric(MUTSPEC$C_G.x[i]) == as.numeric(MUTSPEC$C_G.y[i]))
  }
  if (MUTSPEC$C_T.x[i] == 0){
    MUTSPEC$Check[i] = print(as.numeric(MUTSPEC$C_T.x[i]) == as.numeric(MUTSPEC$C_T.y[i]))
  }
  if (MUTSPEC$G_A.x[i] == 0){
    MUTSPEC$Check[i] = print(as.numeric(MUTSPEC$G_A.x[i]) == as.numeric(MUTSPEC$G_A.y[i]))
  }
  if (MUTSPEC$G_T.x[i] == 0){
    MUTSPEC$Check[i] = print(as.numeric(MUTSPEC$G_T.x[i]) == as.numeric(MUTSPEC$G_T.y[i]))
  }
  if (MUTSPEC$G_C.x[i] == 0){
    MUTSPEC$Check[i] = print(as.numeric(MUTSPEC$G_C.x[i]) == as.numeric(MUTSPEC$G_C.y[i]))
  }
  if (MUTSPEC$T_A.x[i] == 0){
    MUTSPEC$Check[i] = print(as.numeric(MUTSPEC$T_A.x[i]) == as.numeric(MUTSPEC$T_A.y[i]))
  }
  if (MUTSPEC$T_G.x[i] == 0){
    MUTSPEC$Check[i] = print(as.numeric(MUTSPEC$T_G.x[i]) == as.numeric(MUTSPEC$T_G.y[i]))
  }
  if (MUTSPEC$T_C.x[i] == 0){
    MUTSPEC$Check[i] = print(as.numeric(MUTSPEC$T_C.x[i]) == as.numeric(MUTSPEC$T_C.y[i]))
  }
}

table(MUTSPEC$Check)
MUTSPEC=MUTSPEC[,-27]
str(MUTSPEC)
write.table(MUT, file = '../../Body/3Results/VertebratePolymorphisms.MutSpecData.with.AbsoluteNumbers.txt', quote = FALSE)
