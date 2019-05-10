rm(list=ls(all=TRUE))

Som = read.table("../../Body/2Derived/HealthySomaticHumanMutations.GTEx.PatientSpecificDerive.txt", header = TRUE, sep = '\t')

########################################
######### 2: ANALYSE TABLE Som: compare mut specs of pairs of tissues within the same individuum
########################################

summary(Som$AF)
# Som = Som[Som$AF >= quantile(Som$AF,0.75),]
Som = Som[Som$AF <= quantile(Som$AF,0.5),]

# Som = Som[Som$Substitution == 'T_C' | Som$Substitution == 'G_A',]  # !!!!!!!!!!!!

VecOfPatients = unique(as.character(Som$subject)); length(VecOfPatients) # 435
flag = 0;

for (i in 1:length(VecOfPatients))
{ # i = 3
temp = Som[Som$subject == VecOfPatients[i],]
VecOfTissues = as.character(unique(temp$TissueShortName))
if (length(VecOfTissues) > 1)
{
  for (tissue1 in 1:(length(VecOfTissues)-1))
  { # tissue1 = 2
    for (tissue2 in (tissue1+1):length(VecOfTissues))
    { # tissue2 = 3
      mut1 = temp[temp$TissueShortName == VecOfTissues[tissue1],]
      mut2 = temp[temp$TissueShortName == VecOfTissues[tissue2],]
      Number1 = nrow(mut1)
      Number2 = nrow(mut2)
      MutSpek1 = nrow(mut1[mut1$Substitution == 'G_A',])/nrow(mut1) # G_A, T_C  #!!!!!!!!!!!!!!!!!!
      MutSpek2 = nrow(mut2[mut2$Substitution == 'G_A',])/nrow(mut2) # G_A, T_C  #!!!!!!!!!!!!!!!!!!
      TurnOver1 = mut1$TurnOverRate[1]
      TurnOver2 = mut2$TurnOverRate[1]
      OneLine = c(VecOfPatients[i],Number1,Number2,MutSpek1,MutSpek2,TurnOver1,TurnOver2)
      OneLine = as.data.frame(t(OneLine))
      names(OneLine)=c('SubjId','Number1','Number2','MutSpek1','MutSpek2','TurnOver1','TurnOver2')
      
      if (flag == 1) {Final = rbind(Final,OneLine);}
      if (flag == 0) {Final = OneLine; flag = 1;}
      }
  }
}
}

str(Final)
Final$SubjId = as.character(Final$SubjId)
Final[2:7] <- lapply(Final[2:7], function(x) { if(is.factor(x)) as.numeric(as.character(x)) else x }) 
str(Final)
Final$DeltaTurnOver = Final$TurnOver1-Final$TurnOver2
Final$DeltaMutSpec  = Final$MutSpek1 - Final$MutSpek2

### all pairs:
cor.test(Final$DeltaMutSpec,Final$DeltaTurnOver, method = 'spearman') # negative!
plot(Final$DeltaTurnOver,Final$DeltaMutSpec)

### pairs with at least 2 mutations:
cor.test(Final[Final$Number1 >= 2 & Final$Number2 >= 2,]$DeltaMutSpec,Final[Final$Number1 >= 2 & Final$Number2 >= 2,]$DeltaTurnOver, method = 'spearman') # negative!
plot(Final[Final$Number1 >= 2 & Final$Number2 >= 2,]$DeltaTurnOver,Final[Final$Number1 >= 2 & Final$Number2 >= 2,]$DeltaMutSpec)

###### if we split by individuals (many pairs versus a few):
Som$Number = 1; Agg = aggregate(Som$Number, by = list(Som$subject), FUN = sum); summary(Agg$x)
UpperQuartile = Agg[Agg$x >= 8,]$Group.1; length(UpperQuartile)

cor.test(Final[Final$SubjId %in% UpperQuartile,]$DeltaMutSpec,Final[Final$SubjId %in% UpperQuartile,]$DeltaTurnOver, method = 'spearman') # negative!
cor.test(Final[!Final$SubjId %in% UpperQuartile,]$DeltaMutSpec,Final[!Final$SubjId %in% UpperQuartile,]$DeltaTurnOver, method = 'spearman') # negative!



