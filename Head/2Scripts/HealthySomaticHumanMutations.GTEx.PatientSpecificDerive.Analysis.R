
rm(list=ls(all=TRUE))

Som = read.table("../../Body/2Derived/HealthySomaticHumanMutations.GTEx.PatientSpecificDerive.txt", header = TRUE, sep = '\t')

########################################
######### 2: ANALYSE TABLE Som: compare mut specs of pairs of tissues within the same individuum
########################################

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
      MutSpek1 = nrow(mut1[mut1$Substitution == 'G_A',])/nrow(mut1)
      MutSpek2 = nrow(mut2[mut2$Substitution == 'G_A',])/nrow(mut2)
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


