rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

wd <- getwd()
wd_in = gsub('Head/2Scripts','Body/2Derived/', wd)  #SWITCH to MP-prefix if MPanalysis data
wd_out = gsub('Head/2Scripts','Body/3Results/', wd) 
setwd(wd_in)
wd_in



table <- read.table("1.txt", header=TRUE, colClasses=c('character'))

FirstCodon = table$AncestorCodon
SecondCodon = table$DescendantCodon
data = data.frame(FirstCodon,SecondCodon)
data$FirstSecond = paste(data$FirstCodon,data$SecondCodon,sep = '_')
COMPAR<-function(x)	{
  cod1 <- unlist(strsplit(x,'_'))[1]
  cod2 <- unlist(strsplit(x,'_'))[2]
  #table <- c("FirstCOD", "SecondCOD")
  NumSub <- 0
  unlist(strsplit(cod1,''))[1]
  if (unlist(strsplit(cod1,''))[1] != unlist(strsplit(cod2,''))[1]){
    NumSub = 1 + NumSub
    FirstC <- unlist(strsplit(cod1,''))[1]
    SecC <- unlist(strsplit(cod2,''))[1]
  } 
  if (unlist(strsplit(cod1,''))[2] != unlist(strsplit(cod2,''))[2]){
    NumSub = 1 + NumSub
    FirstC <- unlist(strsplit(cod1,''))[2]
    SecC <- unlist(strsplit(cod2,''))[2]
  }
  if (unlist(strsplit(cod1,''))[3] != unlist(strsplit(cod2,''))[3])
  {
    NumSub = 1 + NumSub
    FirstC <- unlist(strsplit(cod1,''))[3]
    SecC <- unlist(strsplit(cod2,''))[3]
  }
  if (NumSub == 1) {nucleotides <- paste(FirstC,SecC,sep='_')}
  if (NumSub > 1)  {nucleotides  = 'MoreThanOne_SUBST'}
  return(nucleotides);
}

data$Subst = apply(as.matrix(data$FirstSecond),1,COMPAR)  
table$Subs = data$Subst
table = table[table$Subs != 'MoreThanOne_SUBST',]
data = data[data$Subst != 'MoreThanOne_SUBST',]

setwd(wd_out)
write.table(table, file="Mutational_spectra_in_Chordata.txt", quote = FALSE, row.names = FALSE)
