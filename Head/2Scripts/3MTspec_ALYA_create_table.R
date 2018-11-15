rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

#USER = 'KOSTYA';
USER = 'ALYA';
#USER = 'TOY'

WDA_IN <- normalizePath("C:\\Users\\polar\\Documents\\BODY\\POLARIZEDBR_DATA\\CodonTable", winslash = "\\")
if (USER == 'KOSTYA') {setwd('/home/kostya/konstantin/SCIENCE_PROJECTS_HEAD/MITOCHONDRIA/MutSpectrum/Results/re/');}  
if (USER == 'ALYA') {setwd (WDA_IN);} 
if (USER == 'TOY') {setwd (WDA_TOY);}


filevector <- list.files(pattern=".*\\.SUBS\\.txt")


for (i in 1:length(filevector)){
  temp <- read.table(as.character(filevector[i]), header = TRUE, colClasses=c('character'))
  bob <- data.frame(lapply(temp, as.character), stringsAsFactors=FALSE)
  if (i == 1){
    final <- bob
  }
  if (i > 1){
    final <- rbind(final, bob)
  }
}
write.table(final, file="1.txt", quote = FALSE, row.names = FALSE)



str(temp)

table(final$AncestralAA)
err3 = final[final$AncestralAA == '*',]
table(final$DescendantAA)
table(final$AncestorCodon)
table(final$AncestorCodon)

