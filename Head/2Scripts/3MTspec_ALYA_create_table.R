rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

#USER = 'KOSTYA';
USER = 'ALYA';
#USER = 'TOY'

wd <- getwd()
wd_in = gsub('Head/2Scripts','Body/2Derived/TOTAL_SUBS_ML', wd)  #SWITCH to MP-prefix if MPanalysis data
wd_out = gsub('Head/2Scripts','Body/2Derived/', wd) 
setwd(wd_in)
wd_in


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
setwd(wd_out)
write.table(final, file="1.txt", quote = FALSE, row.names = FALSE)



str(temp)

table(final$AncestralAA)
err3 = final[final$AncestralAA == '*',]
table(final$DescendantAA)
table(final$AncestorCodon)
table(final$AncestorCodon)

