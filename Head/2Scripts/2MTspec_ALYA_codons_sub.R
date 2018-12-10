detach("package:Biostrings", unload=TRUE)
library(seqinr)
rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

wd <- getwd()
wd_in = gsub('Head/2Scripts','Body/2Derived/POLARIZEDBR_DATA_ML', wd)  #SWITCH to MP-prefix if MPanalysis data
wd_out = gsub('Head/2Scripts','Body/2Derived/TOTAL_SUBS_ML', wd) #SWITCH to MP-prefix if MPanalysis data
dir.create (wd_out)
setwd(wd_in)
wd_in




polardata <- list.files(pattern=".*\\.txt")

for (pfiles in polardata){
#pfiles = "Bubalus_bubalis.CytB.POLARISED.txt"
  polartable <- read.table(pfiles, header = TRUE)
  for (y in 1:nrow(polartable)){
    mdns <- as.character(polartable$MoreDeepNodeSeq[y])
    msns <- as.character(polartable$MoreShallowNodeSeq[y])
    s1 <- s2c(mdns)       #conversion of a string into a vector of chars
    s2 <-  s2c(msns)
    cdns1 <- splitseq(s1, 0, 3) #split a sequence into sub-sequences (codons)
    cdns2 <- splitseq(s2, 0, 3)
    
    codonposition <-c(1:length(cdns1))
    
    final <- c()
    allframes <- c()
    for (x in codonposition){
      if (cdns1[x] != cdns2[x]){
        taa1 <- translate(s2c(cdns1[x]), frame = 0, sens = "F", numcode = 2, NAstring = "X", ambiguous = FALSE)
        taa2 <- translate(s2c(cdns2[x]), frame = 0, sens = "F", numcode = 2, NAstring = "X", ambiguous = FALSE)
        el <- strsplit(pfiles, "\\.")
        species <- el [[1]][1]
        gene <- el[[1]][2]
        node1 <- as.character(polartable$MoreDeepNode[y])
        node2 <- as.character(polartable$MoreShallowNode[y])
        branch <- as.character(polartable$BranchPosition[y])
        if (x > 1 & x < length(codonposition))
        {
          final <- rbind(final, c(species, gene, node1, node2, branch, x, cdns1[x-1], cdns1[x], cdns1[x+1], cdns2[x-1], cdns2[x], cdns2[x+1], taa1, taa2))
        }
        if (x == 1)
        {
          final <- rbind(final, c(species, gene, node1, node2, branch, x, '---', cdns1[x], cdns1[x+1], '---', cdns2[x], cdns2[x+1], taa1, taa2))
        }
        if (x == length(codonposition))
        {
          final <- rbind(final, c(species, gene, node1, node2, branch, x, cdns1[x-1], cdns1[x], '---', cdns2[x-1], cdns2[x], '---', taa1, taa2))
        }
        da <- data.frame(final)
        colnames(da) <- c("Species", "Gene", "AncestralSeqName", "DescendantSeqName", "Branch", "CodonPosition", "PreviousAncCodon", "AncestorCodon", "NextAncCodon", "PreviousDesCodon", "DescendantCodon", "NextDesCodon", "AncestralAA", "DescendantAA")
        str(da)
        setwd(wd_out) 
        write.table(da, file=paste(gsub("POLARISED.txt", "", pfiles), node1, ".", node2, ".SUBS.txt", sep = ""), quote = FALSE, row.names = FALSE)
        setwd(wd_in)
        }
      }
  }
}


