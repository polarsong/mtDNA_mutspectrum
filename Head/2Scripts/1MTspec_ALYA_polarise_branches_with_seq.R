source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")

###################################
###### 1: Derive From MAtrix File paired node names: MoreDeepNode, MoreShallowNode
###################################

library("Biostrings")
library("seqinr")
rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 


wd <- getwd()
wd = gsub('Head/2Scripts','Body/1Raw/pipeline_trees', wd)
setwd(wd)
wd

### READ THREE FILES FOR EACH GENE-SPECIES:

filelist <- list.files(pattern=".*\\.MLancestors\\.fa")           #SWITCH to MP-prefix if MPanalysis data
filelist = gsub('.MLancestors.fa','',filelist)                   #SWITCH to MP-prefix if MPanalysis data
for (f in filelist) {
  #f = 'Abbottina_rivularis.ATP6';
  matrixname <- paste(f, ".tree.matrix", sep = "")
  Matrix <- read.table(matrixname, header=TRUE, sep="\t")
  
  terminalname <- paste(f, ".terminals.nuc.fa", sep = "")
  Terminal <- readDNAStringSet(terminalname)
  Terminal = data.frame(Terminal); names(Terminal) = c('TerminalSequence')
  Terminal$TerminalNodeNumber = rownames(Terminal)
  str(Terminal)
  
  ancestralname <- paste(f, ".ancestors.fa", sep = "")
  Ancestral <- read.table(ancestralname, header = FALSE, sep = ' '); names(Ancestral)=c('InternalNodeNumber', 'InternalSequence')
  Ancestral$InternalNodeNumber = as.character(Ancestral$InternalNodeNumber)
  Ancestral$InternalSequence = as.character(Ancestral$InternalSequence)
  str(Ancestral)
  
  names(Terminal)=c('Seq','Name')
  names(Ancestral)=c('Name','Seq')
  Seq = rbind(Ancestral,Terminal)
  
  
  ### FIND ALL PAIRS OF ANCESTOR-DERIVED BRANCHES AND SAVE THEM AS MoreDeepNode AND MoreShallowNode
  
  Final = data.frame('PleaseDeleteMe','PleaseDeleteMe'); names(Final) = c('MoreDeepNode','MoreShallowNode')
  row.names(Matrix) = Matrix$X  # name each raw as first column 
  Matrix = Matrix[-1]           # delete first column
  OneBeforeLast = nrow(Matrix); 
  TempMatrix = Matrix;
  TempMatrix = TempMatrix[-OneBeforeLast,]; TempMatrix = TempMatrix[,-OneBeforeLast]
  Max = nrow(TempMatrix); 
  # TempMatrix is square matrix which we cut down step by step: 
  for (OneBeforeLast in Max:1)
  { # OneBeforeLast = 6
    LINE = TempMatrix[OneBeforeLast,];
    for (j in 1:ncol(LINE))
    {#  j = 2
      if (LINE[j] == 1) {OneLine = data.frame(rownames(LINE),colnames(LINE)[j]); names(OneLine) = c('MoreDeepNode','MoreShallowNode'); Final = rbind(Final,OneLine);}
    }
    TempMatrix = TempMatrix[-OneBeforeLast,]; TempMatrix = TempMatrix[,-OneBeforeLast];
    IfOnlyTerminalBranchesLeft = TempMatrix[grep("RN_|OUTGRP",row.names(TempMatrix), invert=TRUE),];  # OUTGRP'
    if (nrow(IfOnlyTerminalBranchesLeft) == 0) 
    {break;}
  }
  Final = Final[Final$MoreShallowNode != 'PleaseDeleteMe',]
  Final$MoreShallowNode = gsub('X','',Final$MoreShallowNode)  # delete 'X'
  Final$MoreDeepNode = gsub('X','',Final$MoreDeepNode)        # delete 'X'
  
  Final$BranchPosition = 'Internal';
  for (i in 1:nrow(Final))
  { # i = 1
    if (length(as.character(grep('RN_',Final$MoreShallowNode[i])))) {Final$BranchPosition[i] = 'External'}
  }
  
  ####CREATE TABLE WITH ANCESTOR-DERIVED NAMES AND SEQUENSES
  
  Final$MoreDeepNodeSeq = 'NN';
  Final$MoreShallowNodeSeq = 'NN';
  for (i in 1:nrow(Final))
  { # i = 1
    TempDeepNode = Final$MoreDeepNode[i]
    Final$MoreDeepNodeSeq[i] = Seq[Seq$Name == TempDeepNode,]$Seq
    TempShallowNode = Final$MoreShallowNode[i]
    Final$MoreShallowNodeSeq[i] = Seq[Seq$Name == TempShallowNode,]$Seq
    
  }  
  
  if (USER == 'ALYA') {setwd (WDA_OUT);} 
  write.table(Final, file=paste(f, ".POLARISED.txt", sep = ""), quote = FALSE, row.names = FALSE)
  
  if (USER == 'ALYA') {setwd (WDA_IN);} 
}


