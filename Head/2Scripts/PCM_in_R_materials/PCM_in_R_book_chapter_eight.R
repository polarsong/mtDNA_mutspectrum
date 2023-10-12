rm(list = ls(all=TRUE))
## load libraries
library(phytools)
## read tree from file
eel.tree<-read.tree("elopomorph.txt")
print(eel.tree,printlen=2)
## read data
eel.data<-read.csv("elopomorph.csv",row.names=1,
                   stringsAsFactors=TRUE)
head(eel.data)
## extract total body length and log-transform
lnTL<-setNames(log(eel.data$Max_TL_cm),rownames(eel.data))
head(lnTL)
## estimate ancestral states using fastAnc
fit.lnTL<-fastAnc(eel.tree,lnTL,vars=TRUE,CI=TRUE)
print(fit.lnTL,printlen=10)
## plot eel phylogeny using plotTree
plotTree(eel.tree,ftype="i",fsize=0.5,lwd=1)
## add node labels for reference
labelnodes(1:eel.tree$Nnode+Ntip(eel.tree),
           1:eel.tree$Nnode+Ntip(eel.tree),
           interactive=FALSE,cex=0.5)
## compute "contMap" object
eel.contMap<-contMap(eel.tree,lnTL,
                     plot=FALSE,lims=c(2.7,5.8))
## change the color gradient to a custom gradient
eel.contMap<-setMap(eel.contMap,
                    c("white","orange","black"))
## plot "contMap" object
plot(eel.contMap,sig=2,fsize=c(0.4,0.7),
     lwd=c(2,3),leg.txt="log(total length cm)")
## identify the tips descended from node 102
tips<-extract.clade(eel.tree,102)$tip.label
tips
## prune "contMap" object to retain only these tips
pruned.contMap<-keep.tip.contMap(eel.contMap,tips)
## plot object
plot(pruned.contMap,xlim=c(-2,90),lwd=c(3,4),
     fsize=c(0.7,0.8))
## add error bars
errorbar.contMap(pruned.contMap,lwd=8)
