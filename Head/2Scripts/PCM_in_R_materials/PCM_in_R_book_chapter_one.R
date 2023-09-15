rm(list = ls(all=TRUE))
#To go threw this script - read book Phylogenetic Comparative Methods in R Liam J. Revell and Luke J. Harmon
library(ape)
library(phangorn)
library(phytools)
library(geiger)
#Chapter 1
text.string<-
  "(((((Robin,Iguana),((((Cow,Whale),Pig),Bat),
(Lemur,Human))),Coelacanth),Goldfish),Shark);"
vert.tree<-read.tree(text=text.string)
plot(vert.tree,no.margin=TRUE)
#Cool visualization
par(mfrow=c(2,2),mar=c(1.1,1.1,3.1,1.1))
plot(vert.tree)
mtext("(a)",line=1,adj=0)
plot(vert.tree,type="cladogram")
mtext("(b)",line=1,adj=0)
plot(unroot(vert.tree),type="unrooted",
     lab4ut="axial",x.lim=c(-2,6.5),
     y.lim=c(-3,7.5))
mtext("(c)",line=1,adj=0)
#Nodes
plotTree(vert.tree,offset=1,type="cladogram")
labelnodes(1:(Ntip(vert.tree)+vert.tree$Nnode),
           1:(Ntip(vert.tree)+vert.tree$Nnode),
           interactive=FALSE,cex=0.8)
#reading tree files
anolis.tree<-read.tree(file="Anolis.tre.txt")
anolis.tree
plotTree(anolis.tree,ftype="i",fsize=0.4,lwd=1)
Ntip(anolis.tree)
#writing tree files
write.tree(vert.tree,file="example.tre")
#Add arrows to a tree (pruning - extracting smaller tree)
pr.species<-c("cooki","poncensis",
              "gundlachi","pulchellus","stratulus",
              "krugi","evermanni","occultus","cuvieri",
              "cristatellus")
nodes<-sapply(pr.species,grep,x=anolis.tree$tip.label)
nodes
plotTree(anolis.tree,type="fan",fsize=0.6,lwd=1,
         ftype="i")
add.arrow(anolis.tree,tip=nodes,arrl=0.15,col="red",
          offset=2)
anolis.noPR<-drop.tip(anolis.tree,pr.species)
plotTree(anolis.noPR,type="fan",fsize=0.6,lwd=1,
         ftype="i")
#Extracting a clade
#MRCA - most recent common ancestor
node<-getMRCA(anolis.tree,pr.species[
  -which(pr.species%in%c("cuvieri","occultus"))])
node
#Also painting tree
plot(paintSubTree(anolis.tree,node,"b","a"), 
  type="fan",fsize=0.6,lwd=2,
  colors=setNames
  (
    c("gray","blue"),
    c("a","b")),
  ftype="i")
arc.cladelabels(anolis.tree,"clade to extract",node,
                1.2,1.25,mark.node=FALSE,cex=0.6)
pr.clade<-extract.clade(anolis.tree,node)
pr.clade
pr.tree<-keep.tip(anolis.tree,pr.species)
pr.tree
par(mfrow=c(1,2))
plotTree(pr.clade,ftype="i",mar=c(1.1,1.1,3.1,1.1),
         cex=1.1) #specific clad tree
mtext("(a)",line=0,adj=0)
plotTree(pr.tree,ftype="i",mar=c(1.1,1.1,3.1,1.1),
         cex=1.1) #pruning tree
mtext("(b)",line=0,adj=0)
#Interactive pruning
anolis.pruned<-collapseTree(anolis.tree)
#multyphylo
anolis.trees<-c(anolis.tree,anolis.noPR,pr.clade,
                pr.tree)
print(anolis.trees,details=TRUE)
#csv tree
anole.data<-read.csv(file="anole.data.csv",row.names=1,
                     header=TRUE)
ecomorph<-read.csv(file="ecomorph.csv",row.names=1,
                   header=TRUE,stringsAsFactors=TRUE)
#comparing files then pruning
library(geiger)
name.check(anolis.tree,anole.data)
chk<-name.check(anolis.tree,ecomorph)
chk
ecomorph.tree<-drop.tip(anolis.tree,chk$tree_not_data)
ecomorph.tree
ecomorph.data<-anole.data[ecomorph.tree$tip.label,]
head(ecomorph.data)
name.check(ecomorph.tree,ecomorph.data)
#phylopca!!!
ecomorph.pca<-phyl.pca(ecomorph.tree,ecomorph.data)
ecomorph.pca
par(mar=c(4.1,4.1,2.1,1.1),las=1) ## set margins
plot(ecomorph.pca,main="")
ecomorph.pca$Evec[,1]<--ecomorph.pca$Evec[,1]
ecomorph.pca$L[,1]<--ecomorph.pca$L[,1]
ecomorph.pca$S<-scores(ecomorph.pca,
                       newdata=ecomorph.data)

par(cex.axis=0.8,mar=c(5.1,5.1,1.1,1.1))
phylomorphospace(ecomorph.tree,
                 scores(ecomorph.pca)[,1:2],
                 ftype="off",node.size=c(0,1),bty="n",las=1,
                 xlab="PC1 (overall size)",
                 ylab=expression(paste("PC2 ("%up%"lamellae number, "
                                       %down%"tail length)")))
eco<-setNames(ecomorph[,1],rownames(ecomorph))
ECO<-to.matrix(eco,levels(eco))
tiplabels(pie=ECO[ecomorph.tree$tip.label,],cex=0.5)
legend(x="bottomright",legend=levels(eco),cex=0.8,pch=21,
       pt.bg=rainbow(n=length(levels(eco))),pt.cex=1.5)
