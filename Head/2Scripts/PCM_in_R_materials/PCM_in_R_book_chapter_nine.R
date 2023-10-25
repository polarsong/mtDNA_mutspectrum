rm(list = ls(all=TRUE))
## load the phytools package
library(phytools)
## simulate a pure-birth (Yule) tree using pbtree
tree<-pbtree(n=12,scale=100)
## split our plotting area in two
par(mfrow=c(2,1))
## graph our phylogeny
plotTree(tree,ftype="off",mar=c(4.1,4.1,2.1,1.1))
## compute the lineages through time using ltt
obj<-ltt(tree,plot=FALSE)
## draw vertical lines at each lineage accumulation
## event
abline(v=obj$times,lty="dotted",
       col=make.transparent("blue",0.5))
## add a horizontal axis and plot label
axis(1,cex.axis=0.8)
mtext("(a)",line=1,at=-10)
## create a second plot graphing our LTT
plot(obj,mar=c(5.1,4.1,2.1,1.1),bty="n",
     log.lineages=FALSE,las=1,cex.axis=0.8)
## add the same vertical lines as in panel a)
abline(v=obj$times,lty="dotted",
       col=make.transparent("blue",0.5))
## label our plot
mtext("(b)",line=1,at=-10)

darter.tree<-read.tree("etheostoma_percina_chrono.txt")
## plot our tree in fan style
plotTree(darter.tree,ftype="i",
         fsize=0.4,type="fan",lwd=1,part=0.88)
## compute the total height of the tree
h<-max(nodeHeights(darter.tree))
## graph a temporal axis without labeling
obj<-axis(1,pos=-2,at=h-c(0,5,10,15,20),
          cex.axis=0.5,labels=FALSE)
## add labels, but going backwards from the
## present day
text(obj,rep(-5,length(obj)),h-obj,
     cex=0.6)
## add a text label to the axis
text(mean(obj),-8,"time (mybp)",
     cex=0.8)

## compute "ltt" object
darter.ltt<-ltt(darter.tree, plot=FALSE)
## modify the figure margins
par(mar=c(5.1,4.1,2.1,2.1))
## plot "ltt" object
plot(darter.ltt,log.lineages=FALSE,log="y",
     col="blue",lwd=2,bty="n",las=1,
     cex.axis=0.8)
print(darter.ltt)

## resolve polytomies using multi2di
darter.tree<-multi2di(darter.tree)
## recompute "ltt" object
darter.ltt<-ltt(darter.tree,plot=FALSE)
darter.ltt

darter.mccr<-mccr(darter.ltt,rho=201/216,
                  nsim=500)
darter.mccr
par(mar=c(5.1,4.1,2.1,2.1))
plot(darter.mccr,main="",las=1,cex.axis=0.8)
