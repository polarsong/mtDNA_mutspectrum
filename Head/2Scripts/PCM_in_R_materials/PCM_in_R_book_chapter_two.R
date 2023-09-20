rm(list = ls(all=TRUE))
## load packages
library(phytools)
## read in tree
tree<-read.tree(
  text="((A,B,C,D,E,F,G,H,I,J,K,L,M),
(N,O,P,Q,R,S,T,U,V,W,X,Y,Z));")
## set branch lengths on the tree
tree<-compute.brlen(tree,power=1.8)
## simulate data, independently for x & y
x<-fastBM(tree)
y<-fastBM(tree)
## plot the results with clades A & B labeled
## split plotting area
par(mfrow=c(1,2))
## graph tree
plotTree(tree,type="cladogram",ftype="off",
         mar=c(5.1,4.1,3.1,2.1),color="darkgray",
         xlim=c(0,1.3),ylim=c(1,Ntip(tree)))
## add points at the tips of the tree to match those
## on our scatterplot
points(rep(1,13),1:13,pch=21,bg="lightgray",
       cex=1.2)
points(rep(1,13),14:26,pch=22,bg="black",cex=1.2)
## add clade labels to the tree
cladelabels(tree,"A",node=28,offset=2)
cladelabels(tree,"B",node=29,offset=2)
mtext("(a)",line=1,adj=0,cex=1.5)
## create scatterplot of x & y
par(mar=c(5.1,4.1,3.1,2.1))
plot(x,y,bty="n",las=1)
points(x[1:13],y[1:13],pch=21,bg="lightgray",
       cex=1.2)
points(x[14:26],y[14:26],pch=22,bg="black",cex=1.2)
mtext("(b)",line=1,adj=0,cex=1.5)
#PIC
mammalHR<-read.csv("mammalHR.csv",row.names=1)
## set margins of the plot
par(mar=c(5.1,5.1,1.1,1.1))
## create scatterplot
plot(homeRange~bodyMass,data=mammalHR,
     xlab="body mass (kg)",
     ylab=expression(paste("home range (km"^"2",")")),
     pch=21,bg="gray",cex=1.2,log="xy",las=1,cex.axis=0.7,
     cex.lab=0.9,bty="n")
#OLS
fit.ols<-lm(log(homeRange)~log(bodyMass),data=mammalHR)
fit.ols
summary(fit.ols)
## set margins and graph scatterplot
par(mar=c(5.1,5.1,1.1,1.1))
plot(homeRange~bodyMass,data=mammalHR,
     xlab="body mass (kg)",
     ylab=expression(paste("home range (km"^"2",")")),
     pch=21,bg="gray",cex=1.2,log="xy",las=1,
     cex.axis=0.7,cex.lab=0.9,bty="n")
## add the line of best fit from lm
lines(mammalHR$bodyMass,exp(predict(fit.ols)),lwd=2,
      col="darkgray")
#PIC
mammal.tree<-read.tree("../../2Scripts/PCM_in_R_materials/mammalHR.phy.txt")
## plot phylogeny of mammals
plotTree(mammal.tree,ftype="i",fsize=0.7,lwd=1)
## add node labels to the plotted tree
nodelabels(bg="white",cex=0.5,frame="circle")
## pull our home range and body mass as
## numeric vectors
homeRange<-setNames(mammalHR[,"homeRange"],
                    rownames(mammalHR))
bodyMass<-setNames(mammalHR[,"bodyMass"],
                   rownames(mammalHR))
## compute PICs for home range and body size
pic.homerange<-pic(log(homeRange),mammal.tree)
pic.bodymass<-pic(log(bodyMass),mammal.tree)
head(pic.homerange,n=20)
## fit linear model to PICs without intercept
fit.pic<-lm(pic.homerange~pic.bodymass+0)
fit.pic
summary(fit.pic)
## set margins
par(mar=c(5.1,5.1,1.1,1.1))
## graph scatterplot of contrasts
plot(pic.homerange~pic.bodymass,
     xlab="PICs for log(body mass)",
     ylab="PICs for log(range size)",
     pch=21,bg="gray",cex=1.2,las=1,
     cex.axis=0.7,cex.lab=0.9,bty="n")
## add gridlines to the plot
abline(h=0,lty="dotted")
abline(v=0,lty="dotted")
## reset graphing limits of the plot to the
## x/y range of our PICs
clip(min(pic.bodymass),max(pic.bodymass),
     min(pic.homerange),max(pic.homerange))
## graph our fitted line
abline(fit.pic,lwd=2,col="darkgray")
