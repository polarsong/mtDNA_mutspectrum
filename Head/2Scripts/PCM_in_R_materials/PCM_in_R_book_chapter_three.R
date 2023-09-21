rm(list = ls(all=TRUE))
## load packages
library(phytools)
## read data from file
primate.data<-read.csv("primateEyes.csv",row.names=1,
                       stringsAsFactors=TRUE)
## inspect data
head(primate.data,4)
## read tree from file and inspect
primate.tree<-read.tree("primateEyes.phy.txt")
print(primate.tree,printlen=2)
## extract orbit area from our data frame and add names
orbit.area<-setNames(primate.data[,"Orbit_area"],
                     rownames(primate.data))
## extract skull length from our data frame and add names
skull.length<-setNames(primate.data[,"Skull_length"],
                       rownames(primate.data))
## compute PICs on the log-transformed values of both traits
pic.orbit.area<-pic(log(orbit.area),primate.tree)
pic.skull.length<-pic(log(skull.length),
                      primate.tree)
## fit a linear regression to orbit area as a function of
## skull length, without an intercept term
pic.primate<-lm(pic.orbit.area~pic.skull.length+0)
summary(pic.primate)
## set plotting parameters
par(mfrow=c(1,2),
    mar=c(5.1,4.6,2.1,1.1))
## plot our raw data in the original space
plot(orbit.area~skull.length,log="xy",
     pch=21,bg=palette()[4],cex=1.2,
     bty="n",xlab="skull length (cm)",
     ylab=expression(paste("orbit area (",mm^2,")")),
     cex.lab=0.8,cex.axis=0.7,las=1)
mtext("(a)",line=0,adj=0,cex=0.8)
## plot our phylogenetic contrasts
plot(pic.orbit.area~pic.skull.length,pch=21,
     bg=palette()[4],cex=1.2,
     bty="n",xlab="PICs for log(skull length)",
     ylab="PICs for log(orbit area)",
     cex.lab=0.8,cex.axis=0.7,las=1)
mtext("(b)",line=0,adj=0,cex=0.8)
## limit the plotting area to the range of our two traits
clip(min(pic.skull.length),max(pic.skull.length),
     min(pic.orbit.area),max(pic.orbit.area))
## add our fitted contrasts regression line
abline(pic.primate,lwd=2)
## set plotting parameters
par(mfrow=c(1,2),
mar=c(5.1,4.6,2.1,1.1))
## plot our raw data in the original space
plot(orbit.area~skull.length,log="xy",
pch=21,bg=palette()[4],cex=1.2,
bty="n",xlab="skull length (cm)",
ylab=expression(paste("orbit area (",mm^2,")")),
cex.lab=0.8,cex.axis=0.7,las=1)
mtext("(a)",line=0,adj=0,cex=0.8)
## plot our phylogenetic contrasts
plot(pic.orbit.area~pic.skull.length,pch=21,
bg=palette()[4],cex=1.2,
bty="n",xlab="PICs for log(skull length)",
ylab="PICs for log(orbit area)",
cex.lab=0.8,cex.axis=0.7,las=1)
mtext("(b)",line=0,adj=0,cex=0.8)
## limit the plotting area to the range of our two traits
clip(min(pic.skull.length),max(pic.skull.length),
min(pic.orbit.area),max(pic.orbit.area))
## add our fitted contrasts regression line
abline(pic.primate,lwd=2)
library(nlme)
spp<-rownames(primate.data)
corBM<-corBrownian(phy=primate.tree,form=~spp)
corBM
pgls.primate<-gls(log(Orbit_area)~log(Skull_length),
                  data=primate.data,correlation=corBM)
summary(pgls.primate)


coef(pic.primate)
coef(pgls.primate)
abs(coef(pic.primate)[1]-coef(pgls.primate)[2])

## set the random number generator seed
set.seed(88)
## simulate a random 5-taxon tree
tree<-pbtree(n=5,scale=10,tip.label=LETTERS[5:1])
## subdivide our plotting area into two panels
par(mfrow=c(2,1))
## plot the tree
plotTree(tree,mar=c(3.1,1.1,4.1,1.1),fsize=1.25,
         ylim=c(0.5,5.4))
## add a horizontal axis
axis(1)
## add edge labels giving the branch lengths
edgelabels(round(tree$edge.length,2),pos=3,
           frame="none",cex=0.9)
mtext("(a)",line=1,adj=0)
## switch to the second panel
plot.new()
## set new plot margins and plot dimensions
par(mar=c(3.1,1.1,4.1,1.1))
plot.window(xlim=c(0,6),ylim=c(0,6))
## add a grid of lines for our correlation matrix
lines(c(0,6,6,0,0),c(0,0,6,6,0))
for(i in 1:5) lines(c(i,i),c(0,6))
for(i in 1:5) lines(c(0,6),c(i,i))
## compute the assumed correlation structure
V<-cov2cor(vcv(tree)[LETTERS[1:5],LETTERS[1:5]])
## print it into the boxes of our grid
for(i in 1:5) text(i+0.5,5.5,LETTERS[i],cex=1.1)
for(i in 1:5) text(0.5,5.5-i,LETTERS[i],cex=1.1)
for(i in 1:5) for(j in 1:5) text(0.5+i,5.5-j,
                                 round(V[i,j],2),cex=1.1)
mtext("(b)",line=1,adj=0)


corLambda<-corPagel(value=1,phy=primate.tree,form=~spp)
corLambda
pgls.Lambda<-gls(log(Orbit_area)~log(Skull_length),
                 data=primate.data,correlation=corLambda)
summary(pgls.Lambda)
primate.ancova<-gls(log(Orbit_area)~log(Skull_length)+
                      Activity_pattern,data=primate.data,
                    correlation=corBM)
anova(primate.ancova)
## set the margins of our plot using par
par(mar=c(5.1,5.1,2.1,2.1))
## set the point colors for the different levels
## of our factor
pt.cols<-setNames(c("#87CEEB","#FAC358","black"),
                  levels(primate.data$Activity_pattern))
## plot the data
plot(Orbit_area~Skull_length,data=primate.data,pch=21,
     bg=pt.cols[primate.data$Activity_pattern],
     log="xy",bty="n",xlab="skull length (cm)",
     ylab=expression(paste("orbit area (",mm^2,")")),
     cex=1.2,cex.axis=0.7,cex.lab=0.8)
## add a legend
legend("bottomright",names(pt.cols),pch=21,pt.cex=1.2,
       pt.bg=pt.cols,cex=0.8)
## create a common set of x values to plot our
## different lines for each level of the factor
xx<-seq(min(primate.data$Skull_length),
        max(primate.data$Skull_length),length.out=100)
## add lines for each level of the factor
lines(xx,exp(predict(primate.ancova,
                     newdata=data.frame(Skull_length=xx,
                                        Activity_pattern=as.factor(rep("Cathemeral",100))))),
      lwd=2,col=pt.cols["Cathemeral"])
lines(xx,exp(predict(primate.ancova,
                     newdata=data.frame(Skull_length=xx,
                                        Activity_pattern=as.factor(rep("Diurnal",100))))),
      lwd=2,col=pt.cols["Diurnal"])
lines(xx,exp(predict(primate.ancova,
                     newdata=data.frame(Skull_length=xx,
                                        Activity_pattern=as.factor(rep("Nocturnal",100))))),
      lwd=2,col=pt.cols["Nocturnal"])
