rm(list = ls(all=TRUE))
# set values for time steps and sigma squared parameter
t<-0:100
sig2<-0.01
## simulate a set of random changes
x<-rnorm(n=length(t)-1,sd=sqrt(sig2))
## compute their cumulative sum
x<-c(0,cumsum(x))
# create a plot with nice margins
par(mar=c(5.1,4.1,2.1,2.1))
plot(t,x,type="l",ylim=c(-2,2),bty="n",
     xlab="time",ylab="trait value",las=1,
     cex.axis=0.8)

# set number of simulations
nsim<-100
# create matrix of random normal deviates
X<-matrix(rnorm(n=nsim*(length(t)-1),sd=sqrt(sig2)),
          nsim,length(t)-1)
# calculate the cumulative sum of these deviates
# this is now a simulation of Brownian motion
X<-cbind(rep(0,nsim),t(apply(X,1,cumsum)))
# plot the first one
par(mar=c(5.1,4.1,2.1,2.1))
plot(t,X[1,],ylim=c(-2,2),type="l",bty="n",
     xlab="time",ylab="trait value",las=1,
     cex.axis=0.8)
# plot the rest
invisible(apply(X[2:nsim,],1,function(x,t) lines(t,x),
                t=t))

# create matrix of random normal deviates
# but with a smaller sd
X<-matrix(rnorm(n=nsim*(length(t)-1),sd=sqrt(sig2/10)),
          nsim,length(t)-1)
# calculate the cumulative sum of these changes
# this is now a simulation of Brownian motion
X<-cbind(rep(0,nsim),t(apply(X,1,cumsum)))
# plot as above
par(mar=c(5.1,4.1,2.1,2.1))
plot(t,X[1,],ylim=c(-2,2),type="l",bty="n",
     xlab="time",ylab="trait value",las=1,
     cex.axis=0.8)
invisible(apply(X[2:nsim,],1,function(x,t) lines(t,x),
                t=t))

# calculate variance of columns
v<-apply(X,2,var)
# plot the results
par(mar=c(5.1,4.1,2.1,2.1))
plot(t,v,ylim=c(0,0.1),type="l",xlab="time",
     ylab="variance",bty="n",las=1,
     cex.axis=0.8)
lines(t,t*sig2/10,lwd=3,col=rgb(0,0,0,0.1))
legend("topleft",c("observed variance","expected variance"),
       lwd=c(1,3),col=c("black",rgb(0,0,0,0.1)),
       bty="n",cex=0.8)
# find variance at the end of the simulations
var(X[,length(t)])

## load phytools package
library(phytools)
## simulate a tree and Brownian evolution on that
## tree using simBMphylo
object<-simBMphylo(n=6,t=100,sig2=0.01,
                   fsize=0.8,cex.axis=0.6,cex.lab=0.8,
                   las=1)
## pull the phylogeny out of the object we simulated
## for figure 4.5 using simBMphylo
tree<-object$tree
## simulate 1000 instance of Brownian evolution on that
## tree
X<-fastBM(tree,nsim=1000)
## set the orientation of the axis labels to be
## horizontal
par(las=1)
## create a scatterplot matrix from our simulated
## data using pairs
pairs(t(X)[,tree$tip.label[6:1]],pch=19,
      col=make.transparent("blue",0.05),
      cex.axis=0.9)

## read bacterial data from file
bacteria.data<-read.csv("bac_rates.csv", row.names=1)
head(bacteria.data,3)
bacteria.tree<-read.tree("bac_rates.txt")
print(bacteria.tree,printlen=2)
## graph phylogeny using plotTree
plotTree(bacteria.tree,ftype="i",fsize=0.5,
         lwd=1,mar=c(2.1,2.1,0.1,1.1))
## add a horizontal axis to our plot
axis(1,at=seq(0,1,length.out=5),cex.axis=0.8)

library(geiger)
name.check(bacteria.tree,bacteria.data)
genome_size<-bacteria.data[,"Genome_Size_Mb"]
genome_size
names(genome_size)<-rownames(bacteria.data)
head(genome_size)

## fit Brownian motion model using fitContinuous
fitBM_gs<-fitContinuous(bacteria.tree,genome_size)
fitBM_gs

## pull our mutation accumulation rate as a named vector
mutation<-setNames(bacteria.data[,"Accumulation_Rate"],
                   rownames(bacteria.data))
head(mutation)

## set up for side-by-side plots
par(mfrow=c(1,2),mar=c(6.1,4.1,2.1,1.1))
## histogram of mutation accumulation rates on original scale
hist(mutation,main="",las=2,xlab="",
     cex.axis=0.7,cex.lab=0.9,
     breaks=seq(min(mutation),max(mutation),
                length.out=12))
mtext("(a)",adj=0,line=1)
mtext("rate",side=1,line=4,cex=0.9)
## histogram of mutation accumulation rates on log scale
ln_mutation<-log(mutation)
hist(ln_mutation,main="",las=2,xlab="",
     cex.axis=0.7,cex.lab=0.9,
     breaks=seq(min(ln_mutation),max(ln_mutation),
                length.out=12))
mtext("(b)",adj=0,line=1)
mtext("ln(rate)",side=1,line=4,cex=0.9)
## fit Brownian motion model to log(mutation accumulation)
fitBM_ar<-fitContinuous(bacteria.tree,ln_mutation)
fitBM_ar
phylosig(bacteria.tree, genome_size)

## test for significant phylogenetic signal using
## Blomberg’s K
K_gs<-phylosig(bacteria.tree,genome_size,
               test=TRUE,nsim=10000)
K_gs
## set plot margins and font size
par(cex=0.8,mar=c(5.1,4.1,2.1,2.1))
## plot null-distribution and observed value of K
plot(K_gs,las=1,cex.axis=0.9)

## test for phylogenetic signal in mutation accumulation
## rate
K_ar<-phylosig(bacteria.tree,ln_mutation,
               test=TRUE,nsim=10000)
K_ar
## plot the results
par(cex=0.8,mar=c(5.1,4.1,2.1,2.1))
plot(K_ar,las=1,cex.axis=0.9)

## simulate 10000 datasets
nullX<-fastBM(bacteria.tree,nsim=10000)
## for each, carry out a test for phylogenetic signal
## and accumulate these into a vector using sapply
nullK<-apply(nullX,2,phylosig,tree=bacteria.tree)
## calculate P-values
Pval_gs<-mean(nullK<=K_gs$K)
Pval_gs
## 
Pval_ar<-mean(nullK<=K_ar$K)
Pval_ar

## set up for side-by-side plots
par(mfrow=c(1,2))
## plot for Genome size
## null distribution
hist(c(nullK,K_gs$K),breaks=30,col="lightgray",
     border="lightgray",main="",xlab="K",las=1,
     cex.axis=0.7,cex.lab=0.9,ylim=c(0,4000))
## actual value as an arrow
arrows(x0=K_gs$K,y0=par()$usr[4],y1=0,length=0.12,
       col=make.transparent("blue",0.5),lwd=2)
text(K_gs$K,0.96*par()$usr[4],
     paste("observed value of K (P = ",
           round(Pval_gs,4),")",sep=""),
     pos=4,cex=0.8)
mtext("(a)",line=1,adj=0)
## plot for mutation accumulation rate
## null distribution
hist(c(nullK,K_ar$K),breaks=30,col="lightgray",
     border="lightgray",main="",xlab="K",las=1,
     cex.axis=0.7,cex.lab=0.9,ylim=c(0,4000))
## actual value as an arrow
arrows(x0=K_ar$K,y0=par()$usr[4],y1=0,length=0.12,
       col=make.transparent("blue",0.5),lwd=2)
text(K_ar$K,0.96*par()$usr[4],
     paste("observed value of K (P = ",
           round(Pval_ar,4),")",sep=""),
     pos=4,cex=0.8)
mtext("(b)",line=1,adj=0)

## compute phylogenetic signal, lambda, for genome size
## and mutation accumulation rate
phylosig(bacteria.tree,genome_size,method="lambda")
phylosig(bacteria.tree,ln_mutation,method="lambda")

## test for significant phylogenetic signal, lambda,
## in each of our two traits
lambda_gs<-phylosig(bacteria.tree,genome_size,
                    method="lambda",test=TRUE)
lambda_gs

lambda_ar<-phylosig(bacteria.tree, ln_mutation,
  method="lambda",test=TRUE)
lambda_ar

## our plot area into 1 column and two rows
par(mfrow=c(2,1),mar=c(5.1,4.1,2.1,2.1),
    cex=0.8)
## plot the likelihood surfaces of lambda for each of our
## two traits
plot(lambda_gs,las=1,cex.axis=0.9,bty="n",
     xlim=c(0,1.1))
mtext("(a)",line=1,adj=0)
plot(lambda_ar,las=1,cex.axis=0.9,bty="n",
     xlim=c(0,1.1))
mtext("(b)",line=1,adj=0)