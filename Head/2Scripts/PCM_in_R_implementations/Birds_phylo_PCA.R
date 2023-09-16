#birds phyloPCA
rm(list = ls(all=TRUE))
library(ape)
library(phangorn)
library(phytools)
library(geiger)
#Reading and plotting tree
birds_tree<-read.tree(file="anc_kg.treefile")
birds_tree
plotTree(birds_tree,ftype="i", type = 'fan', fsize=0.4,lwd=1)
#Preparing data
df_mtdna = read.csv('Birds_dataset_paper.csv')
df_eco = read.csv('flying_birds.csv')
df_names = as.data.frame(unique(df_mtdna$species_name))
names(df_names) = c('species_name')
df_eco = df_eco[,c(2,3,4)]
names(df_eco) = c('species_name', 'Flight','Diving')
df_eco = merge(df_names, df_eco)
df_eco$species_name = gsub(' ', '_', df_eco$species_name)
typeof(df_eco$species_name)
df_eco = df_eco[df_eco$species_name != "Agapornis_pullarius",]
df_eco = df_eco[df_eco$species_name != "Mergus_squamatus",]
row.names(df_eco) = df_eco$species_name

name.check(birds_tree,df_eco)
df_eco = na.omit(df_eco)
name.check(birds_tree,df_eco)
unique(df_eco$Flight)
df_eco0 = df_eco[df_eco$Flight == '0',]
df_eco1 = df_eco[df_eco$Flight == "Sphenisciformes",]
df_eco2 = df_eco[df_eco$Flight == "Apterygiformes",]
df_eco3 = df_eco[df_eco$Flight == "Gruiformes",]
df_eco4 = df_eco[df_eco$Flight == "Casuariiformes",]
df_eco5 = df_eco[df_eco$Flight == "Tinamiformes",]
df_eco6 = df_eco[df_eco$Flight == "Columbiformes",]
df_eco7 = df_eco[df_eco$Flight == "Rheiformes",]
df_eco8 = df_eco[df_eco$Flight == "Eurypygiformes",]
df_eco9 = df_eco[df_eco$Flight == "Psittaciformes",]
df_eco10 = df_eco[df_eco$Flight == "Struthioniformes",]

df_eco1$Flight = 1
df_eco2$Flight = 2
df_eco3$Flight = 3
df_eco4$Flight = 4
df_eco5$Flight = 5
df_eco6$Flight = 6
df_eco7$Flight = 7
df_eco8$Flight = 8
df_eco9$Flight = 9
df_eco10$Flight = 10

df_eco_cooler = rbind(df_eco0, df_eco1, df_eco2, df_eco3, df_eco4, df_eco5, df_eco6, df_eco7, df_eco8, df_eco9, df_eco10)

unique(df_eco_cooler$Diving)

df_eco0 = df_eco_cooler[df_eco_cooler$Diving == '0',]
df_eco1 = df_eco_cooler[df_eco_cooler$Diving == "Charadriiformes",]
df_eco2 = df_eco_cooler[df_eco_cooler$Diving == "Anseriformes",]
df_eco3 = df_eco_cooler[df_eco_cooler$Diving == "Coraciiformes",]
df_eco4 = df_eco_cooler[df_eco_cooler$Diving == "Suliformes",]
df_eco5 = df_eco_cooler[df_eco_cooler$Diving == "Passeriformes",]
df_eco6 = df_eco_cooler[df_eco_cooler$Diving == "Procellariiformes",]
df_eco7 = df_eco_cooler[df_eco_cooler$Diving == "Gruiformes",]
df_eco8 = df_eco_cooler[df_eco_cooler$Diving == "Gaviiformes",]
df_eco9 = df_eco_cooler[df_eco_cooler$Diving == "Podicipediformes",]
df_eco10 = df_eco_cooler[df_eco_cooler$Diving == "Sphenisciformes",]

df_eco1$Diving = 1
df_eco2$Diving = 2
df_eco3$Diving = 3
df_eco4$Diving = 4
df_eco5$Diving = 5
df_eco6$Diving = 6
df_eco7$Diving = 7
df_eco8$Diving = 8
df_eco9$Diving = 9
df_eco10$Diving = 10


df_eco_cooler1 = rbind(df_eco0, df_eco1, df_eco2, df_eco3, df_eco4, df_eco5, df_eco6, df_eco7, df_eco8, df_eco9, df_eco10)

df_eco_cooler = df_eco_cooler[,c(2,3)]
df_eco_cooler$Diving = 0
df_eco_cooler$Flight = as.numeric(as.character(df_eco_cooler$Flight))
df_eco_cooler$Diving = as.numeric(as.character(df_eco_cooler$Diving))
name.check(birds_tree, df_eco_cooler)

#PhyloPCA
birds_pca<-phyl.pca(birds_tree,df_eco_cooler)
birds_pca
par(mar=c(4.1,4.1,2.1,1.1),las=1) ## set margins
plot(birds_pca,main="")
birds_pca$Evec[,1]<--birds_pca$Evec[,1]
birds_pca$L[,1]<--birds_pca$L[,1]
birds_pca$S<-scores(birds_pca,
                       newdata=df_eco_cooler)

par(cex.axis=0.8,mar=c(5.1,5.1,1.1,1.1))
#everything works until here
phylomorphospace(birds_tree,
                 scores(birds_pca)[,1:2],
                 ftype="off",node.size=c(0,1),bty="n",las=1,
                 xlab="PC1 (overall size)",
                 ylab=expression(paste("PC2 ("%up%"lamellae number, "
                                       %down%"tail length)")))
eco<-setNames(ecomorph[,1],rownames(ecomorph))
ECO<-to.matrix(eco,levels(eco))
tiplabels(pie=ECO[ecomorph.tree$tip.label,],cex=0.5)
legend(x="bottomright",legend=levels(eco),cex=0.8,pch=21,
       pt.bg=rainbow(n=length(levels(eco))),pt.cex=1.5)