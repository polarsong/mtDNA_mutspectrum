rm(list = ls(all=TRUE))
library(ape); library(phytools);  library(geiger)
#root - Node494
birds_tree = read.tree('anc_kg.tre')
birds_tree$node.label <- NULL
is.rooted(birds_tree)
str(birds_tree)
plot(birds_tree)
birds_tree_correct = read.nexus('Correct_rooted_tree.txt')
is.rooted(birds_tree_correct)  
#feathertree <- root(feathertree, outgroup = "Struthio_camelus", resolve.root = TRUE)
#feathertree = FullTree
#mycalibration <- makeChronosCalib(feathertree, node="root", age.max=94)
#feathertree <- chronos(feathertree, lambda=0, model = 'correlated', calibration = mycalibration)  