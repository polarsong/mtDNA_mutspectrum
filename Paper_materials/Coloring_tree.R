# devtools::install_github('r-lib/systemfonts')
# BiocManager::install("svglite", dependencies = TRUE)
#mut_data = read.table("C:/Users/User/Desktop/Birds mutspec results from Bogdan/mutspec12.tsv", header = TRUE, fill = TRUE)
#color_file = mut_data[,c(1,5,7)]
#color_file = color_file[color_file$Mut == 'T>C',]
#color_file$Mut = 'A>G'
#write.csv(color_file, 'color_file.csv', row.names = FALSE)
#rm(list=ls(all=T))
#setwd("birds")
rm(list = ls(all=TRUE))
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")

library(ggtree)
library(ggplot2)
library(treeio)
library(dplyr)
library(ape)

# function_name <- function(arg_1, arg_2, ...) {
#   Function body 
# }


my_tree = read.tree("../paper_materials/anc_kg.treefile")

plot_sbs_on_tree <- function(label, indir='../paper_materials', outdir='../paper_materials') {
  grad_val = read.table('../paper_materials/color_file.csv', sep=",", header = T)
  names(grad_val) = c('label', 'MutSpec', 'RefNode')  # this label is not the label in input args
  
  test = full_join(as_tibble(my_tree), tibble(grad_val[,c(1,2)]), by = 'label')
  der = as.treedata(test)
  
  ext <- 'pdf'
  out_image_file <- paste(outdir, label, 'bird_colored_tree', ext, sep = '')
  der <- groupOTU(der, c("Aptenodytes_forsteri", "Pygoscelis_adeliae"))
  ggtree(der, aes(color = der@data$MutSpec)) + geom_tiplab(size = 3, colour = "darkgray") +
    scale_color_continuous(low="green", high="red") +
    theme(legend.position="bottom") +
    labs(col=label) +
    geom_tippoint(aes(alpha = group), col = "red") +
    scale_color_manual(values = c(0,1), aesthetics = "alpha")
  
  ggsave(out_image_file,  width = 2000, height = 2000, units = 'px', scale=2.2)
}
label <- 'A>G'
plot_sbs_on_tree(label)

paste
for (n1 in c('A', 'C', 'G', 'T')) {
  for (n2 in c('A', 'C', 'G', 'T')) {
    if (n1 == n2) {
      next
    }
    label <- paste(n1, n2, sep = ">")
    print(label)
    plot_sbs_on_tree(label)
  }
}

# A>T / T>A
label = 'A>T_T>A_ratio'
plot_sbs_on_tree(label)




# T on tree
label <- 'T share'
mutspec_file <- './data/processed/percent_of_T.tsv'
# out_image_file <- './figures/T_ff_tree.svg'
out_image_file <- './figures/T_tree.png'


grad_val = read.table(mutspec_file, sep="\t", header = T)
names(grad_val) = c('RefNode', 'label', 'MutSpec')

test = full_join(as_tibble(my_tree), tibble(grad_val[,c(2,3)]), by = 'label')
der = as.treedata(test)

# svg(out_image_file, width = 20, height = 20)
png(out_image_file, width = 2000, height = 2000)
der <- groupOTU(der, c("Halicephalobus_mephisto", "Caenorhabditis_elegans"))
ggtree(der, aes(color = der@data$MutSpec)) + geom_tiplab(size = 4, colour = "darkgray") +
  scale_color_continuous(low="green", high="red") +
  theme(legend.position="bottom")+
  labs(col=label)+
  geom_tippoint(aes(alpha = group), col = "purple")+
  scale_color_manual(values = c(0,1), aesthetics = "alpha")

dev.off()