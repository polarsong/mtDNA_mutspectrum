# devtools::install_github('r-lib/systemfonts')
# BiocManager::install("svglite", dependencies = TRUE)


rm(list=ls(all=T))

library(ggtree)
library(ggplot2)
library(treeio)
library(dplyr)
library(ape)
library(geiger)
birds_tree = read.tree('anc_kg.treefile')
df_mtdna = read.csv('Birds_dataset_paper.csv')

df_need = data.frame()
for (i in unique(df_mtdna$species_name))
{
  a = df_mtdna[df_mtdna$species_name == i,]
  b = sum(a$ghahSkew)/12
  ab = c(i,b)
  df_need = rbind(df_need, ab)
}

names(df_need) = c('species_name', 'GhAhSkew')
df_need$species_name = gsub(' ', '_', df_need$species_name)



df_mutspec = read.table('C:/Users/User/Desktop/mutspec12.tsv', header = TRUE, fill = TRUE)
df_syn = df_mutspec[df_mutspec$Label == 'syn',]
df_syn = df_syn[,c(1,5,7)]
df_syn1 = reshape(data = df_syn, idvar = 'AltNode',
                  timevar = 'Mut',
                  direction = 'wide')
names(df_syn1) = c('species_name', 'Mutation_TG_syn', 'Mutation_TC_syn', 'Mutation_TA_syn', 'Mutation_GT_syn', 'Mutation_GC_syn', 
                   'Mutation_GA_syn', 'Mutation_CT_syn', 'Mutation_CG_syn', 'Mutation_CA_syn', 'Mutation_AT_syn', 'Mutation_AG_syn', 'Mutation_AC_syn')
df_syn1 = df_syn1[,c(1,12)]
df_syn2 = df_syn1[!grepl('Node', df_syn1$species_name),]
df_syn3 = df_syn1[grepl('Node', df_syn1$species_name),]

need_species = setdiff(df_need$species_name, df_syn2$species_name)
correct_need_species = setdiff(df_need$species_name, need_species)
birds_ms_tree<-keep.tip(birds_tree,correct_need_species)
df_syn1$Mutation_AG_syn = as.numeric(df_syn1$Mutation_AG_syn)
typeof(df_syn1$Mutation_AG_syn)


#plot_sbs_on_tree <- function(label, indir='./data/processed/sbs_on_tree/', outdir='./figures/') {
#  mutspec_file <- paste(indir, label, '_ff.tsv', sep="")
#  print(mutspec_file)
#  grad_val = read.table(mutspec_file, sep="\t", header = T)
#  names(grad_val) = c('RefNode', 'label', 'MutSpec')  # this label is not the label in input args
names(df_syn1) = c('label', 'MutSpec')
test = full_join(as_tibble(birds_tree), tibble(df_syn1[,c(1,2)]), by = 'label')
#test = na.omit(test)
der = as.treedata(test)
label = 'Mutspec'
pdf('Birds.pdf',         
    width = 40,height = 300)         
ggtree(der, aes(color = der@data$MutSpec)) + geom_tiplab(size = 2, colour = "darkgray") +
  scale_color_continuous(low="green", high="red")
dev.off()
