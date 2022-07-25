rm(list = ls(all=TRUE))
#installing packages
install.packages('ape')
install.packages('geiger')
install.packages('caper')
#end

#checking packages
library(ape)
library(geiger)
library(caper)
#end

#reading tree, birds table, doing some shenanigans
df_all = read.csv('../../Head/2Scripts/Table_for_PGLS.csv')
tree = read.tree("../../Body/3Results/phylo.treefile")
phy=multi2di(tree)
row.names(df_all) = df_all$species_name
#end

#preparing files for PGLS
tree_all = treedata(phy, df_all, sort = T, warnings = T)$phy
all_pgls = as.data.frame(treedata(tree_all, df_all, sort = T, warnings = T)$data)
MutComp_all = comparative.data(tree_all, all_pgls, species_name, vcv = TRUE)
vec1 = names(all_pgls)
vec1 = vec1[2:58]
vec2 = names(all_pgls)
vec2 = vec2[2:58]
pgls_finale = data.frame()
#end

#PGLS loop
skip_to_next = FALSE
for (i in vec1)
{
  for (j in vec2) 
  {
    if (i != j)
    {
      skip_to_next = FALSE
      model_all = tryCatch(pgls(as.numeric(as.character(MutComp_all$data[,i])) ~ as.numeric(as.character(MutComp_all$data[,j])),  MutComp_all, lambda = 'ML'), error = function(e) { skip_to_next <<- TRUE})
      if (skip_to_next == FALSE)
      {
        p = as.numeric(summary(model_all)$coefficients[,4][2])
        r = as.numeric(summary(model_all)$r.squared)
        ef = as.numeric(summary(model_all)$coefficients[,1][2])
        pgls_finale = rbind(pgls_finale, c(i, j, p, r, ef))
      }
      if (skip_to_next == TRUE)
      {
        p = NA
        r = NA
        ef = NA
        pgls_finale = rbind(pgls_finale, c(i, j, p, r, ef))
      }
    }
  }
}
#end

#writing table we got
write.csv(pgls_finale, file = 'Birds_PGLS.csv')