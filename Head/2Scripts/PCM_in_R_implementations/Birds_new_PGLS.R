#new PGLS for birds
rm(list = ls(all=TRUE))
library(phytools)
library(nlme)
library(geiger)
#flying birds
df_flight = read.csv('flight_and_gene.csv')
df_flight = df_flight[,c(2,3,4,5,6)]
df_flight$flight_num = 0
df_flight[df_flight$ability_to_fly != 'Flying',]$flight_num = 1
rownames(df_flight) = df_flight$species_name
flying_tree = read.tree('flying_birds_tree.tre')
name.check(flying_tree, df_flight)
df_flight$ability_to_fly = as.character(df_flight$ability_to_fly)
df_flight$GhAhSkew = as.numeric(as.character(df_flight$GhAhSkew))
df_flight$fAn = as.numeric(as.character(df_flight$fAn))
df_flight$fGn = as.numeric(as.character(df_flight$fGn))
df_flight$flight_num = as.numeric(as.character(df_flight$flight_num))
spp = rownames(df_flight)
corBM = corBrownian(phy=flying_tree,form=~spp)
corBM
pgls_flying = gls(GhAhSkew~ability_to_fly,
                  data=df_flight,correlation=corBM)
summary(pgls_flying)

anova(pgls_flying)

pgls_flying_1 = gls(GhAhSkew~flight_num,
                  data=df_flight,correlation=corBM)
summary(pgls_flying_1)


anova(pgls_flying_1)

pgls_flying_2 = gls(fGn~ability_to_fly,
                    data=df_flight,correlation=corBM)
summary(pgls_flying_2)


anova(pgls_flying_2)

pgls_flying_3 = gls(fGn~flight_num,
                    data=df_flight,correlation=corBM)
summary(pgls_flying_3)


anova(pgls_flying_3)


pgls_flying_4 = gls(fAn~ability_to_fly,
                    data=df_flight,correlation=corBM)
summary(pgls_flying_4)


anova(pgls_flying_4)

pgls_flying_4 = gls(fAn~flight_num,
                    data=df_flight,correlation=corBM)
summary(pgls_flying_4)


anova(pgls_flying_4)


#diving birds
df_dive = read.csv('dive_and_gene.csv')
df_dive = df_dive[,c(2,3,4,5,6)]
df_dive$dive_num = 0
df_dive[df_dive$ability_to_dive != 'Non-diving',]$dive_num = 1
rownames(df_dive) = df_dive$species_name
diving_tree = read.tree('diving_birds_tree.tre')
name.check(diving_tree, df_dive)
df_dive$ability_to_dive = as.character(df_dive$ability_to_dive)
df_dive$GhAhSkew = as.numeric(as.character(df_dive$GhAhSkew))
df_dive$fAn = as.numeric(as.character(df_dive$fAn))
df_dive$fGn = as.numeric(as.character(df_dive$fGn))
df_dive$dive_num = as.numeric(as.character(df_dive$dive_num))
spp = rownames(df_dive)
corBM = corBrownian(phy=diving_tree,form=~spp)
corBM
pgls_dive = gls(GhAhSkew~ability_to_dive,
                  data=df_dive,correlation=corBM)
summary(pgls_dive)

anova(pgls_dive)

pgls_dive_1 = gls(GhAhSkew~dive_num,
                    data=df_dive,correlation=corBM)
summary(pgls_dive_1)


anova(pgls_dive_1)

pgls_dive_2 = gls(fGn~ability_to_dive,
                    data=df_dive,correlation=corBM)
summary(pgls_dive_2)


anova(pgls_dive_2)

pgls_dive_3 = gls(fGn~dive_num,
                    data=df_dive,correlation=corBM)
summary(pgls_dive_3)


anova(pgls_dive_3)


pgls_dive_4 = gls(fAn~ability_to_dive,
                    data=df_dive,correlation=corBM)
summary(pgls_dive_4)


anova(pgls_dive_4)

pgls_dive_4 = gls(fAn~dive_num,
                    data=df_dive,correlation=corBM)
summary(pgls_dive_4)


anova(pgls_dive_4)