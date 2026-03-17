rm(list = ls(all=TRUE))
df_nd6 = read.csv('../Birds_mtDNA_data.csv')
df_nd6$GhAhSkew = (df_nd6$neutral_c- df_nd6$neutral_T)/(df_nd6$neutral_c + df_nd6$neutral_T)
df_nd6$ThChSkew = (df_nd6$neutral_A - df_nd6$neutral_g)/(df_nd6$neutral_A + df_nd6$neutral_g)

df_need = data.frame()
for (i in unique(df_nd6$species_name))
{
  a = df_nd6[df_nd6$species_name == i,]
  b = sum(a$GhAhSkew)/12
  ab = c(i,b)
  df_need = rbind(df_need, ab)
}
names(df_need) = c('species_name', 'GhAhSkew')
df_fly = read.csv('../flying_birds.csv')
df_fly1 = df_fly[df_fly$Flightless != 'Almost_flightless',]
df_fly1 = df_fly1[df_fly1$Flightless != 'Galliformes',]
df_fly1 = df_fly1[df_fly1$Flightless != 'Flightless',]
#merge data, grab tree and go
need_species = setdiff(df_flight2$species_name, df_flight4$species_name)
correct_need_species = setdiff(df_flight2$species_name, need_species)
midori_data_tree<-keep.tip(birds_ms_and_temp_tree,correct_need_species)
df_avonet = read.csv('../../../Body/1Raw/Avonet_data.csv')
df_birds_names = df_avonet[,c(1,2)]
names(df_birds_names) = c('species_name', 'x')
df_birds_names$species_name = gsub(' ', '_', df_birds_names$species_name)
df_birds_names1 = merge(df_birds_names, df_need1)

row.names(df_flight4) = df_flight4$species_name
df_flight4$Midori_AG_mutspec = as.numeric(as.character(df_flight4$Midori_AG_mutspec))
df_flight4$GhAhSkew = as.numeric(as.character(df_flight4$GhAhSkew))
lnTL<-setNames(df_flight4$GhAhSkew,rownames(df_flight4))
head(lnTL)
## estimate ancestral states using fastAnc
fit.lnTL<-fastAnc(midori_data_tree,lnTL,vars=TRUE,CI=TRUE)
print(fit.lnTL,printlen=10)
## compute "contMap" object
birds_contMap<-contMap(midori_data_tree,lnTL,
                       plot=FALSE)
## plot "contMap" object
plot(birds_contMap,sig=2,fsize=c(0.45,0.9),
     lwd=c(2,3))
tips<-extract.clade(midori_data_tree,'Node690')$tip.label #699 - peng, 690 peng + ant 582 non-flying 496 ducks
tips
## prune "contMap" object to retain only these tips
pruned.contMap<-keep.tip.contMap(birds_contMap,tips)
## plot object
plot(pruned.contMap)