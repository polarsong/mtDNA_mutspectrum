rm(list = ls(all=TRUE))
library(ggplot2)
library(heatmaply)

#reading table
df = read.csv('/../../Head/2Scripts/Birds_all_pgls.csv')

#discarding p-values
df_disc = df[is.na(df$p_value) == FALSE,]
df_disc = df_disc[df_disc$p_value <= 0.05,]
df_disc = df_disc[df_disc$p_value != 0,]

#adjusting p-values
df_disc$p_adj_bon = p.adjust(df_disc$p_value, method = 'bonferroni')
df_disc$p_adj_holm = p.adjust(df_disc$p_value, method = 'holm')
df_disc$p_adj_hochberg = p.adjust(df_disc$p_value, method = 'hochberg')
df_disc$p_adj_hommel = p.adjust(df_disc$p_value, method = 'hommel')
df_disc$p_adj_bh = p.adjust(df_disc$p_value, method = 'BH')
df_disc$p_adj_fdr = p.adjust(df_disc$p_value, method = 'fdr')
df_disc$p_adj_by = p.adjust(df_disc$p_value, method = 'BY')
df_disc$p_adj_none = p.adjust(df_disc$p_value, method = 'none')

ggplot(df_disc, aes(x = Ecology1, y = p_adj_bon))+
  geom_boxplot()


df_p = df_disc[c(1:3)]
row.names(df_p) = df_p$Ecology1
heatmaply(df_p)
