rm(list=ls(all=TRUE))

pdf("../../Body/4Figures/Humanz.ZaidiAnalysesR.pdf")

library(dplyr)

suppl = read.csv('../../Body/1Raw/Zaidi2019/pnas.1906331116.sd01.csv')
de_novo = read.csv('../../Body/1Raw/Zaidi2019/MakovaDeNovoMutations.csv')

names(de_novo)[2] = 'mother_id'
names(de_novo)[1] = 'family_id'
names(suppl)[2] = 'position'

mut = left_join(de_novo, suppl, by = c('family_id', 'position'))

table(mut$mutation)

# C>A C>G C>T T>A T>C T>G 
# 0   4   164   0 196   0 

mut$age_birth..yrs. = mut$age_birth..days. / 365

# hist(mut$age_birth..days.)
hist(mut$age_birth..yrs.)

# mut$age_birth..yrs. == Age of the mother at birth?

table(mut$Individual.id)

mut_fr = mut %>%
  group_by(Individual.id) %>%
  mutate(T_C_fr = sum(mutation == 'T>C') / length(mutation)
  ) %>%
  select(Individual.id, age_birth..days., age_birth..yrs., T_C_fr) %>%
  summarise(T_C = mean(T_C_fr), age = mean(age_birth..yrs.))

plot(mut_fr$T_C, log2(mut_fr$age))

cor.test(mut_fr$T_C, log2(mut_fr$age), method = 'spearman')

# data:  mut_fr$T_C and log2(mut_fr$age)
# S = 107756, p-value = 0.6305
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# -0.05291718 

mut_fr = mut %>%
  group_by(Individual.id) %>%
  mutate(T_C_fr = sum(ancestral_allele == 'T') / length(mutation)
  ) %>%
  select(Individual.id, age_birth..days., age_birth..yrs., T_C_fr) %>%
  summarise(T_C = mean(T_C_fr), age = mean(age_birth..yrs.))

plot(mut_fr$T_C, log2(mut_fr$age))

cor.test(mut_fr$T_C, log2(mut_fr$age), method = 'spearman')

# data:  mut_fr$T_C and log2(mut_fr$age)
# S = 99540, p-value = 0.8037
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.02736107 

dev.off()