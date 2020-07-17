rm(list=ls(all=TRUE))

library(dplyr)

suppl = read.csv('../../Body/1Raw/Zaidi2019/pnas.1906331116.sd01.csv')
de_novo = read.csv('../../Body/1Raw/Zaidi2019/MakovaDeNovoMutations.csv')

names(de_novo)[2] = 'mother_id'
names(de_novo)[1] = 'family_id'
names(suppl)[2] = 'position'

mut = left_join(de_novo, suppl, by = c('family_id', 'position'))

setdiff(unique(mut$Individual.id), unique(de_novo$Individual.id))
# character(0)

table(mut$mutation)

# C>A C>G C>T T>A T>C T>G 
# 0   4   164   0 196   0 

mut$age_birth..yrs. = mut$age_birth..days. / 365

# hist(mut$age_birth..days.)
hist(mut$age_birth..yrs.)

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
  mutate(T_C_fr = sum(ancestral_allele == 'T') / length(mutation),
         notT_C = sum(ancestral_allele != 'T') / length(mutation),
         C_T_fr = sum(ancestral_allele == 'C') / length(mutation)
  ) %>%
  select(Individual.id, age_birth..days., age_birth..yrs., T_C_fr, notT_C, C_T_fr) %>%
  summarise(T_C = mean(T_C_fr), notT_C = mean(notT_C), C_T = mean(C_T_fr), 
            age = mean(age_birth..yrs.))

plot(mut_fr$T_C, log2(mut_fr$age))

cor.test(mut_fr$T_C, log2(mut_fr$age), method = 'spearman')

# data:  mut_fr$T_C and log2(mut_fr$age)
# S = 99540, p-value = 0.8037
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.02736107 


hist(mut_fr[mut_fr$T_C > 0,]$age, breaks = 20)
hist(mut_fr[mut_fr$notT_C > 0,]$age, breaks = 20)

summary(mut_fr[mut_fr$T_C > 0,]$age)
summary(mut_fr[mut_fr$notT_C > 0,]$age)
summary(mut_fr[mut_fr$C_T > 0,]$age)

# for each mutation in mut: position_anc_derived
# contrasts - presence of T>C

# mother id, kid id, number of de novo, number of T_C, C_T, etc, age of birth 

one_line = c()

for(i in unique(mut$mother_id.x)){
  # i = "F253m1"
  
  temp = mut[mut$mother_id.x == i,]
  
  for(j in unique(temp$Individual.id)){
    # j = 'F253m1c1'
    
    temp_kids = temp[temp$Individual.id == j,]
    numberMut = length(unique(temp_kids$position)) #number of de nove in table
    
    birthAge = unique(temp_kids$age_birth..days.) / 365 # mother' age of birth
    
    if(numberMut == 1){
      T_C_num = length(which(any(temp_kids$ancestral_allele == 'T')))
      G_A_num = length(which(any(temp_kids$ancestral_allele == 'G')))
      others_mut = length(which(any(temp_kids$ancestral_allele %in% c('A', 'C'))))
      C_T_num = length(which(any(temp_kids$ancestral_allele == 'C')))
      A_G_num = length(which(any(temp_kids$ancestral_allele == 'A')))
      
      one_line = rbind(one_line, c(i, j, numberMut, T_C_num, G_A_num, others_mut, 
                                   C_T_num, A_G_num, birthAge))
    }
    if(numberMut > 1){
      for(pos in unique(temp_kids$position)){
        temp_pos = temp_kids[temp_kids$position == pos,]
        T_C_num = length(which(any(temp_pos$ancestral_allele == 'T')))
        G_A_num = length(which(any(temp_pos$ancestral_allele == 'G')))
        others_mut = length(which(any(temp_pos$ancestral_allele %in% c('A', 'C'))))
        C_T_num = length(which(any(temp_pos$ancestral_allele == 'C')))
        A_G_num = length(which(any(temp_pos$ancestral_allele == 'A')))
        
        one_line = rbind(one_line, c(i, j, numberMut, T_C_num, G_A_num, others_mut, 
                                     C_T_num, A_G_num, birthAge))
      }
    }
  }
}

a = as.data.frame(one_line)
names(a) = c('mother_id', 'kid_id', 'mutNumber', 'T_C', 'G_A', 'othersMut', 'C_T',
             'A_G', 'motherAge')

to_merge = a 

data = a %>% 
  group_by(kid_id) %>%
  select(T_C:A_G) %>%
  mutate_all(as.character) %>%
  mutate_all(as.integer) %>%
  summarise_all(sum) 

final = inner_join(data, to_merge %>% select(mother_id, kid_id, mutNumber, motherAge) %>%
                     distinct)

##############################################################################
# siblings contrasts

final$motherAge = as.numeric(as.character(final$motherAge))
final$mutNumber = as.numeric(as.character(final$mutNumber))

one_line = c()
for(i in unique(final$mother_id)){
  # i = 'F163m1'
  temp = final[final$mother_id == i,]
  
  if(length(unique(temp$motherAge)) >= 2){
    older_kid = temp[temp$motherAge == max(temp$motherAge),]
    younger_kid = temp[temp$motherAge == min(temp$motherAge),]
    mutContrast = older_kid$mutNumber - younger_kid$mutNumber
    T_C_contrast = older_kid$T_C - younger_kid$T_C
    G_A_contrast = older_kid$G_A - younger_kid$G_A
    othContrastr = older_kid$othersMut - younger_kid$othersMut
    one_line = rbind(one_line, c(i, mutContrast, T_C_contrast, G_A_contrast, othContrastr))
  }
  if(length(unique(temp$motherAge)) > 2){print(i)}
}

contrastTable = as.data.frame(one_line)
names(contrastTable) = c('mother_id', 'mutNumber', 'T_C', 'G_A', 'others')

table(contrastTable$mutNumber)
table(contrastTable$T_C)
table(contrastTable$G_A)
table(contrastTable$others)


summary(final[final$T_C > 0,]$motherAge)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 18.65   26.52   31.41   30.15   34.74   38.74 

summary(final[final$G_A > 0,]$motherAge)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 21.05   25.89   30.92   30.28   33.79   38.74 

summary(final[final$C_T > 0,]$motherAge)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 18.65   29.92   31.30   29.97   33.28   37.50 

summary(final[final$A_G > 0,]$motherAge)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 22.46   25.93   29.72   29.55   32.87   39.08 

summary(final[final$G_A > 0 | final$C_T > 0 | final$A_G > 0,]$motherAge)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 18.65   26.78   30.61   30.04   33.54   39.08 

wilcox.test(final[final$G_A > 0 | final$C_T > 0 | final$A_G > 0,]$motherAge,
            final[final$T_C > 0,]$motherAge)


summary(final[final$othersMut > 0,]$motherAge)

nrow(final[final$T_C > 0,]) # 31
nrow(final[final$G_A > 0,]) # 32
nrow(final[final$C_T > 0,]) # 14
nrow(final[final$A_G > 0,]) # 27
nrow(final[final$G_A > 0 | final$C_T > 0 | final$A_G > 0,]) # 66
nrow(final[final$othersMut > 0,]) # 41

# wilcox.test(final[final$T_C > 0,]$motherAge, final[final$G_A > 0,]$motherAge)

